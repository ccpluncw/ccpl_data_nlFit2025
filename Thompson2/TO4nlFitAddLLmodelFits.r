library(chutils)
library(dplyr)
library(tidyr)
library(readxl)

mainDir <- getwd()
nlFitDatafile.a <- "./TO4_final.v12/Feb_12_2025_22-27_ONL_snFits.txt"
df.nlFit <-read.table(nlFitDatafile.a, header=T, sep="\t", quote="\"", check.names = FALSE)

#remove previous LL model
df.nlFit$r2.ll <- NULL
df.nlFit$BIC.ll <- NULL

rawDatafile <- "TO4data.xlsx"
sheets <- excel_sheets(rawDatafile)

df.data.all <- list()
df.data.all <- data.frame(read_excel(rawDatafile, sheet = sheets[1], .name_repair = "minimal"), check.names = FALSE)
df.data.all$st <-  sheets[1]

data.in <- df.data.all
columnNames <- colnames(data.in)

targetColumns <- columnNames[grepl("^[[:digit:]]+", columnNames)]

data.long <- data.in %>% pivot_longer (cols = all_of(targetColumns), names_to = "target", values_to = "estimate")
data.long$target <- as.numeric(as.character(data.long$target))

df.out <- NULL
for(idx in 1:nrow(df.nlFit)) {
	st <- df.nlFit$st[idx]
	tm <- df.nlFit$tm[idx]
	gd <- df.nlFit$gd[idx]
	sn <- df.nlFit$sn[idx]

	lowerBound <- 0
	upperBound <- 1000 * tm

	if(sn == "all") {
		df.tmp <- data.long[data.long$st == st & data.long$targetMultiplier == tm & data.long$grade == gd, ] %>% group_by(target) %>% summarize(estimate = mean(estimate, na.rm = T))
	} else {
		df.tmp <- data.long[data.long$st == st & data.long$targetMultiplier == tm & data.long$grade == gd & data.long$sn == sn, ]		
	}
	
	#add constrained LL model
	df.tmp$logTarget <- (upperBound/log(upperBound)) * log(df.tmp$target) 
	####### LL without the a
	fit.ll.1 <- tryCatch ({
		 nls(estimate ~ ( ((1-lambda)*target) + (lambda*logTarget)), data=df.tmp, start=list(lambda = 1), algorithm="port", lower=c(0), upper=c(1), control = nls.control(minFactor=1/10000000, maxiter=10000, warnOnly = FALSE))
	}, error = function(x) {
	            print(paste("Constrained LL Model Failed on ", sn, st, pc, nm))								
	            # Choose a return value in case of error
	            NULL
	})
	
	if(!is.null(fit.ll.1)) {
		df.tmp$fit.ll.1 <- predict(fit.ll.1, newdata = df.tmp)
		coef.lam.ll.1 <- round(coef(fit.ll.1)[[1]], 2)
		r2.ll.1 <- round(ch.R2(df.tmp$estimate, df.tmp$fit.ll.1, standardize = FALSE),2)
		BIC.ll.1 <- round(ch.BIC(df.tmp$estimate, df.tmp$fit.ll.1, 1, standardize = FALSE), 0)
	} else {
		coef.lam.ll.1 <- NA
		r2.ll.1 <- NA
		BIC.ll.1 <- NA
	}

	### LL model with the a
	####### LL without the a
	fit.ll <- tryCatch ({
		 nls(estimate ~ a*( ((1-lambda)*target) + (lambda*logTarget)), data=df.tmp, start=list(a=1, lambda = 1), 					algorithm="port", lower=c(-Inf, 0), upper=c(Inf, 1), control = nls.control(minFactor=1/10000000, maxiter=10000, warnOnly = FALSE))
	}, error = function(x) {
	            print(paste ("Full LL Model Failed on ", sn, st, pc, nm))								
	            # Choose a return value in case of error
	            NULL
	})

	if(!is.null(fit.ll)) {
		df.tmp$fit.ll <- predict(fit.ll, newdata = df.tmp)
		coef.int.ll <- round(coef(fit.ll)[[1]],2)
		coef.lam.ll <- round(coef(fit.ll)[[2]], 2)
		r2.ll <- round(ch.R2(df.tmp$estimate, df.tmp$fit.ll, standardize = FALSE),2)
		BIC.ll <- round(ch.BIC(df.tmp$estimate, df.tmp$fit.ll, 2, standardize = FALSE), 0)
	} else {
		coef.int.ll <- NA
		coef.lam.ll <- NA
		r2.ll <- NA
		BIC.ll <- NA
		
	}

	df.out.tmp <- cbind(df.nlFit[idx, ], data.frame(coef.int.ll, coef.lam.ll, r2.ll,BIC.ll, coef.lam.ll.1, r2.ll.1, BIC.ll.1))
	
	df.out <- ch.rbind(df.out, df.out.tmp)
	
}

nlFitDir <- getwd()
ch.newDir (nlFitDir, "TO4_SummaryStats")
fileTag <- paste(format(Sys.time(), "%b_%d_%Y_%H-%M"), sep="_")
fileName <- paste(fileTag, "TO4_ParameterData.txt", sep="_")
write.table(df.out, fileName, append = FALSE, col.names=T, row.names=F, quote=F, sep="\t")

setwd(mainDir)

