library(chutils)
library(dplyr)
library(tidyr)
library(readxl)

mainDir <- getwd()
nlFitDatafile.a <- "./TO3.a_final.v12/Jan_23_2025_00-23_NLfit_snFits.txt"
nlFitDatafile.b <- "./TO3.b_final.v12/Jan_23_2025_00-23_NLfit_snFits.txt"
df.nlFit.a <-read.table(nlFitDatafile.a, header=T, sep="\t", quote="\"", check.names = FALSE)
df.nlFit.b <-read.table(nlFitDatafile.b, header=T, sep="\t", quote="\"", check.names = FALSE)

df.nlFit <- rbind(df.nlFit.a, df.nlFit.b)
#remove previous LL model
df.nlFit$r2.ll <- NULL
df.nlFit$BIC.ll <- NULL

rawDatafile <- "TO3data.xlsx"
sheets <- excel_sheets(rawDatafile)
data.long <- NULL
df.data.all <- list()
for(i in 1:length(sheets)) {
  #read in the clinicians billed units
  df.data.all[[sheets[i]]] <- data.frame(read_excel(rawDatafile, sheet = sheets[i], .name_repair = "minimal"), check.names = FALSE)
	
	columnNames <- colnames(df.data.all[[sheets[i]]])	
	targetColumns <- columnNames[grepl("^[[:digit:]]+", columnNames)]
	data.long.tmp <- df.data.all[[sheets[i]]] %>% pivot_longer (cols = all_of(targetColumns), names_to = "target", values_to = "estimate")
	data.long.tmp$st <- sheets[i]
	data.long <- ch.rbind(data.long, data.long.tmp)
}
data.long$target <- as.numeric(as.character(data.long$target))

lowerBound <- 0
upperBound <- 1000

df.out <- NULL
for(idx in 1:nrow(df.nlFit)) {
	st <- df.nlFit$st[idx]
	pc <- df.nlFit$pc[idx]
	nm <- df.nlFit$nm[idx]
	sn <- df.nlFit$sn[idx]
	age <- df.nlFit$age[idx]

	if(sn == "all") {
		df.tmp <- data.long[data.long$st == st & data.long$proc == pc & data.long$num == nm, ] %>% group_by(target) %>% summarize(estimate = mean(estimate, na.rm = T), Age = mean(Age, na.rm = t))
	} else {
		df.tmp <- data.long[data.long$st == st & data.long$proc == pc & data.long$num == nm & data.long$sn == sn, ]		
	}
	
	#add constrained LL model
	df.tmp$logTarget <- (upperBound/log(upperBound)) * log(df.tmp$target) 
	####### LL without the a
	fit.ll.1 <- tryCatch ({
		 nls(estimate ~ ( ((1-lambda)*target) + (lambda*logTarget)), data=df.tmp, start=list(lambda = 1), 					algorithm="port", lower=c(0), upper=c(1), control = nls.control(minFactor=1/10000000, maxiter=10000, warnOnly = FALSE))
	}, error = function(x) {
	            print(paste("Constrained LL Model Failed on ", sn, st, pc, nm))								
	            # Choose a return value in case of error
	            NULL
	})
	
#	fit.ll.1 <- nls(estimate ~ ( ((1-lambda)*target) + (lambda*logTarget)), data=df.tmp, start=list(lambda = 1), 					algorithm="port", lower=c(0), upper=c(1), control = nls.control(minFactor=1/10000000, maxiter=10000, warnOnly = FALSE))
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
#	fit.ll <- nls(estimate ~ a*( ((1-lambda)*target) + (lambda*logTarget)), data=df.tmp, start=list(a=1, lambda = 1), 					algorithm="port", lower=c(0), upper=c(1), control = nls.control(minFactor=1/10000000, maxiter=10000, warnOnly = FALSE))
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
ch.newDir (nlFitDir, "TO3_SummaryStats")
fileTag <- paste(format(Sys.time(), "%b_%d_%Y_%H-%M"), sep="_")
fileName <- paste(fileTag, "TO_OT_ParameterData.txt", sep="_")
write.table(df.out, fileName, append = FALSE, col.names=T, row.names=F, quote=F, sep="\t")

setwd(mainDir)

