library(chutils)
library(dplyr)
library(tidyr)

screenWidth <- 1920
rightMargin <- 100
unitErrThreshold <- 3

mainDir <- getwd()
nlFitDatafile.a <- "./unbounded_0_final.v12/Jan_21_2025_16-28_ONL_snFits.txt"
nlFitDatafile.b <- "./bounded_0_final.v12/Jan_21_2025_16-28_ONL_snFits.txt"
df.nlFit.a <-read.table(nlFitDatafile.a, header=T, sep="\t", quote="\"", check.names = FALSE)
df.nlFit.b <-read.table(nlFitDatafile.b, header=T, sep="\t", quote="\"", check.names = FALSE)

df.nlFit <- rbind(df.nlFit.a, df.nlFit.b)
#remove previous LL model
df.nlFit$r2.ll <- NULL
df.nlFit$BIC.ll <- NULL

rawDatafile <- "alldataEP.txt"
#read in the sheets of the clinician production excel file
df.in <- read.table(rawDatafile, header = T, sep="\t")
#remove Practice
df.in.1 <- df.in[df.in$pract > 0, ]
#filter high unitErr
df.in.2 <- df.in.1[df.in.1$unitErr < unitErrThreshold, ]
#get only important columns
relevantColumns <- c("sn","trial","session","leftMargin", "width", "bUnitValue","tUnitValue", "static", "userRespValue", "estTask", "estStimTime")
df.data <-df.in.2[,c(relevantColumns)]
names(df.data)[names(df.data) == 'tUnitValue'] <- 'target'
names(df.data)[names(df.data) == 'userRespValue'] <- 'estimate'

df.data$cond <- ifelse(df.data$static == TRUE, paste("bounded", df.data$estStimTime, sep = "_"), paste("unbounded", df.data$estStimTime, sep = "_"))
df.data$lowerBound <- 0
df.data$unitPixels <- df.data$width/df.data$bUnitValue
df.data$upperBound <- ifelse(df.data$static == TRUE, df.data$bUnitValue, (screenWidth - rightMargin - df.data$leftMargin) / df.data$unitPixels )
#start with production task:
df.data <- df.data[df.data$estTask == FALSE, ]

data.long <- df.data
df.out <- NULL
for(idx in 1:nrow(df.nlFit)) {
	cnd <- df.nlFit$cnd[idx]
	sn <- df.nlFit$sn[idx]
	age <- "adult"

	if(sn == "all") {
		df.tmp.1 <- data.long[data.long$cond == cnd, ]
	} else {
		df.tmp.1 <- data.long[data.long$cond == cnd & data.long$sn == sn, ] 
	}

	upperBound <- mean(unique(df.tmp.1$upperBound), na.rm=T)
	lowerBound <- 0
	df.tmp <- df.tmp.1 %>% group_by(target) %>% summarize(estimate = mean(estimate, na.rm = T))
	
	df.tmp$logTarget <- (22/log(22)) * log(df.tmp$target)
	
	####### LL without the a
	fit.ll.1 <- tryCatch ({
		 nls(estimate ~ ( ((1-lambda)*target) + (lambda*logTarget)), data=df.tmp, start=list(lambda = 1), 					algorithm="port", lower=c(0), upper=c(1), control = nls.control(minFactor=1/10000000, maxiter=10000, warnOnly = FALSE))
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
ch.newDir (nlFitDir, "CO_SummaryStats")
fileTag <- paste(format(Sys.time(), "%b_%d_%Y_%H-%M"), sep="_")
fileName <- paste(fileTag, "CO_ParameterData.txt", sep="_")
write.table(df.out, fileName, append = FALSE, col.names=T, row.names=F, quote=F, sep="\t")

setwd(mainDir)

