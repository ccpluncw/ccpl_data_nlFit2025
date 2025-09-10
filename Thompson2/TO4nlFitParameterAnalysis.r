library(chutils)
library(dplyr)
library(tidyr)
library(effectsize)

mainDir <- getwd()
ch.newDir (mainDir, "TO4_SummaryStats")

datafile <- "./Feb_13_2025_13-57_TO4_ParameterData.txt"
df.data <-read.table(datafile, header=T, sep="\t", quote="\"")

df.data <- df.data[df.data$sn != "all", ]
uniqueSNs <- length(unique(df.data$sn))

df.data$r2 <- ifelse(df.data$r2 < 0, 0, df.data$r2)
df.data$r2.ll <- ifelse(df.data$r2.ll < 0, 0, df.data$r2.ll)
df.data$r2.ll.1 <- ifelse(df.data$r2.ll.1 < 0, 0, df.data$r2.ll.1)

t.test.r2.nlVll <- with(df.data, t.test(r2, r2.ll, paired = T))
t.test.r2.nlVll.1 <- with(df.data, t.test(r2 , r2.ll.1, paired = T))

t.test.BIC.nlVll <- with(df.data, t.test(BIC , BIC.ll, paired = T))
t.test.BIC.nlVll.1 <-with(df.data, t.test(BIC , BIC.ll.1, paired = T))

cohensD.r2.nlVll <- with(df.data, cohens_d(r2, r2.ll, paired = T))
cohensD.r2.nlVll.1 <- with(df.data, cohens_d(r2 , r2.ll.1, paired = T))

cohensD.BIC.nlVll <- with(df.data, cohens_d(BIC , BIC.ll, paired = T))
cohensD.BIC.nlVll.1 <-with(df.data, cohens_d(BIC , BIC.ll.1, paired = T))


df.sum <- data.frame(df.data %>% group_by(gd, tm) %>% summarize(mR2 = mean(r2, na.rm = T), sdR2 = sd(r2, na.rm = T), mR2.ll = mean(r2.ll, na.rm = T), sdR2.ll = sd(r2.ll, na.rm = T), mR2.cll = mean(r2.ll.1, na.rm = T), sdR2.cll = sd(r2.ll.1, na.rm = T), mAcc = mean(accuracyPercent, na.rm = T), sdAcc = sd(accuracyPercent, na.rm = T), mNS = mean(numberSensitivity, na.rm = T), sdNS = sd(numberSensitivity, na.rm = T), mPI = mean(pIncludeConceptualPoints, na.rm = T), sdPI = sd(pIncludeConceptualPoints, na.rm = T)))

df.sum.gd <- data.frame(df.data %>% group_by(gd) %>% summarize(mR2 = mean(r2, na.rm = T), sdR2 = sd(r2, na.rm = T), mR2.ll = mean(r2.ll, na.rm = T), sdR2.ll = sd(r2.ll, na.rm = T), mR2.cll = mean(r2.ll.1, na.rm = T), sdR2.cll = sd(r2.ll.1, na.rm = T), mAcc = mean(accuracyPercent, na.rm = T), sdAcc = sd(accuracyPercent, na.rm = T), mNS = mean(numberSensitivity, na.rm = T), sdNS = sd(numberSensitivity, na.rm = T), mPI = mean(pIncludeConceptualPoints, na.rm = T), sdPI = sd(pIncludeConceptualPoints, na.rm = T)))

df.sum.tm <- data.frame(df.data %>% group_by(tm) %>% summarize(mR2 = mean(r2, na.rm = T), sdR2 = sd(r2, na.rm = T), mR2.ll = mean(r2.ll, na.rm = T), sdR2.ll = sd(r2.ll, na.rm = T), mR2.cll = mean(r2.ll.1, na.rm = T), sdR2.cll = sd(r2.ll.1, na.rm = T), mAcc = mean(accuracyPercent, na.rm = T), sdAcc = sd(accuracyPercent, na.rm = T), mNS = mean(numberSensitivity, na.rm = T), sdNS = sd(numberSensitivity, na.rm = T), mPI = mean(pIncludeConceptualPoints, na.rm = T), sdPI = sd(pIncludeConceptualPoints, na.rm = T)))

df.sum.all <- data.frame(df.data %>% summarize(mR2 = mean(r2, na.rm = T), sdR2 = sd(r2, na.rm = T), mR2.ll = mean(r2.ll, na.rm = T), sdR2.ll = sd(r2.ll, na.rm = T), mR2.cll = mean(r2.ll.1, na.rm = T), sdR2.cll = sd(r2.ll.1, na.rm = T), mAcc = mean(accuracyPercent, na.rm = T), sdAcc = sd(accuracyPercent, na.rm = T), mNS = mean(numberSensitivity, na.rm = T), sdNS = sd(numberSensitivity, na.rm = T), mPI = mean(pIncludeConceptualPoints, na.rm = T), sdPI = sd(pIncludeConceptualPoints, na.rm = T)))

df.sum.gd$gd <- factor(df.sum.gd$gd , levels=c("G2", "G3", "G6", "adult"))
df.sum.tm$tm <- factor(df.sum.tm$tm , levels=c("1", "10", "100"))


df.data$gd <- factor(df.data$gd , levels=c("G2", "G3", "G6", "adult"))
df.data$tm <- factor(df.data$tm , levels=c("1", "10", "100"))
contrasts(df.data$tm, how.many=1) <- contr.poly(3)
contrasts(df.data$gd, how.many=1) <- contr.poly(4)
ap.aov <- aov(accuracyPercent ~ gd*tm, data = df.data)
eta.ap <- eta_squared(ap.aov)

ns.aov <- aov(numberSensitivity ~ gd*tm, data = df.data)
eta.ns <- eta_squared(ns.aov)

sink("TO4_ParameterAnalysis.txt")
	cat("\nNumber of Unique Subjects:", uniqueSNs, "\n\n")
	cat("\n********** summary stats for nlFit **********\n\n")
	print(df.sum.all)
	print(df.sum.gd)
	print(df.sum.tm)
	print(df.sum)

	cat("\n********** Accuracy Percent: ANOVA with linear contrasts Grade and Range **********\n\n")
	print(anova(ap.aov))
	print(eta.ap)
	cat("\n********** Number Sensitivity: ANOVA with linear contrasts Grade and Range **********\n\n")
	print(anova(ns.aov))
	print(eta.ns)
	cat("\n********** t-test comparing r2 for nlFit vs LL **********\n\n")
	print(t.test.r2.nlVll)
	print(cohensD.r2.nlVll)
	cat("\n\n********** t-test comparing r2 for nlFit vs constrained LL **********\n\n")
	print(t.test.r2.nlVll.1)
	print(cohensD.r2.nlVll.1)
	cat("\n********** t-test comparing BIC for nlFit vs LL **********\n\n")
	print(t.test.BIC.nlVll)
	print(cohensD.BIC.nlVll)
	cat("\n\n********** t-test comparing BIC for nlFit vs constrained LL **********\n\n")
	print(t.test.BIC.nlVll.1)
	print(cohensD.BIC.nlVll.1)
sink(NULL)

pdf("TO4_ParameterFigures.pdf")
	with(df.data, plot(r2 ~ r2.ll, xlim = c(0,1), ylim = c(0,1), frame = F))
	abline(0,1)
	with(df.data, plot(r2 ~ r2.ll.1, xlim = c(0,1), ylim = c(0,1), frame = F))
	abline(0,1)

	with(df.data, plot(BIC ~ BIC.ll, frame = F))
	abline(0,1)
	with(df.data, plot(BIC ~ BIC.ll.1, frame = F))
	abline(0,1)

	with(df.data, boxplot(r2 , r2.ll, ylim = c(0,1), names = c("nlFit", "LL"), ylab = "r2", frame = F))
	with(df.data, boxplot(r2 , r2.ll.1, ylim = c(0,1), names = c("nlFit", "Constrained LL"), ylab = "r2", frame = F))
	with(df.data, boxplot(BIC , BIC.ll, names = c("nlFit", "LL"), ylab = "BIC", frame = F))
	with(df.data, boxplot(BIC , BIC.ll.1, names = c("nlFit", "Constrained LL"), ylab = "BIC", frame = F))

	with(df.data, plot(accuracyPercent ~ numberSensitivity, xlim = c(0,1), ylim = c(0,1), frame = F))
	with(df.data, plot(accuracyPercent ~ pIncludeConceptualPoints, xlim = c(0,1), ylim = c(0,1), frame = F))

	with(df.data, plot(accuracyPercent ~ coef.lam.ll), frame = F)
	with(df.data, plot(accuracyPercent ~ coef.int.ll, frame = F))

	with(df.data, plot(numberSensitivity ~ coef.lam.ll, frame = F))
	with(df.data, plot(numberSensitivity ~ coef.int.ll, frame = F))

	df.data$gd <- factor(df.data$gd , levels=c("G2", "G3", "G6", "adult"))
	
	with(df.data, boxplot(accuracyPercent ~ tm, ylim = c(0,1), frame = F))
	with(df.data, boxplot(accuracyPercent ~ gd, ylim = c(0,1), frame = F))
	
	tms <- unique(df.data$tm) 
	for(tmv in tms) {
		with(df.data[df.data$tm==tmv, ], boxplot(accuracyPercent ~ gd, ylim = c(0,1), frame = F, main = tmv))
	}
	with(df.data, boxplot(accuracyPercent ~ gd*tm, ylim = c(0,1), frame = F))

	with(df.data, boxplot(numberSensitivity ~ tm, ylim = c(0,1), frame = F))
	with(df.data, boxplot(numberSensitivity ~ gd, ylim = c(0,1), frame = F))

	for(tmv in tms) {
		with(df.data[df.data$tm==tmv, ], boxplot(numberSensitivity ~ gd, ylim = c(0,1), frame = F, main = tmv))
	}
	with(df.data, boxplot(numberSensitivity ~ gd*tm, ylim = c(0,1), frame = F))

	with(df.sum.gd, barplot(mAcc ~ gd, ylim = c(0,1)))
	with(df.sum.gd, barplot(mNS ~ gd, ylim = c(0,1)))
	with(df.sum.gd, barplot(mR2 ~ gd, ylim = c(0,1)))
	with(df.sum.tm, barplot(mAcc ~ tm, ylim = c(0,1)))
	with(df.sum.tm, barplot(mNS ~ tm, ylim = c(0,1)))
	with(df.sum.tm, barplot(mR2 ~ tm, ylim = c(0,1)))
	
dev.off()
setwd(mainDir)


		