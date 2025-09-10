library(chutils)
library(dplyr)
library(tidyr)
library(effectsize)

mainDir <- getwd()
ch.newDir (mainDir, "TO3_SummaryStats")

datafile <- "./Feb_13_2025_13-56_TO_OT_ParameterData.txt"
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

t.test.ap.pc <- with(df.data, t.test(accuracyPercent ~ pc, paired = F))
t.test.ap.nm <- with(df.data, t.test(accuracyPercent ~ nm, paired = F))
t.test.ns.pc <- with(df.data, t.test(numberSensitivity ~ pc, paired = F))
t.test.ns.nm <- with(df.data, t.test(numberSensitivity ~ nm, paired = F))

cohensD.ap.pc <- with(df.data, cohens_d(accuracyPercent ~ pc, paired = F))
cohensD.ap.nm <- with(df.data, cohens_d(accuracyPercent ~ nm, paired = F))
cohensD.ns.pc <- with(df.data, cohens_d(numberSensitivity ~ pc, paired = F))
cohensD.ns.nm <- with(df.data, cohens_d(numberSensitivity ~ nm, paired = F))

df.sum.all <- data.frame(df.data %>% summarize(mR2 = mean(r2, na.rm = T), sdR2 = sd(r2, na.rm = T), mR2.ll = mean(r2.ll, na.rm = T), sdR2.ll = sd(r2.ll, na.rm = T), mR2.cll = mean(r2.ll.1, na.rm = T), sdR2.cll = sd(r2.ll.1, na.rm = T), mAcc = mean(accuracyPercent, na.rm = T), sdAcc = sd(accuracyPercent, na.rm = T), mNS = mean(numberSensitivity, na.rm = T), sdNS = sd(numberSensitivity, na.rm = T), mPI = mean(pIncludeConceptualPoints, na.rm = T), sdPI = sd(pIncludeConceptualPoints, na.rm = T)))

df.sum.int <- data.frame(df.data %>% group_by(pc, nm) %>% summarize(mR2 = mean(r2, na.rm = T), sdR2 = sd(r2, na.rm = T), mR2.ll = mean(r2.ll, na.rm = T), sdR2.ll = sd(r2.ll, na.rm = T), mR2.cll = mean(r2.ll.1, na.rm = T), sdR2.cll = sd(r2.ll.1, na.rm = T), mAcc = mean(accuracyPercent, na.rm = T), sdAcc = sd(accuracyPercent, na.rm = T), mNS = mean(numberSensitivity, na.rm = T), sdNS = sd(numberSensitivity, na.rm = T), mPI = mean(pIncludeConceptualPoints, na.rm = T), sdPI = sd(pIncludeConceptualPoints, na.rm = T), mFE = mean(firstEstimate, na.rm = T), sdFE = sd(firstEstimate, na.rm = T)))

df.sum.pc <- data.frame(df.data %>% group_by(pc) %>% summarize(mR2 = mean(r2, na.rm = T), sdR2 = sd(r2, na.rm = T), mR2.ll = mean(r2.ll, na.rm = T), sdR2.ll = sd(r2.ll, na.rm = T), mR2.cll = mean(r2.ll.1, na.rm = T), sdR2.cll = sd(r2.ll.1, na.rm = T), mAcc = mean(accuracyPercent, na.rm = T), sdAcc = sd(accuracyPercent, na.rm = T), mNS = mean(numberSensitivity, na.rm = T), sdNS = sd(numberSensitivity, na.rm = T), mPI = mean(pIncludeConceptualPoints, na.rm = T), sdPI = sd(pIncludeConceptualPoints, na.rm = T), mFE = mean(firstEstimate, na.rm = T), sdFE = sd(firstEstimate, na.rm = T)))

df.sum.nm <- data.frame(df.data %>% group_by(nm) %>% summarize(mR2 = mean(r2, na.rm = T), sdR2 = sd(r2, na.rm = T), mR2.ll = mean(r2.ll, na.rm = T), sdR2.ll = sd(r2.ll, na.rm = T), mR2.cll = mean(r2.ll.1, na.rm = T), sdR2.cll = sd(r2.ll.1, na.rm = T), mAcc = mean(accuracyPercent, na.rm = T), sdAcc = sd(accuracyPercent, na.rm = T), mNS = mean(numberSensitivity, na.rm = T), sdNS = sd(numberSensitivity, na.rm = T), mPI = mean(pIncludeConceptualPoints, na.rm = T), sdPI = sd(pIncludeConceptualPoints, na.rm = T), mFE = mean(firstEstimate, na.rm = T), sdFE = sd(firstEstimate, na.rm = T)))


sink("TO3_ParameterAnalysis.txt")
	cat("\nNumber of Unique Subjects:", uniqueSNs, "\n\n")
	cat("\n********** summary stats for nlFit **********\n\n")
	print(df.sum.all)
	print(df.sum.pc)
	print(df.sum.nm)
	print(df.sum.int)

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
	cat("\n\n********** t-test comparing Accuracy Percent for Barth vs Thompson Procedure **********\n\n")
	print(t.test.ap.pc)
	print(cohensD.ap.pc)
	cat("\n\n********** t-test comparing Accuracy Percent for Barth vs Thompson Number Set **********\n\n")
	print(t.test.ap.nm)
	print(cohensD.ap.nm)
	cat("\n\n********** t-test comparing Number Sensitivity for Barth vs Thompson Procedure **********\n\n")
	print(t.test.ns.pc)
	print(cohensD.ns.pc)
	cat("\n\n********** t-test comparing Number Sensitivity for Barth vs Thompson Number Set **********\n\n")
	print(t.test.ns.nm)
	print(cohensD.ns.nm)
sink(NULL)

### examine overall data
pdf("TO3_ParameterFigures.pdf")
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

	with(df.data, plot(accuracyPercent ~ age, frame = F))
	with(df.data, plot(numberSensitivity ~ age, frame = F))
	with(df.data, plot(coef.lam.ll ~ age, frame = F))

	with(df.data, boxplot(accuracyPercent ~ pc, ylim = c(0,1), frame = F))
	with(df.data, boxplot(accuracyPercent ~ nm, ylim = c(0,1), frame = F))
	with(df.data, boxplot(accuracyPercent ~ pc*nm, ylim = c(0,1), frame = F))

	with(df.data, boxplot(numberSensitivity ~ pc, ylim = c(0,1), frame = F))
	with(df.data, boxplot(numberSensitivity ~ nm, ylim = c(0,1), frame = F))
	with(df.data, boxplot(numberSensitivity ~ pc*nm, ylim = c(0,1), frame = F))

dev.off()
setwd(mainDir)

		