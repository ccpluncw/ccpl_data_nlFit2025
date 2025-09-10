library(chutils)
library(dplyr)
library(tidyr)
library(effectsize)

mainDir <- getwd()
ch.newDir (mainDir, "TO_OT_SummaryStats")

datafile <- "./Feb_13_2025_13-53_TO_OT_ParameterData.txt"
df.data <-read.table(datafile, header=T, sep="\t", quote="\"")

df.data <- df.data[df.data$sn != "all", ]
uniqueSNs <- length(unique(df.data$sn))

df.data$r2 <- ifelse(df.data$r2 < 0, 0, df.data$r2)
df.data$r2.ll <- ifelse(df.data$r2.ll < 0, 0, df.data$r2.ll)
df.data$r2.ll.1 <- ifelse(df.data$r2.ll.1 < 0, 0, df.data$r2.ll.1)

#### Examine pre and post influence
df.data.wide <- df.data  %>% pivot_wider (names_from = c(ev),
    values_from = c(firstEstimate, accuracyPercent, numberSensitivity, pIncludeConceptualPoints, r2, BIC, coef.int.ll, coef.lam.ll, r2.ll, BIC.ll, coef.lam.ll.1, r2.ll.1, BIC.ll.1))
		
df.data.wide$accPer_preVsPost <- df.data.wide$accuracyPercent_post - df.data.wide$accuracyPercent_pre
df.data.wide$numSen_preVsPost <- df.data.wide$numberSensitivity_post - df.data.wide$numberSensitivity_pre
df.data.wide$firstEst_preVsPost <- df.data.wide$firstEstimate_post - df.data.wide$firstEstimate_pre

### examine overall data
t.test.r2.nlVll <- with(df.data, t.test(r2, r2.ll, paired = T))
t.test.r2.nlVll.1 <- with(df.data, t.test(r2 , r2.ll.1, paired = T))

t.test.BIC.nlVll <- with(df.data, t.test(BIC , BIC.ll, paired = T))
t.test.BIC.nlVll.1 <-with(df.data, t.test(BIC , BIC.ll.1, paired = T))

cohensD.r2.nlVll <- with(df.data, cohens_d(r2, r2.ll, paired = T))
cohensD.r2.nlVll.1 <- with(df.data, cohens_d(r2 , r2.ll.1, paired = T))

cohensD.BIC.nlVll <- with(df.data, cohens_d(BIC , BIC.ll, paired = T))
cohensD.BIC.nlVll.1 <-with(df.data, cohens_d(BIC , BIC.ll.1, paired = T))

t.test.numSens.preVpost <- with(df.data.wide, t.test(numSen_preVsPost ~ tr, paired = F))
t.test.accPerc.preVpost <-with(df.data.wide, t.test(accPer_preVsPost ~ tr, paired = F))
t.test.firstEst.preVpost <-with(df.data.wide, t.test(firstEst_preVsPost ~ tr, paired = F))
t.test.pi <- with(df.data[df.data$tr == 1 & df.data$ev == "post", ], t.test(pIncludeConceptualPoints, mu = 0, paired = F))

cohensD.numSens.preVpost <- with(df.data.wide, cohens_d(numSen_preVsPost ~ tr, paired = F))
cohensD.accPerc.preVpost <-with(df.data.wide, cohens_d(accPer_preVsPost ~ tr, paired = F))
cohensD.firstEst.preVpost <-with(df.data.wide, cohens_d(firstEst_preVsPost ~ tr, paired = F))
cohensD.pi <- with(df.data[df.data$tr == 1 & df.data$ev == "post", ], cohens_d(pIncludeConceptualPoints, mu = 0, paired = F))


df.sum.all <- data.frame(df.data %>% summarize(mR2 = mean(r2, na.rm = T), sdR2 = sd(r2, na.rm = T), mR2.ll = mean(r2.ll, na.rm = T), sdR2.ll = sd(r2.ll, na.rm = T), mR2.cll = mean(r2.ll.1, na.rm = T), sdR2.cll = sd(r2.ll.1, na.rm = T), mAcc = mean(accuracyPercent, na.rm = T), sdAcc = sd(accuracyPercent, na.rm = T), mNS = mean(numberSensitivity, na.rm = T), sdNS = sd(numberSensitivity, na.rm = T), mPI = mean(pIncludeConceptualPoints, na.rm = T), sdPI = sd(pIncludeConceptualPoints, na.rm = T), mFE = mean(firstEstimate, na.rm = T), sdFE = sd(firstEstimate, na.rm = T)))

df.sum.treat <- data.frame(df.data %>% group_by(tr, ev) %>% summarize(mR2 = mean(r2, na.rm = T), sdR2 = sd(r2, na.rm = T), mR2.ll = mean(r2.ll, na.rm = T), sdR2.ll = sd(r2.ll, na.rm = T), mR2.cll = mean(r2.ll.1, na.rm = T), sdR2.cll = sd(r2.ll.1, na.rm = T), mAcc = mean(accuracyPercent, na.rm = T), sdAcc = sd(accuracyPercent, na.rm = T), mNS = mean(numberSensitivity, na.rm = T), sdNS = sd(numberSensitivity, na.rm = T), mPI = mean(pIncludeConceptualPoints, na.rm = T), sdPI = sd(pIncludeConceptualPoints, na.rm = T), mFE = mean(firstEstimate, na.rm = T), sdFE = sd(firstEstimate, na.rm = T)))

df.sum.p1 <- data.frame(df.data[df.data$tr == 1 & df.data$ev == "post", ] %>% summarize(mPI = mean(pIncludeConceptualPoints, na.rm = T), sdPI = sd(pIncludeConceptualPoints, na.rm = T)))

df.sum.ev <- data.frame(df.data %>% group_by(ev) %>% summarize(mR2 = mean(r2, na.rm = T), sdR2 = sd(r2, na.rm = T), mR2.ll = mean(r2.ll, na.rm = T), sdR2.ll = sd(r2.ll, na.rm = T), mR2.cll = mean(r2.ll.1, na.rm = T), sdR2.cll = sd(r2.ll.1, na.rm = T), mAcc = mean(accuracyPercent, na.rm = T), sdAcc = sd(accuracyPercent, na.rm = T), mNS = mean(numberSensitivity, na.rm = T), sdNS = sd(numberSensitivity, na.rm = T), mPI = mean(pIncludeConceptualPoints, na.rm = T), sdPI = sd(pIncludeConceptualPoints, na.rm = T), mFE = mean(firstEstimate, na.rm = T), sdFE = sd(firstEstimate, na.rm = T)))

sink("TO1_ParameterAnalysis.txt")
	cat("\nNumber of Unique Subjects:", uniqueSNs, "\n\n")

	cat("\n********** summary stats for nlFit **********\n\n")
	print(df.sum.all)
	print(df.sum.p1)
	print(df.sum.treat)
	print(df.sum.ev)
	cat("\n********** t-test comparing pIncludeConceptualPoints to 0 **********\n\n")
	print(t.test.pi)
	print(cohensD.pi)
	cat("\n\n********** t-test comparing Number Sensitivity pre vs post treatement **********\n\n")
	print(t.test.numSens.preVpost)
	print(cohensD.numSens.preVpost)
	cat("\n\n********** t-test comparing Accuracy Percent pre vs post treatement **********\n\n")
	print(t.test.accPerc.preVpost)
	print(cohensD.accPerc.preVpost)
	cat("\n\n********** t-test comparing First Estimate pre vs post treatement **********\n\n")
	print(t.test.firstEst.preVpost)
	print(cohensD.firstEst.preVpost)

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

### examine overall data
pdf("TO1_ParameterFigures.pdf")
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


	### within subject analysis
	with(df.data.wide[df.data.wide$tr == 1, ], plot(accuracyPercent_post ~ accuracyPercent_pre, xlim = c(0, 1), ylim = c(0,1), frame = F, main = "treatment 1"))
	abline(0,1)
	with(df.data.wide[df.data.wide$tr == 0, ], plot(accuracyPercent_post ~ accuracyPercent_pre, xlim = c(0, 1), ylim = c(0,1), frame = F, main = "treatment 0"))
	abline(0,1)
		
	with(df.data.wide, boxplot(accuracyPercent_post ~ tr, ylim = c(0,1), frame = F))
	with(df.data.wide, boxplot(accPer_preVsPost ~ tr, ylim = c(-1,1), frame = F))
	abline(h = 0, col = "grey", lwd = 2)
		

	with(df.data.wide[df.data.wide$tr == 1, ], plot(numberSensitivity_post ~ numberSensitivity_pre, xlim = c(0, 1), ylim = c(0,1), frame = F, main = "treatment 1"))
	abline(0,1)
	with(df.data.wide[df.data.wide$tr == 0, ], plot(numberSensitivity_post ~ numberSensitivity_pre, xlim = c(0, 1), ylim = c(0,1), frame = F, main = "treatment 0"))
	abline(0,1)
		
	with(df.data.wide, boxplot(numberSensitivity_post ~ tr, ylim = c(0,1), frame = F))
		
	with(df.data.wide, boxplot(numSen_preVsPost ~ tr, ylim = c(-1,1), frame = F))
	abline(h = 0, col = "grey", lwd = 2)


		
	with(df.data.wide, boxplot(firstEstimate_post ~ tr, frame = F))
		
	with(df.data.wide, boxplot(firstEst_preVsPost ~ tr, frame = F))
	abline(h = 0, col = "grey", lwd = 2)
	
dev.off()
setwd(mainDir)
		