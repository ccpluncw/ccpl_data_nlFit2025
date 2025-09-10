library(chutils)
library(dplyr)
library(tidyr)

mainDir <- getwd()
ch.newDir (mainDir, "CO_SummaryStats")

datafile <- "Feb_13_2025_13-54_CO_ParameterData.txt"
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

t.test.acc <- with(df.data, t.test(accuracyPercent ~ cnd, paired = F))
t.test.ns <- with(df.data, t.test(numberSensitivity ~ cnd, paired = F))
t.test.pi <- with(df.data[df.data$cnd == "bounded_0", ], t.test(pIncludeConceptualPoints, mu = 0, paired = F))

cohensD.acc <- with(df.data, cohens_d(accuracyPercent ~ cnd, paired = F))
cohensD.ns <- with(df.data, cohens_d(numberSensitivity ~ cnd, paired = F))
cohensD.pi <- with(df.data[df.data$cnd == "bounded_0", ], cohens_d(pIncludeConceptualPoints, mu = 0, paired = F))

df.sum <- data.frame(df.data %>% group_by(cnd) %>% summarize(mR2 = mean(r2, na.rm = T), sdR2 = sd(r2, na.rm = T), mR2.ll = mean(r2.ll, na.rm = T), sdR2.ll = sd(r2.ll, na.rm = T), mR2.cll = mean(r2.ll.1, na.rm = T), sdR2.cll = sd(r2.ll.1, na.rm = T), mAcc = mean(accuracyPercent, na.rm = T), sdAcc = sd(accuracyPercent, na.rm = T), mNS = mean(numberSensitivity, na.rm = T), sdNS = sd(numberSensitivity, na.rm = T), mPI = mean(pIncludeConceptualPoints, na.rm = T), sdPI = sd(pIncludeConceptualPoints, na.rm = T)))


sink("CO_ParameterAnalysis.txt")
	cat("\nNumber of Unique Subjects:", uniqueSNs, "\n\n")
	cat("\n********** summary stats for nlFit **********\n\n")
	print(df.sum)
	cat("\n********** t-test comparing Accuracy Percent by task **********\n\n")
	print(t.test.acc)
	print(cohensD.acc)
	cat("\n********** t-test comparing Number Sensitivity by task **********\n\n")
	print(t.test.ns)
	print(cohensD.ns)
	cat("\n********** t-test comparing pIncludeConceptualPoints to 0 **********\n\n")
	print(t.test.pi)
	print(cohensD.pi)
	
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
pdf("CO_ParameterFigures.pdf")
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

	with(df.data, boxplot(accuracyPercent ~ cnd, ylim = c(0,1), ylab = "accuracyPercent", frame = F))
	with(df.data, boxplot(numberSensitivity ~ cnd, ylim = c(0,1),ylab = "numberSensitivity", frame = F))
dev.off()

setwd(mainDir)

		
		