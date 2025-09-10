library(chutils)
library(dplyr)
library(tidyr)

library(emmeans)

mainDir <- getwd()
ch.newDir (mainDir, "all_SummaryStats")

datafile <- "All_ParameterData.txt"
df.data <-read.table(datafile, header=T, sep="\t", quote="\"")

#remove collapsed data
df.SNs <- df.data[, c("sn", "exp")]
groups <-sum(startsWith(df.SNs$sn, "all") )
df.SNs <- df.SNs[!startsWith(df.SNs$sn, "all"), ]
uniqueSN <- nrow(unique(df.SNs))
duplicateSN <- nrow(df.SNs[duplicated(df.SNs),])
totalObs <- nrow(df.data)

#remove group data before analysis
df.data <- df.data[!startsWith(df.data$sn, "all"), ]

### examine overall data
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

lm.ap_ns <- with(df.data, lm(accuracyPercent ~ numberSensitivity))
lm.ap_pi <- with(df.data[df.data$pIncludeConceptualPoints > 0, ], lm(accuracyPercent ~ pIncludeConceptualPoints))

lm.ap_lam.ll <- with(df.data, lm(accuracyPercent ~ coef.lam.ll))
lm.ap_int.ll <- with(df.data, lm(accuracyPercent ~ coef.int.ll))

lm.ns_lam.ll <- with(df.data, lm(numberSensitivity ~ coef.lam.ll))
lm.ns_int.ll <- with(df.data, lm(numberSensitivity ~ coef.int.ll))


sink("All_ParameterAnalysis.txt")
	cat("\n Total Observation:", totalObs, "\n")
	cat("\n Unique Participants:", uniqueSN, "\n")
	cat("\n Participants in 2 conditions:", duplicateSN, "\n")
	cat("\n grouped condition:", groups, "\n")
	
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
	cat("\n\n********** regression relating accuracy percent and number sensitivity **********\n\n")
	print(summary(lm.ap_ns))
	cat("\n\n********** regression relating accuracy percent and pIncludeConceptualPoints **********\n\n")
	print(summary(lm.ap_pi))
	cat("\n\n********** regression relating accuracy percent and lambda from LL **********\n\n")
	print(summary(lm.ap_lam.ll))
	cat("\n\n********** regression relating accuracy percent and constant from LL **********\n\n")
	print(summary(lm.ap_int.ll))
	cat("\n\n********** regression relating number sensitivity and lambda from LL **********\n\n")
	print(summary(lm.ns_lam.ll))
	cat("\n\n********** regression relating number sensitivity and constant from LL **********\n\n")
	print(summary(lm.ns_int.ll))
sink(NULL)

### examine overall data
pdf("All_ParameterFigures.pdf")
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
	with(df.data, boxplot(r2 , r2.ll, r2.ll.1, ylim = c(0,1), names = c("nlFit", "Less", "More"), ylab = "r2", frame = F))
	with(df.data, boxplot(BIC , BIC.ll, names = c("nlFit", "LL"), ylab = "BIC", frame = F))
	with(df.data, boxplot(BIC , BIC.ll.1, names = c("nlFit", "Constrained LL"), ylab = "BIC", frame = F))
	with(df.data, boxplot(BIC , BIC.ll, BIC.ll.1, names = c("nlFit", "Less", "More"), ylab = "BIC", frame = F))

	with(df.data, plot(accuracyPercent ~ numberSensitivity, xlim = c(0,1), ylim = c(0,1), frame = F))
	abline(lm.ap_ns, col = "grey")
	
	with(df.data[df.data$pIncludeConceptualPoints > 0, ], plot(accuracyPercent ~ pIncludeConceptualPoints, xlim = c(0,1), ylim = c(0,1), frame = F))
	abline(lm.ap_pi, col = "grey")

	with(df.data, plot(accuracyPercent ~ coef.lam.ll, frame = F))
	abline(lm.ap_lam.ll, col = "grey")
	with(df.data, plot(accuracyPercent ~ coef.int.ll, frame = F))
	abline(lm.ap_int.ll, col = "grey")

	with(df.data, plot(numberSensitivity ~ coef.lam.ll, frame = F))
	abline(lm.ns_lam.ll, col = "grey")
	with(df.data, plot(numberSensitivity ~ coef.int.ll, frame = F))
	abline(lm.ns_int.ll, col = "grey")
	
	plot(NULL, frame = F, xlim = c(0,100), ylim = c(-1000, 1000))
	xVal <- seq (0, 100, 1)
	const <- seq(-8, 8, 2)
	lambdas <- c(0, 0.5, 1)
	for(ll in lambdas) {
		for (a in const) {
			yVal <- a * ( ((1 - ll) * xVal) + (ll*(100/log(100)) * log(xVal)))
			lines(xVal, yVal, col = "black")
		}
	}

	plot(NULL, frame = F, xlim = c(0,100), ylim = c(0, 100))
	xVal <- seq (0, 100, 1)
	const <- seq(1)
	lambdas <- seq(0,1, 0.1)
	for(ll in lambdas) {
		for (a in const) {
			yVal <- a * ( ((1 - ll) * xVal) + (ll*(100/log(100)) * log(xVal)))
			lines(xVal, yVal, col = "black")
		}
	}
	
	
dev.off()
setwd(mainDir)


