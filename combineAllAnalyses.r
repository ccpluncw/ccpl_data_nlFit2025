library(chutils)
library(dplyr)
library(tidyr)

mainDir <- getwd()

statFileNames <- c("./Thompson2/TO3_SummaryStats/Feb_13_2025_13-56_TO_OT_ParameterData.txt", "./Thompson2/TO4_SummaryStats/Feb_13_2025_13-57_TO4_ParameterData.txt", "./Cohen/CO_SummaryStats/Feb_13_2025_13-54_CO_ParameterData.txt", "./T&O data/TO_OT_SummaryStats/Feb_13_2025_13-53_TO_OT_ParameterData.txt")

keepColumns <- c("sn", "firstEstimate", "accuracyPercent", "numberSensitivity", "pIncludeConceptualPoints", "rangeLength", "memoryLength", "r2", "BIC", "coef.int.ll", "coef.lam.ll", "r2.ll", "BIC.ll", "coef.lam.ll.1", "r2.ll.1", "BIC.ll.1")

df.data <- NULL
for(fn in statFileNames) {
	df.data.tmp <-read.table(fn, header=T, sep="\t", quote="\"")
	if(basename(dirname(fn)) == "TO4_SummaryStats") {
		df.data.tmp$sn <- ifelse(df.data.tmp$sn == "all", paste(df.data.tmp$sn, df.data.tmp$gd, sep="_"), df.data.tmp$sn)		
#		df.data.tmp$sn <- paste(df.data.tmp$sn, df.data.tmp$gd, sep="_")	
	}
	if(basename(dirname(fn)) == "TO_OT_SummaryStats") {
		df.data.tmp$sn <- ifelse(df.data.tmp$sn == "all", paste(df.data.tmp$sn, df.data.tmp$ex, df.data.tmp$ev, df.data.tmp$tr, sep="_"), df.data.tmp$sn)		
#		df.data.tmp$sn <- paste(df.data.tmp$sn, df.data.tmp$ex, df.data.tmp$ev, df.data.tmp$tr, sep="_")
	}
	if(basename(dirname(fn)) == "CO_SummaryStats") {
		df.data.tmp$sn <- ifelse(df.data.tmp$sn == "all", paste(df.data.tmp$sn, df.data.tmp$cnd, sep="_"), df.data.tmp$sn)		
#		df.data.tmp$sn <- paste(df.data.tmp$sn, df.data.tmp$cnd, sep="_")
	}
	if(basename(dirname(fn)) == "TO3_SummaryStats") {
		df.data.tmp$sn <- ifelse(df.data.tmp$sn == "all", paste(df.data.tmp$sn, df.data.tmp$pc,df.data.tmp$nm, sep="_"), df.data.tmp$sn)		
#		df.data.tmp$sn <- paste(df.data.tmp$sn, df.data.tmp$pc,df.data.tmp$nm, sep="_")
	}
	
	df.data.tmp <- df.data.tmp[,keepColumns]
	df.data.tmp$exp <- basename(dirname(fn))
	df.data <- ch.rbind(df.data, df.data.tmp)
}

ch.newDir (mainDir, "all_SummaryStats")

#fileTag <- paste(format(Sys.time(), "%b_%d_%Y_%H-%M"), sep="_")
fileName <- "All_ParameterData.txt"
write.table(df.data, fileName, append = FALSE, col.names=T, row.names=F, quote=F, sep="\t")

setwd(mainDir)
