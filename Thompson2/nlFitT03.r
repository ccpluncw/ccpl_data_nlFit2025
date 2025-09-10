library(chutils)
library(dplyr)
library(smartGridSearch)
library(nlFit)
library(dplyr)
library(tidyr)
library(readxl)

runGroupData <- FALSE
loops <- 500
finalLoops <- 500

test <- TRUE
useMultiCore <- TRUE
multicorePackages <- "nlFit"


#number of intervals to split each bounded parameter into in order to get the random values between the hi and low values
numGridIntervals <- 200
#the number of times to loop with grid search before looking for the top 10 best fits.
numGridLoops <- ifelse (test, 250, 2500)
#paramtable confirmation loops
optBoundLoops <- 10
#The number of best fit runs from which to resize the parameter bounds.
#Here, I set it to the sqrt(total number of runs per grid loop) or 10 at minimum.
  #optParamListN <- ifelse ( trunc(sqrt(numGridLoops)) < 10, 10, trunc(sqrt(numGridLoops)) )
optParamListN <- 10
#the number of simulations to average using the final parameters.
numSimsToAverage <- ifelse (test, 20, 40)

numLoops <- numGridLoops
numIntervals <- numGridIntervals
optParamListN <- optParamListN
optBoundLoops <- optBoundLoops

dataTargetCol <- "target"
dataEstimateCol <- "estimate"
freeMemoryLength <- FALSE
freeRangeLength <- FALSE
freeNumberSensitivity <- TRUE
directoryTag <- "final.v12"

defaultRangeLength <- 1
defaultNumberSensitivity <- 1

start <- Sys.time()

fileTag <- paste(format(Sys.time(), "%b_%d_%Y_%H-%M"), "NLfit", sep="_")

datafile <- "TO3data.xlsx"
#read in the sheets of the clinician production excel file
sheets <- excel_sheets(datafile)

#combine all sheets into one large dataset
df.data.all <- list()
for(i in 1:length(sheets)) {
  #read in the clinicians billed units
  df.data.all[[sheets[i]]] <- data.frame(read_excel(datafile, sheet = sheets[i], .name_repair = "minimal"), check.names = FALSE)
}

mainDir <- getwd()

for(st in sheets) {
	data.in <- df.data.all[[st]]
	columnNames <- colnames(data.in)
	
	targetColumns <- columnNames[grepl("^[[:digit:]]+", columnNames)]
	
	procs <- unique(data.in$proc)
	nums <- unique(data.in$num)

	data.long <- data.in %>% pivot_longer (cols = all_of(targetColumns), names_to = dataTargetCol, values_to = dataEstimateCol)
	data.long$target <- as.numeric(as.character(data.long$target))

	df.out <- NULL
	df.dat <- NULL
	
	lowerBound <- 0
	upperBound <- 1000
	visibleReferencePoints <- c(lowerBound, upperBound)
	
	if(test) {
		ch.newDir (mainDir, paste(st, directoryTag, sep="_"))
	} else {
		ch.newDir (mainDir, st)
	}
	expDir <- getwd()
	
	ch.newDir (expDir, "nlFit")
		sinkFileName <- paste(fileTag, st, "nlFitout.txt", sep="_")
		sink(sinkFileName, append = F)
			cat("\n\n\n************* NEW RUN **************\n\n")
		sink(NULL)

	nlFitDir <- getwd()
	ch.newDir (nlFitDir, "plot")
	plotDir <- getwd()
	setwd(nlFitDir)
	ch.newDir (nlFitDir, "smartGridOutput")
	sgOutDir <- getwd()

	for(pc in procs) {
		df.tmp.1 <- data.long[data.long$proc == pc, ]
		if(pc == "S") {
			# add midpoint
			conceptualReferencePoints <- c(0.5*upperBound)
			freepIncludeConceptualPoints <- TRUE
		} else {
			conceptualReferencePoints <- NULL
			freepIncludeConceptualPoints <- FALSE
		}

		for(nm in nums) {
			df.tmp.2 <- df.tmp.1[df.tmp.1$num == nm, ]
	
			if(runGroupData) {
				subs <- "all"
			} else {
				subs <- c(unique(df.tmp.2$sn), "all")
			}
	
			for(sn in subs) {
				if(sn == "all") {
					df.tmp <- df.tmp.2 %>% group_by(target) %>% summarize(estimate = mean(estimate, na.rm = T), Age = mean(Age, na.rm = t))
				} else {
					df.tmp <- df.tmp.2[df.tmp.2$sn == sn, ]
				}
				numTargets <- length(unique(df.tmp$target))
				
				setwd(sgOutDir)
				
				age = round(mean(df.tmp$Age, na.rm = T), 2)
				
				### nlFit
				pUpper <- list(firstEstimate = upperBound, accuracyPercent = 1)
				pLower <- list(firstEstimate = lowerBound, accuracyPercent = 0)
				pInt <- list(firstEstimate = 0.1, accuracyPercent = 0.001)
				
				otherParameterList <- list(data = df.tmp, upperBound = upperBound, lowerBound = lowerBound,visibleReferencePoints = visibleReferencePoints, conceptualReferencePoints = conceptualReferencePoints, loops = loops, targetOrder = "random", dataTargetCol = "target", dataEstimateCol = "estimate", minimizeStat = 'BIC')
				
				if(freeNumberSensitivity) {
					pUpper[["numberSensitivity"]] <- 1
					pLower[["numberSensitivity"]] <- 0
					pInt[["numberSensitivity"]] <- 0.001
				} else {
					otherParameterList[["numberSensitivity"]] <- defaultNumberSensitivity
				}
				if(freeRangeLength) {
					pUpper[["rangeLength"]] <- numTargets
					pLower[["rangeLength"]] <- 1
					pInt[["rangeLength"]] <- 1
				} else {
					otherParameterList[["rangeLength"]] <- defaultRangeLength
				}
				if(freeMemoryLength) {
					pUpper[["memoryLength"]] <- numTargets
					pLower[["memoryLength"]] <- 1
					pInt[["memoryLength"]] <- 1
				} else {
					otherParameterList[["memoryLength"]] <- numTargets
				}
				if(freepIncludeConceptualPoints) {
					pUpper[["pIncludeConceptualPoints"]] <- 1
					pLower[["pIncludeConceptualPoints"]] <- 0
					pInt[["pIncludeConceptualPoints"]] <- 0.001
				} else {
					otherParameterList[["pIncludeConceptualPoints"]] <- 0
				}
				pars.n <- length(pUpper)
				otherParameterList[["pars.n"]] <- pars.n
				
				x.grid <- smartGridSearch(getOrdinalNumberlineFit, pUpper, pLower, pInt, otherParameterList, numLoops = numLoops, numIntervals = numIntervals, optParamListN = optParamListN, optBoundLoops = optBoundLoops, multicore = useMultiCore, multicorePackages = multicorePackages)
				
				statList <- x.grid$final
				statList$upperBound <- upperBound
				statList$lowerBound <- lowerBound
				if(!freeMemoryLength) statList$memoryLength <- numTargets
				if(!freeRangeLength) statList$rangeLength <- defaultRangeLength
				if(!freepIncludeConceptualPoints) statList$pIncludeConceptualPoints <- 0
				if(!freeNumberSensitivity) statList$numberSensitivity <- defaultNumberSensitivity
				statList$targetOrder <- "random"
				statList$visibleReferencePoints <- visibleReferencePoints
				statList$conceptualReferencePoints <- conceptualReferencePoints
				statList$loops <- finalLoops

				setwd(plotDir)
			  sink(paste("..", sinkFileName, sep = "/"), append = TRUE)
			    cat("\n\n ***********************", st, sn, pc, nm, age, "*********************** \n\n")
				sink(NULL)
				outList <- outputNLfitStats(df.tmp, statList, dataTargetCol = "target", dataEstimateCol = "estimate", pars.n = pars.n, loops = 1000, sinkFilename = paste("..", sinkFileName, sep = "/"), appendSinkFile = TRUE, plotFileName = paste(st, sn, pc, nm,"Plot.pdf", sep="_"))
								
				##########

				####### compare to LL #######
				df.tmp$logTarget <- (visibleReferencePoints[2]/log(visibleReferencePoints[2])) * log(df.tmp$target) 
				fit.ll <- nls(estimate ~ a*( ((1-lambda)*target) + (lambda*logTarget)), data=df.tmp, start=list(a=1, lambda = 1), 					algorithm="port", lower=c(-Inf,0), upper=c(Inf,1), control = nls.control(minFactor=1/10000000, maxiter=10000, warnOnly = FALSE))
				df.tmp$fit.ll <- predict(fit.ll, newdata = df.tmp)
				r2.ll <- round(ch.R2(df.tmp$estimate, df.tmp$fit.ll, standardize = FALSE),2)
				BIC.ll <- round(BIC(fit.ll),1)
				
				############################

				df.out.tmp <- data.frame(sn, st, pc, nm, age, firstEstimate = outList$runStats$parameters$firstEstimate, accuracyPercent = outList$runStats$parameters$accuracyPercent, numberSensitivity = outList$runStats$parameters$numberSensitivity, pIncludeConceptualPoints = outList$runStats$parameters$pIncludeConceptualPoints, rangeLength = outList$runStats$parameters$rangeLength, memoryLength = outList$runStats$parameters$memoryLength, r2 = outList$runStats$fitStats$r2, BIC = outList$runStats$fitStats$BIC, r2.ll, BIC.ll)
				
				df.out <- ch.rbind(df.out, df.out.tmp)				
			}
		}
	}
	setwd(expDir)
	fileName <- paste(fileTag, "snFits.txt", sep="_")
	write.table(df.out, fileName, append = FALSE, col.names=T, row.names=F, quote=F, sep="\t")
	setwd(mainDir)
}

