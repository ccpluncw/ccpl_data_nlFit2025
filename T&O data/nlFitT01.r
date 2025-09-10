library(chutils)
library(dplyr)
library(smartGridSearch)
library(nlFit)
library(dplyr)
library(tidyr)

runGroupData <- FALSE
upperBound <- 1000
lowerBound <- 0
visibleReferencePoints <- c(lowerBound, upperBound)
conceptualReferencePoints <- NULL

loops <- 500
finalLoops <- 500

testCode <- TRUE

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

start <- Sys.time()

fileTag <- paste(format(Sys.time(), "%b_%d_%Y_%H-%M"), "NLFit", sep="_")


datafile <- "TOdataAll.txt"
data.in <-read.table(datafile, header=T, sep="\t", quote="\"", check.names = FALSE)

exps <- unique(data.in$exp)
events <- unique(data.in$event)
treatments <- unique(data.in$Treatment)

dataTargetCol <- "target"
dataEstimateCol <- "estimate"
freeMemoryLength <- FALSE
freeRangeLength <- FALSE
freeNumberSensitivity <- TRUE
freepIncludeConceptualPoints <- FALSE
directoryTag <- "final.v12"

lowerBound <- 0
upperBound <- 1000

defaultRangeLength <- 1
defaultNumberSensitivity <- 1

data.long <- data.in %>% pivot_longer (cols = c(2:23), names_to = "target", values_to = "estimate")
data.long[[dataTargetCol]] <- as.numeric(as.character(data.long[[dataTargetCol]]))


df.out <- NULL
df.dat <- NULL

mainDir <- getwd()


if(testCode) {
	ch.newDir (mainDir, paste( "nlFit", directoryTag, sep="_"))
} else {
	ch.newDir (mainDir, "nlFit")
}
	sinkFileName <- paste(fileTag, "TO1_nlFitout.txt", sep="_")
	
	sink(sinkFileName, append = F)
		cat("\n\n\n************* NEW RUN **************\n\n")
	sink(NULL)

nlFitDir <- getwd()
ch.newDir (nlFitDir, "plot")
plotDir <- getwd()
setwd(nlFitDir)
ch.newDir (nlFitDir, "smartGridOutput")
sgOutDir <- getwd()


for(ex in exps) {
	df.tmp.1 <- data.long[data.long$exp == ex, ]
	for(ev in events) {
		df.tmp.2 <- df.tmp.1[df.tmp.1$event == ev, ]
		for (tr in treatments) {
			df.tmp.3 <- df.tmp.2[df.tmp.2$Treatment == tr, ]
			
			if(runGroupData) {
				subs <- "all"
			} else {
				subs <- c(unique(df.tmp.3$sn), "all")
			}
			
			numTargets <- length(unique(df.tmp.3[[dataTargetCol]]))
			
			if(ev == "post" & tr == 1) {
				# add midpoint
				conceptualReferencePoints <- c(150, 147, 179, 187, 156,	163, 172)
				freepIncludeConceptualPoints <- TRUE
			} else {
				conceptualReferencePoints <- NULL
				freepIncludeConceptualPoints <- FALSE
			}
			
			for(sn in subs) {
				if(sn == "all") {
					df.tmp <- df.tmp.3 %>% group_by(target) %>% summarize(estimate = mean(estimate, na.rm = T), Age = mean(Age, na.rm = t))
				} else {
					df.tmp <- df.tmp.3[df.tmp.3$sn == sn, ]
				}
				
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
			    cat("\n\n ***********************", sn, ex, ev, tr, age, "*********************** \n\n")
				sink(NULL)
				outList <- outputNLfitStats(df.tmp, statList, dataTargetCol = "target", dataEstimateCol = "estimate", pars.n = pars.n, loops = 1000, sinkFilename = paste("..", sinkFileName, sep = "/"), appendSinkFile = TRUE, plotFileName = paste(sn, ex, ev, tr,"Plot.pdf", sep="_"))				
				##########
				
				####### compare to LL #######
				df.tmp$logTarget <- (visibleReferencePoints[2]/log(visibleReferencePoints[2])) * log(df.tmp$target) 
				fit.ll <- nls(estimate ~ a*( ((1-lambda)*target) + (lambda*logTarget)), data=df.tmp, start=list(a=1, lambda = 1), 					algorithm="port", lower=c(-Inf,0), upper=c(Inf,1), control = nls.control(minFactor=1/10000000, maxiter=10000, warnOnly = FALSE))
				df.tmp$fit.ll <- predict(fit.ll, newdata = df.tmp)
				r2.ll <- round(ch.R2(df.tmp$estimate, df.tmp$fit.ll, standardize = FALSE),2)
				BIC.ll <- round(BIC(fit.ll),1)
				
				############################
							
				df.out.tmp <- data.frame(sn, ex, ev, tr, age, firstEstimate = outList$runStats$parameters$firstEstimate, accuracyPercent = outList$runStats$parameters$accuracyPercent, numberSensitivity = outList$runStats$parameters$numberSensitivity, pIncludeConceptualPoints = outList$runStats$parameters$pIncludeConceptualPoints, rangeLength = outList$runStats$parameters$rangeLength, memoryLength = outList$runStats$parameters$memoryLength, r2 = outList$runStats$fitStats$r2, BIC = outList$runStats$fitStats$BIC, r2.ll, BIC.ll)
				
				df.out <- ch.rbind(df.out, df.out.tmp)				
			}
		}
	}
}

setwd(nlFitDir)
fileName <- paste(fileTag, "snFits.txt", sep="_")
write.table(df.out, fileName, append = FALSE, col.names=T, row.names=F, quote=F, sep="\t")
setwd(mainDir)