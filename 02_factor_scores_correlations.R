
# --------------------------------------------------------------------- #
## factor scores correlations
## simulated data
## script: Daniel Bustamante, PhD
## Harvard University & Broad Institute of MIT and Harvard
## Feb 2025
# --------------------------------------------------------------------- #

rm(list=ls())

# load packages
library(OpenMx)
library(data.table)
library(psych)


mxOption(key='Number of Threads', value=parallel::detectCores()) #now
mxOption(key="Default optimizer", value="SLSQP") ## NPSOL, CSOLNP, SLSQP


# read in the data (simulated data are already residualized, their outliers removed, and standardized)
theDataKe <- as.data.frame(fread("sim_LEC_PCL_FSs_Kenya.csv"))
theDataUg <- as.data.frame(fread("sim_LEC_PCL_FSs_Uganda.csv"))


# ----------------------------------------------------------------------
# correlations via SEM-ML
## change optimizer if needed below
# mxOption(key="Default optimizer", value="CSOLNP") ## change optimizer to NPSOL, CSOLNP, SLSQP


# run SEM
# define confidence interval
ciInterval <- 0.95

# set Starting Values for models 
svMe       <- 1                 # start value for means
svVar      <- 1                 # start value for variance
svCov      <- .1                # start value for covariance


# define model objects
subSample <- "Kenya" # Kenya or Uganda
modelData <- theDataKe # ---------------------------------- change dataframe for corresponding sample 
colnames(modelData)

outDataPath <- "./outSim/"
filename <- "sim_results_"
phenoName <- "LEC_PTEs_PCL_PTSDsxs"
covAdj <- "sexAge"

# create directory for output
dir.create(outDataPath)

# declare variables to run in the loop
vars2run4loop <- grep("rFSs", colnames(modelData), value = TRUE)
vars2run4loop


# start of loop (runs one variable against n other variables)
for (i in 1:length(vars2run4loop)) {  
  x1 <- vars2run4loop[i] 
  
  
  for (k in 1:length(vars2run4loop)) {
    x2 <- vars2run4loop[k]
    
    tryCatch({
      
      # set the OpenMx model  
      manifests <- c(x1, x2)
      manifests
      latents <- c()
      
      model <- mxModel("fs2fs_Model", 
                       type="RAM",
                       manifestVars = manifests,
                       latentVars = latents,
                       mxPath(from="one",to=c(x2, x1), free=c(TRUE,TRUE), value=c(svMe, svMe) , arrows=1, label=c(paste0("const__", x2), paste0("const__", x1)) ),
                       
                       mxPath(from=x1,to=c(x1,x2), free=c(TRUE,TRUE), value=c(svVar, svCov) , arrows=2, label=c(paste0("VAR_", x1), paste0("COV_", x1, "_", x2)) ),
                       mxPath(from=x2,to=c(x2), free=c(TRUE), value=c(svVar) , arrows=2, label=c(paste0("VAR_", x2)) ),
                       
                       CIs <- mxCI( c(paste0("VAR_", x1), paste0("COV_", x1, "_", x2), paste0("VAR_", x2), paste0("const__", x1), paste0("const__", x2) ), interval = ciInterval),
                       mxData(modelData, type = "raw"),
                       mxFitFunctionML()
      );
      
      # run it or TryHard it
      result <- mxTryHard(model, intervals = TRUE, extraTries = 20)
      summary(result)
      
      
      # set objects for the table
      # CIs
      SE <- as.vector(result$output$standardErrors[2,])
      SE
      corrVars <- paste0(x1,"_",x2)
      corrVars
      EstCI <- result$output$confidenceIntervals[2,]
      estCImx <- matrix(EstCI, nrow = 1, ncol = 3)
      estCImx
      estCIdf <- data.frame(estCImx)
      estCIdf
      theLine <- cbind(corrVars, estCIdf, SE)
      theLine
      colnames(theLine) <- c( "vars", paste0("lboundCI0p", substr(round(ciInterval,3),3,5)), "corr", paste0("uboundCI0p", substr(round(ciInterval,3),3,5)), "se")
      theLine
      
      statsLine <- theLine
      statsLine
      
      theColNames <- t(colnames(statsLine))
      theColNames
      
      # save the table, append each iteration per row in the .csv with no colnames 
      write.table(statsLine, 
                  file = paste0(outDataPath, 
                                filename, phenoName, "_", subSample, "_sample_Adj4", covAdj, "_outRmvd", paste0("_CI0p", substr(ciInterval,3,nchar(ciInterval))), "_statsTable.csv"),
                  col.names = FALSE, row.names = FALSE, append = TRUE, sep = ",")
      
      
    }, error=function(e){cat("WHOOPS :", conditionMessage(e), "\n")})  
    
    
  } # end of inner loop
  
  
} # end of outer loop


# results from the loop above will be saved progressively as they append by line in a single file in the directory specified by user in the outDataPath objetc before the loop
# the names hardcoded here correspond to the set strings for each object before the loop, change correspondingly if change any of the objects besides subSample
corrs_SEM_ML_Kenya <- fread(paste0(outDataPath,"./sim_results_LEC_PTEs_PCL_PTSDsxs_Kenya_sample_Adj4sexAge_outRmvd_CI0p95_statsTable.csv"))
corrs_SEM_ML_Uganda <- fread(paste0(outDataPath,"./sim_results_LEC_PTEs_PCL_PTSDsxs_Uganda_sample_Adj4sexAge_outRmvd_CI0p95_statsTable.csv"))

# results from simulated data will not necessarily be exactly the same as those from the article using real data, however the analysis pipeline is the same for both cases.

colnames(corrs_SEM_ML_Kenya) <- theColNames
colnames(corrs_SEM_ML_Uganda) <- theColNames

corrs_SEM_ML_Kenya
corrs_SEM_ML_Uganda











