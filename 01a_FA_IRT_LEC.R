
# --------------------------------------------------------------------- #
## Factor and IRT analysis, factor scores estimation
## simulated data
## script: Daniel Bustamante, PhD
## Harvard University & Broad Institute of MIT and Harvard
## Dec 2024
# --------------------------------------------------------------------- #

rm(list=ls())

# load libraries
library(data.table)
library(OpenMx)
library(polycor)
library(psych)
library(mirt)
library(matrixStats)
library(latticeExtra)
library(lavaan)

# OpenMx extra steps
omxDetectCores() # number of cores
getOption('mxOptions')$"Number of Threads" # cores used by OpenMx
mxOption(key='Number of Threads', value=parallel::detectCores()) #now
getOption('mxOptions')$"Number of Threads" # cores used by OpenMx
mxOption(key="Default optimizer", value="SLSQP") ## change optimizer to NPSOL, CSOLNP, SLSQP

# load data
theData <- fread("lec_sim_data.csv")
theData <- data.frame(theData) # to bypass data.table() and use regular indexing

## make the ordinal variables ordered factors, is_case variable only for measurement invariance tests
ordFacData <- mxFactor(theData[,c(1:15)], levels = 0:1, ordered = TRUE)#, exclude = 888) # if non-level values
ordFacData4MI <- theData
theData4EFA <- theData[,c(1:15)]
ordFacData4MI[,c(1:15)] <- mxFactor(ordFacData4MI[,c(1:15)], levels = 0:1, ordered = TRUE)#, exclude = 888) # if non-level values

## extract polychoric corrs, SEs etc
corrsOut <- hetcor(ordFacData, ML = FALSE, std.err = TRUE, thresholds = TRUE)
corrs <- corrsOut$correlations
corrsSEs <- corrsOut$std.errors

# simulated thresholds for binary data (15 items)
thresholds <- c(1.834267, 1.823803, 1.084342, 1.5855, 2.123654, 0.6831468, 1.29325, 1.764867, 1.979092, 1.232216, 2.228277, 0.6116753, 1.054758, 0.6451828, 2.1627)

## number of thresholds:
str(ordFacData)
nThresholds <- 1

## quickly check corrs and SEs
round(corrs, 3)
round(corrsSEs, 3)

## obtain eigenvalues and eigenvectors of the correlations
eigenStuff <- eigen(corrs)
eigenVals <- eigenStuff$values
eigenVecs <- eigenStuff$vectors

# how many eigenvalues >1?
above1 <- sort(eigenVals[eigenVals > 1], decreasing = TRUE)
above1
length(above1)

# plor eigenvalues on scree plot
plot(eigenVals, ylab = "Eigenvalues", xlab = "Factors", type = "o", pch = 17, main = "Scree Plot & Kaiser Rule")
abline(a = 1, b = 0, col = "tomato", lwd = 1.7)

# parallel analysis scree plot
paRes <- fa.parallel(ordFacData, fa = "fa", fm = "wls", cor = "poly", n.iter = 100, main = "Parallel Analysis")


# exploratory factor analysis with raw data (using integer data, fa() requires numerica data), 1 to 6-factor models based on number of eigenvalues
efa1facVariRaw <- fa(theData4EFA, nfactors = 1, cor = "poly", fm = "wls", rotate = "varimax") 
efa1facObliRaw <- fa(theData4EFA, nfactors = 1, cor = "poly", fm = "wls", rotate = "oblimin") 

efa2facVariRaw <- fa(theData4EFA, nfactors = 2, cor = "poly", fm = "wls", rotate = "varimax") 
efa2facObliRaw <- fa(theData4EFA, nfactors = 2, cor = "poly", fm = "wls", rotate = "oblimin") 

efa3facVariRaw <- fa(theData4EFA, nfactors = 3, cor = "poly", fm = "wls", rotate = "varimax") 
efa3facObliRaw <- fa(theData4EFA, nfactors = 3, cor = "poly", fm = "wls", rotate = "oblimin") 

efa4facVariRaw <- fa(theData4EFA, nfactors = 4, cor = "poly", fm = "wls", rotate = "varimax") 
efa4facObliRaw <- fa(theData4EFA, nfactors = 4, cor = "poly", fm = "wls", rotate = "oblimin") 

efa5facVariRaw <- fa(theData4EFA, nfactors = 5, cor = "poly", fm = "wls", rotate = "varimax") 
efa5facObliRaw <- fa(theData4EFA, nfactors = 5, cor = "poly", fm = "wls", rotate = "oblimin") 

efa6facVariRaw <- fa(theData4EFA, nfactors = 6, cor = "poly", fm = "wls", rotate = "varimax") 
efa6facObliRaw <- fa(theData4EFA, nfactors = 6, cor = "poly", fm = "wls", rotate = "oblimin")


# check loadings and R^2s
# orthogonal
print(efa1facVariRaw$loadings, cutoff =.2)
print(efa2facVariRaw$loadings, cutoff =.2)
print(efa3facVariRaw$loadings, cutoff =.2)
print(efa4facVariRaw$loadings, cutoff =.2)
print(efa5facVariRaw$loadings, cutoff =.2)
print(efa6facVariRaw$loadings, cutoff =.2)

# oblique
print(efa1facObliRaw$loadings, cutoff =.2)
print(efa2facObliRaw$loadings, cutoff =.2)
print(efa3facObliRaw$loadings, cutoff =.2)
print(efa4facObliRaw$loadings, cutoff =.2)
print(efa5facObliRaw$loadings, cutoff =.2)
print(efa6facObliRaw$loadings, cutoff =.2)

# check figures from each EFA
# orthogonal
fa.diagram(efa1facVariRaw, cut = .1)
fa.diagram(efa2facVariRaw, cut = .1)
fa.diagram(efa3facVariRaw, cut = .1)
fa.diagram(efa4facVariRaw, cut = .1)
fa.diagram(efa5facVariRaw, cut = .1)
fa.diagram(efa6facVariRaw, cut = .1)
fa.diagram(efa7facVariRaw, cut = .1)
# oblique
fa.diagram(efa1facObliRaw, cut = .1)
fa.diagram(efa2facObliRaw, cut = .1)
fa.diagram(efa3facObliRaw, cut = .1)
fa.diagram(efa4facObliRaw, cut = .1)
fa.diagram(efa5facObliRaw, cut = .1)
fa.diagram(efa6facObliRaw, cut = .1)
fa.diagram(efa7facObliRaw, cut = .1)



# confirmatory factor analysis
# 4-factor oblique solution (see item separation below)
# OpenMx
numFactors <- 4
numVars <- ncol(ordFacData)
vars <- colnames(ordFacData)
rot <- "oblique"

items1fac <- vars[(1:5)] # 5 non-death accident items, 
items2fac <- vars[c(6:9)] # 4 physical & sexual violence items, 
items3fac <- vars[c(10:12)] # 3 for the third factor with 2 combat & captivity together with, life threatening and severe human suffering
items4fac <- vars[c(13:15)] # 4 for the 4th factor with Serious injury, harm, or death you caused to someone else item, other stressful event, together with the two witnessing death items

load1fac <- vars %in% items1fac
load2fac <- vars %in% items2fac
load3fac <- vars %in% items3fac
load4fac <- vars %in% items4fac

# free and fix the loading paths accordingly
freeFixed <- c(load1fac, load2fac, load3fac, load4fac)

# check them
length(freeFixed[freeFixed == TRUE])
length(freeFixed[freeFixed == FALSE])

# set up OpenMx objects for loadings, errors, thresholds, factor correlations, etc
load              <- mxMatrix("Full", numVars, numFactors, values=0, free=freeFixed, lbound = 0.00001, ubound = 0.9999, name="load")
Errors            <- mxAlgebra(vec2diag(Ones - diag2vec(load %*% t(load))), name="Errors")
Ones              <- mxMatrix("Unit", numVars, 1, name="Ones")
Mean              <- mxMatrix("Full", 1, numVars, free = F, name="Mean")
impliedCovs       <- mxAlgebra((load %*% FacCov %*% t(load)) + Errors, name="impliedCovs")
threshDev         <- mxMatrix("Full", name="threshDev", nrow=nThresholds, ncol=numVars, values=thresholds, free = T, lbound = rep( c(-Inf,rep(.01,(nThresholds-1))) , numVars), dimnames =  list(c(), vars))
unitLower         <- mxMatrix("Lower",nThresholds,nThresholds,values=1,free=F,name="unitLower")
thresholdMatrix   <- mxAlgebra(unitLower %*% threshDev, name="thresholdMatrix")
FacCov            <- mxMatrix("Symm", numFactors, numFactors, values=c(1,0,0,0,1,0,0,1,0,1), free=c(F,T,T,T,F,T,T,F,T,F), lbound = -0.9999, ubound = 0.99999, name="FacCov")

obj <-     mxExpectationNormal(covariance="impliedCovs", means="Mean", dimnames = vars, thresholds="thresholdMatrix")
fitFun <- mxFitFunctionWLS(type = "DWLS")
ordData <- mxData(observed=ordFacData, type='raw')
FullModel4fac <- mxModel("Full4fac", load, Ones, Mean, Errors, impliedCovs, threshDev, unitLower, thresholdMatrix, FacCov, obj, fitFun, ordData)
fitCFA4fac <- mxTryHardOrdinal(FullModel4fac)
summary(fitCFA4fac) 

checkIdentfitCFA4fac <- mxCheckIdentification(fitCFA4fac)

# lavaan
lvModelCFA4obl <- c(paste0("f1 =~ ", paste(items1fac, collapse = " + ")),
                          paste0("f2 =~ ", paste(items2fac, collapse = " + ")),
                          paste0("f3 =~ ", paste(items3fac, collapse = " + ")),
                          paste0("f4 =~ ", paste(items4fac, collapse = " + ")),
                          "f1 ~~ f2", 
                          "f1 ~~ 0.7*f3", # fixed the correlations of these factors of the simulated data covariance matrix since it was non positive definite
                          "f1 ~~ f4",
                          "f2 ~~ f3",
                          "f2 ~~ f4",
                          "f3 ~~ f4"
)
lvFitCFA4obl <- cfa(lvModelCFA4obl, data = ordFacData, std.lv=TRUE, missing = "pairwise")#, auto.cov.lv.x=FALSE)
summlvFitCFA4obl <- summary(lvFitCFA4obl, fit.measures = TRUE, standardized = TRUE)
summlvFitCFA4obl


# Measurement invariance
# model from CFA results
lvModelCFA4obl_MI <- lvModelCFA4obl

# configural
lvFit_config <- cfa(lvModelCFA4obl_MI, data = ordFacData4MI, std.lv = TRUE, missing = "pairwise", group = "is_case")
summary(lvFit_config, fit.measures = TRUE, standardized = TRUE)
fitMeasures(lvFit_config, c("rmsea", "cfi", "tli"))

lavInspect(lvFit_config, "cov.lv")
eigen(lavInspect(lvFit_config, "cov.lv")$`0`)
eigen(lavInspect(lvFit_config, "cov.lv")$`1`)

# metric
lvFit_metric <- cfa(lvModelCFA4obl_MI, data = ordFacData4MI, std.lv=TRUE, missing = "pairwise", group = "is_case", group.equal = "loadings")
summary(lvFit_metric, fit.measures = TRUE, standardized = TRUE)
fitMeasures(lvFit_metric, c("rmsea", "cfi", "tli"))

# scalar
lvFit_scalar <- cfa(lvModelCFA4obl_MI, data = ordFacData4MI, std.lv=TRUE, missing = "pairwise", group = "is_case", group.equal = c("loadings", "intercepts"))
summary(lvFit_scalar, fit.measures = TRUE, standardized = TRUE)
fitMeasures(lvFit_scalar, c("rmsea", "cfi", "tli"))

# strict (constraining residuals too in addition to loadings and intercepts)
lvFit_strict <- cfa(lvModelCFA4obl_MI, data = ordFacData4MI, std.lv=TRUE, missing = "pairwise", group = "is_case", group.equal = c("loadings", "intercepts", "residuals"))
summary(lvFit_strict, fit.measures = TRUE, standardized = TRUE)
fitMeasures(lvFit_strict, c("rmsea", "cfi", "tli"))

# produce a table for MI differences
MI_table <- data.frame(t(sapply(list(configural = lvFit_config, metric = lvFit_metric, scalar = lvFit_scalar, strict = lvFit_strict), fitMeasures, fit.measures = c("rmsea", "cfi", "tli"))))
MI_table$delta_rmsea <- c(NA, diff(MI_table[,1]))
MI_table$delta_cfi <- c(NA, diff(MI_table[,2]))
MI_table$delta_tli <- c(NA, diff(MI_table[,3]))
MI_table


# Item response theory
# multidimensional case binary data
# specify model (from CFA results)

length(unique(c(items1fac, items2fac, items3fac, items4fac)))

modelIRT <- '
F1 = lec_1, lec_2, lec_3, lec_4, lec_5
F2 = lec_6, lec_7, lec_8, lec_9
F3 = lec_10, lec_11, lec_12
F4 = lec_13, lec_14, lec_q15
COV = F1*F2, F1*F3, F1*F4, F2*F3, F2*F4, F3*F4
'
modelIRT
# fit the model
# best fitting model from CFA
modelIRTfit <- mirt(theData4EFA, modelIRT, itemtype = "graded", method = "MHRM") # default EM becomes unstable here with three or more factors
# method = a character object specifying the estimation algorithm to be used. The default is 'EM', for the standard EM algorithm with fixed quadrature, 'QMCEM' for quasi-Monte Carlo EM estimation, or 'MCEM' for Monte Carlo EM estimation. The option 'MHRM' may also be passed to use the MH-RM algorithm, 'SEM' for the Stochastic EM algorithm (first two stages of the MH-RM stage using an optimizer other than a single Newton-Raphson iteration), and 'BL' for the Bock and Lieberman approach (generally not recommended for longer tests).
# The 'EM' is generally effective with 1-3 factors, but methods such as the 'QMCEM', 'MCEM', 'SEM', or 'MHRM' should be used when the dimensions are 3 or more. Note that when the optimizer is stochastic the associated SE.type is automatically changed to SE.type = 'MHRM' by default to avoid the use of quadrature
summary(modelIRTfit)
fit_stats <- M2(modelIRTfit, type = "C2", CI = 0.95, QMC = TRUE)
fit_stats

item_pars <- coef(modelIRTfit, IRTpars=TRUE , simplify=TRUE)
item_pars # a are discrimination for each item-factor, b1 is difficulty
item_params_table <- data.frame(item_pars$items)
item_params_table$items <- rownames(item_params_table)
item_params_table <- item_params_table[,c(ncol(item_params_table),1:(ncol(item_params_table)-1))]
rownames(item_params_table) <- NULL
item_params_table


# estimate factor scores
lvcfaFacFSsEBM <- lavPredict(lvFitCFA4obl, method = "EBM") # EBM for empirical Bayes modal 
describe(lvcfaFacFSsEBM)








