# Jags-Ymet-XmetMulti-Mrobust.R 
# Accompanies the book:
#  Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#  A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

#===============================================================================
# Setup and data
#===============================================================================

library(coda)
library(LaplacesDemon)
library(ggplot2)

setwd("/Users/baste/Documents/GitHub/math2269-bayesian-workforce-analysis")
source("data/DBDA2E-utilities.R")

# data used to fit bayesian model
myData <- read.csv("data/jags_model_data.csv", check.names=FALSE)

# Data used to identify occupations for posterior prediction
modelData <- read.csv("data/model_data.csv", check.names = FALSE)

yName <- "employment_growth_rate"

xName <- c("log_employed", "median_age",	"higher_education_share",	"female_share",	"part_time_share",	"average_full_time_hours")

#select random ANZSCO codes for prediction segment
predictionCodes <- c("2611", "3513", "2411", "3312", "2613")

# keep the selected occupations for later
predictionData <- modelData[modelData$anzsco_code %in% predictionCodes,]
predictionData <- predictionData[match(predictionCodes, predictionData$anzsco_code),] #fix orders

predictionData[, c("anzsco_code", "occupation")] #check

#extract predictors , used by JAGS
xPred <- as.matrix(predictionData[, xName])

# Create lookup so pred[1], pred[2], etc. can be matched back to their ANZSCO code and occupation name
predictionLookup <- predictionData[, c("anzsco_code", "occupation")]
predictionLookup$pred_parameter <- paste0("pred[", seq_len(nrow(predictionLookup)),"]")


#===============================================================================

genMCMC = function( data , xName="x" , yName="y" , 
                    numSavedSteps=10000 , thinSteps=1 , saveName=NULL  ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault ,
                    xPred) { 
  require(runjags)
  #-----------------------------------------------------------------------------
  # THE DATA.
  y = data[,yName]
  x = as.matrix(data[,xName],ncol=length(xName))
  # Do some checking that data make sense:
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
  if ( any( !is.finite(x) ) ) { stop("All x values must be finite.") }
  cat("\nCORRELATION MATRIX OF PREDICTORS:\n ")
  show( round(cor(x),3) )
  cat("\n")
  flush.console()
  # make sure prediction data are a matrix
  xPred = as.matrix(xPred)
  # Check prediction data
  if (any(!is.finite(xPred))) {stop("All prediction values must be finite.")}
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    x = x ,
    y = y ,
    Nx = dim(x)[2] ,
    Ntotal = dim(x)[1] ,
    xPred = xPred ,
    Npred = nrow(xPred)
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString <- "

data {

  # Standardize response
  ym <- mean(y)
  ysd <- sd(y)

  for (i in 1:Ntotal) {
  zy[i] <- (y[i] - ym) / ysd
  }

  # Standardize predictors
  for (j in 1:Nx) {
    xm[j] <- mean(x[,j])
    xsd[j] <- sd(x[,j])

    for (i in 1:Ntotal) {
    zx[i,j] <- (x[i,j] - xm[j]) / xsd[j]
    }
    }
}

model {
  # Student-t regression likelihood
  for (i in 1:Ntotal) {
    mu[i] <- zbeta0 +sum(zbeta[1:Nx] * zx[i,1:Nx])
    zy[i] ~ dt(mu[i], 1 / zsigma^2, nu)
  }


  # Weak prior for standardized intercept
  zbeta0 ~ dnorm(0, 1/2^2)
  
  # Weak priors for standardized regression coefficients
  for (j in 1:Nx) {
    zbeta[j] ~ dnorm(0, 1/2^2)
  }
  
  # Broad prior for standardized residual scale
  zsigma ~ dunif(1.0E-5, 1.0E+1)
  
  # Prior for Student-t degrees of freedom
  nuMinusOne ~ dexp(1/29.0)
  nu <- nuMinusOne + 1


  # Transform parameters back to original scale
  beta[1:Nx] <- (zbeta[1:Nx] / xsd[1:Nx]) * ysd
  beta0 <- zbeta0 * ysd + ym - sum(zbeta[1:Nx] * xm[1:Nx] / xsd[1:Nx]) * ysd
  sigma <- zsigma * ysd
  
  # Posterior predictions
  for (i in 1:Npred) {
    pred[i] <- beta0 + sum(beta[1:Nx] * xPred[i,1:Nx])}
}
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Let JAGS do it...
  
  # Must standardize data first...
  #   lmInfo = lm( zy ~ zx )
  #   initsList = list(
  #     beta0 = lmInfo$coef[1] ,   
  #     beta = lmInfo$coef[-1] ,        
  #     sigma = sqrt(mean(lmInfo$resid^2)) ,
  #     nu = 5
  #   )
  
  
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c("beta0", "beta", "sigma", "zbeta0", "zbeta", "zsigma", "nu", "pred")
  adaptSteps = 500  # Number of steps to "tune" the samplers
  burnInSteps = 1000
  runJagsOut <- run.jags( method=runjagsMethod ,
                          model="TEMPmodel.txt" , 
                          monitor=parameters , 
                          data=dataList ,  
                          #inits=initsList , 
                          n.chains=nChains ,
                          adapt=adaptSteps ,
                          burnin=burnInSteps , 
                          sample=ceiling(numSavedSteps/nChains) ,
                          thin=thinSteps ,
                          summarise=FALSE ,
                          plots=FALSE )
  codaSamples = as.mcmc.list( runJagsOut )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
} # end function

#===============================================================================

smryMCMC = function(  codaSamples , 
                      saveName=NULL ) {
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  paramName = colnames(mcmcMat)
  for ( pName in paramName ) {
    summaryInfo = rbind( summaryInfo , summarizePost( mcmcMat[,pName] ) )
  }
  rownames(summaryInfo) = paramName
  summaryInfo = rbind( summaryInfo , 
                       "log10(nu)" = summarizePost( log10(mcmcMat[,"nu"]) ) )
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================

plotMCMC = function( codaSamples , data , xName="x" , yName="y" ,
                     showCurve=FALSE ,  pairsPlot=FALSE ,
                     saveName=NULL , saveType="jpg" ) {
  # showCurve is TRUE or FALSE and indicates whether the posterior should
  #   be displayed as a histogram (by default) or by an approximate curve.
  # pairsPlot is TRUE or FALSE and indicates whether scatterplots of pairs
  #   of parameters should be displayed.
  #-----------------------------------------------------------------------------
  y = data[,yName]
  x = as.matrix(data[,xName])
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  zbeta0 = mcmcMat[,"zbeta0"]
  zbeta  = mcmcMat[,grep("^zbeta$|^zbeta\\[",colnames(mcmcMat))]
  if ( ncol(x)==1 ) { zbeta = matrix( zbeta , ncol=1 ) }
  zsigma = mcmcMat[,"zsigma"]
  beta0 = mcmcMat[,"beta0"]
  beta  = mcmcMat[,grep("^beta$|^beta\\[",colnames(mcmcMat))]
  if ( ncol(x)==1 ) { beta = matrix( beta , ncol=1 ) }
  sigma = mcmcMat[,"sigma"]
  nu = mcmcMat[,"nu"]
  log10nu = log10(nu)
  #-----------------------------------------------------------------------------
  # Compute R^2 for credible parameters:
  YcorX = cor( y , x ) # correlation of y with each x predictor
  Rsq = as.vector(zbeta %*% matrix(YcorX, ncol = 1)) # convert matrix to numeric vector
  #-----------------------------------------------------------------------------
  if ( pairsPlot ) {
    # Plot the parameters pairwise, to see correlations:
    openGraph()
    nPtToPlot = 1000
    plotIdx = floor(seq(1,chainLength,by=chainLength/nPtToPlot))
    panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...) {
      oldusr <- par("usr")
      on.exit(par(usr = oldusr))
      par(usr = c(0, 1, 0, 1)) # temporarily change coords to fix positioning
      r = (cor(x, y))
      txt = format(c(r, 0.123456789), digits=digits)[1]
      txt = paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      text(0.5, 0.5, txt, cex=1.25 ) # was cex=cex.cor*r
    }
    pairs( cbind( beta0 , beta , sigma , log10nu )[plotIdx,] ,
           labels=c( "beta[0]" , 
                     paste0("beta[",1:ncol(beta),"]\n",xName) , 
                     expression(sigma) ,  expression(log10(nu)) ) , 
           lower.panel=panel.cor , col="skyblue" )
    if ( !is.null(saveName) ) {
      saveGraph( file=paste(saveName,"PostPairs",sep=""), type=saveType)
    }
  }
  #-----------------------------------------------------------------------------
  # Marginal histograms:
  
  decideOpenGraph = function( panelCount , saveName , finished=FALSE , 
                              nRow=2 , nCol=3 ) {
    # If finishing a set:
    if ( finished==TRUE ) {
      if ( !is.null(saveName) ) {
        saveGraph( file=paste0(saveName,ceiling((panelCount-1)/(nRow*nCol))), 
                   type=saveType)
      }
      panelCount = 1 # re-set panelCount
      return(panelCount)
    } else {
      # If this is first panel of a graph:
      if ( ( panelCount %% (nRow*nCol) ) == 1 ) {
        # If previous graph was open, save previous one:
        if ( panelCount>1 & !is.null(saveName) ) {
          saveGraph( file=paste0(saveName,(panelCount%/%(nRow*nCol))), 
                     type=saveType)
        }
        # Open new graph
        openGraph(width=nCol*7.0/3,height=nRow*2.0)
        layout( matrix( 1:(nRow*nCol) , nrow=nRow, byrow=TRUE ) )
        par( mar=c(4,4,2.5,0.5) , mgp=c(2.5,0.7,0) )
      }
      # Increment and return panel count:
      panelCount = panelCount+1
      return(panelCount)
    }
  }
  
  # Original scale:
  panelCount = 1
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( beta0 , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(beta[0]) , main="Intercept" )
  for ( bIdx in 1:ncol(beta) ) {
    panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
    histInfo = plotPost( beta[,bIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                         xlab=bquote(beta[.(bIdx)]) , main=xName[bIdx] )
  }
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( sigma , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(sigma) , main=paste("Scale") )
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( log10nu , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(log10(nu)) , main=paste("Normality") )
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( Rsq , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(R^2) , main=paste("Prop Var Accntd") )
  panelCount = decideOpenGraph( panelCount , finished=TRUE , saveName=paste0(saveName,"PostMarg") )
  
  # Standardized scale:
  panelCount = 1
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMargZ") )
  histInfo = plotPost( zbeta0 , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(z*beta[0]) , main="Intercept" )
  for ( bIdx in 1:ncol(beta) ) {
    panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMargZ") )
    histInfo = plotPost( zbeta[,bIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                         xlab=bquote(z*beta[.(bIdx)]) , main=xName[bIdx] )
  }
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMargZ") )
  histInfo = plotPost( zsigma , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(z*sigma) , main=paste("Scale") )
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMargZ") )
  histInfo = plotPost( log10nu , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(log10(nu)) , main=paste("Normality") )
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMargZ") )
  histInfo = plotPost( Rsq , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(R^2) , main=paste("Prop Var Accntd") )
  panelCount = decideOpenGraph( panelCount , finished=TRUE , saveName=paste0(saveName,"PostMargZ") )
  
  #-----------------------------------------------------------------------------
}
#===============================================================================
# MCMC run 1
#===============================================================================

# create new folder for the run
dir.create("runs/run1", recursive = TRUE, showWarnings = FALSE)

numSavedSteps <- 10000
thinSteps <- 1
nChains <- 3
# adaptSteps = 500
# burnInSteps = 1000

# file name
fileNameRoot <- "runs/run1/Workforce-Run1-"

startTime <- proc.time()

mcmcCoda <- genMCMC(
  data = myData,
  xName = xName,
  yName = yName,
  numSavedSteps = numSavedSteps,
  thinSteps = thinSteps,
  nChains = nChains,
  xPred = xPred,
  saveName = fileNameRoot
)

#duration of the run
stopTime <- proc.time()
duration <- stopTime - startTime
show(duration)

# summarise posterior distributions
summaryInfo <- smryMCMC(mcmcCoda, saveName = fileNameRoot)

# plot MCMCs
plotMCMC(mcmcCoda, data = myData, xName = xName, yName = yName, showCurve = FALSE, pairsPlot = TRUE, saveName = fileNameRoot)

# specify parameters to assess with MCMC diagnostics
diagParams <- c("beta0", paste0("beta[", 1:6, "]"), "sigma", "nu")

for (p in diagParams) {
  diagMCMC(mcmcCoda, parName = p, saveName = fileNameRoot)
}

#===============================================================================
# MCMC run 2
#===============================================================================

# create new folder for the run
dir.create("runs/run2", recursive = TRUE, showWarnings = FALSE)

numSavedSteps <- 30000
thinSteps <- 1
nChains <- 3
# adaptSteps = 500
# burnInSteps = 1000

# file name
fileNameRoot <- "runs/run2/Workforce-Run2-"

startTime <- proc.time()

mcmcCoda_run2 <- genMCMC(
  data = myData,
  xName = xName,
  yName = yName,
  numSavedSteps = numSavedSteps,
  thinSteps = thinSteps,
  nChains = nChains,
  xPred = xPred,
  saveName = fileNameRoot
)

#duration of the run
stopTime <- proc.time()
duration_run2 <- stopTime - startTime
show(duration_run2)

# summarise posterior distributions
summaryInfo_run2 <- smryMCMC(mcmcCoda_run2, saveName = fileNameRoot)

# plot MCMCs
plotMCMC(mcmcCoda_run2, data = myData, xName = xName, yName = yName, showCurve = FALSE, pairsPlot = TRUE, saveName = fileNameRoot)

# specify parameters to assess with MCMC diagnostics
diagParams <- c("beta0", paste0("beta[", 1:6, "]"), "sigma", "nu")

for (p in diagParams) {
  diagMCMC(mcmcCoda_run2, parName = p, saveName = fileNameRoot)
}

#===============================================================================
# MCMC run 3
#===============================================================================

# create new folder for the run
dir.create("runs/run3", recursive = TRUE, showWarnings = FALSE)

numSavedSteps <- 60000
thinSteps <- 1
nChains <- 3
#daptSteps = 1500 
#burnInSteps = 5000

# file name
fileNameRoot <- "runs/run3/Workforce-Run3-"

startTime <- proc.time()

mcmcCoda_run3 <- genMCMC(
  data = myData,
  xName = xName,
  yName = yName,
  numSavedSteps = numSavedSteps,
  thinSteps = thinSteps,
  nChains = nChains,
  xPred = xPred,
  saveName = fileNameRoot
)

#duration of the run
stopTime <- proc.time()
duration_run3 <- stopTime - startTime
show(duration_run3)

# summarise posterior distributions
summaryInfo_run3 <- smryMCMC(mcmcCoda_run3, saveName = fileNameRoot)

# plot MCMCs
plotMCMC(mcmcCoda_run3, data = myData, xName = xName, yName = yName, showCurve = FALSE, pairsPlot = TRUE, saveName = fileNameRoot)

# specify parameters to assess with MCMC diagnostics
diagParams <- c("beta0", paste0("beta[", 1:6, "]"), "sigma", "nu")

for (p in diagParams) {
  diagMCMC(mcmcCoda_run3, parName = p, saveName = fileNameRoot)
}

# standardised parameters too
diagParams_std <- c("zbeta0", paste0("zbeta[", 1:6, "]"), "zsigma")

for (p in diagParams_std) {
  diagMCMC(
    mcmcCoda_run3,
    parName = p,
    saveName = fileNameRoot
  )
}

# Gelman-Rubin shrink factors for final Run 3
gelman_results <- gelman.diag(
  mcmcCoda_run3[, diagParams],
  autoburnin = FALSE,
  multivariate = FALSE)

gelman_results$psrf

#===============================================================================
# Prior sensitivity analysis
# Stronger coefficient prior: zbeta ~ Normal(0, 1)
#===============================================================================

# create new folder for the run
dir.create("runs/sensitivity_var1", recursive = TRUE, showWarnings = FALSE)

numSavedSteps <- 60000
thinSteps <- 1
nChains <- 3

# file name
fileNameRoot <- "runs/sensitivity_var1/Workforce-Sensitivity-Var1-"

startTime <- proc.time()

mcmcCoda_sensitivity <- genMCMC(
  data = myData,
  xName = xName,
  yName = yName,
  numSavedSteps = numSavedSteps,
  thinSteps = thinSteps,
  nChains = nChains,
  xPred = xPred,
  saveName = fileNameRoot
)

#duration of the run
stopTime <- proc.time()
duration_sensitivity <- stopTime - startTime
show(duration_sensitivity)

# summarise posterior distributions
summaryInfo_sensitivity <- smryMCMC(mcmcCoda_sensitivity, saveName = fileNameRoot)

# plot MCMCs
plotMCMC(mcmcCoda_sensitivity, data = myData, xName = xName, yName = yName, showCurve = FALSE, pairsPlot = TRUE, saveName = fileNameRoot)

# specify parameters to assess with MCMC diagnostics
diagParams <- c("beta0", paste0("beta[", 1:6, "]"), "sigma", "nu")

for (p in diagParams) {
  diagMCMC(mcmcCoda_sensitivity, parName = p, saveName = fileNameRoot)
}


#===============================================================================
# Posterior predictions for selected ANZSCO occupations
#===============================================================================

# extract summaries of posteriors for pred[1] to pred[5]
predSummary_run3 <- summaryInfo_run3[grep("^pred\\[", rownames(summaryInfo_run3)),]
predSummary_run3

# combine summaries with the anzsco names and codes
predictionResults <- cbind(
  predictionLookup,
  predSummary_run3[, c("Mean",  "Median", "Mode", "HDIlow", "HDIhigh")])
predictionResults

# Names of prediction parameters saved from JAGS
predNames <- paste0("pred[", 1:5, "]")

# Convert Run 3 MCMC samples to a matrix
mcmcMat_run3 <- as.matrix(mcmcCoda_run3)

# Open figure
openGraph(width = 12, height = 7)

# Arrange five plots in a 2 x 3 layout
par(mfrow = c(2, 3))

# Plot posterior distribution for each selected occupation
for (i in 1:5) {
  
  plotPost(
    mcmcMat_run3[, predNames[i]],
    main = paste0(predictionLookup$anzsco_code[i], " - ", predictionLookup$occupation[i]),
    xlab = "Expected employment growth rate (%)",
    showCurve = FALSE
  )
}

plot.new()

# Save figure
saveGraph(file = "runs/run3/Workforce-Run3-PredictionPosteriors", type = "jpg")

#===============================================================================
# Posterior predictive check
#===============================================================================

# Extract posterior modes from final Run 3
coefficients <- summaryInfo_run3[c("beta0", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", "beta[6]"), "Mode"]
sigma_mode <- summaryInfo_run3["sigma", "Mode"]
nu_mode    <- summaryInfo_run3["nu", "Mode"]

# Predictor and observed outcome data
x <- myData[, xName]
y <- myData[, yName]

# Expected growth for each observed occupation
mu_hat <- as.matrix(
  cbind(Intercept = 1, x)
) %*% as.vector(coefficients)

# Generate replicated growth values from fitted Student-t model
set.seed(123)

randomData <- rstp(
  n = nrow(myData),
  mu = as.vector(mu_hat),
  tau = 1 / sigma_mode^2,
  nu = nu_mode
)

# Combine observed and replicated data
observed <- data.frame(growth = y, type = "Observed")
predicted <- data.frame(growth = randomData, type = "Posterior predictive")
ppcData <- rbind(observed, predicted)

# Density comparison
ppc_plot <- ggplot(ppcData, aes(x = growth, fill = type)) +
  geom_density(alpha = 0.3) +
  labs(
    x = "Annual employment growth rate (%)",
    y = "Density",
    fill = NULL,
    title = "Posterior Predictive Check") +
  theme_minimal()

ppc_plot

# save plot
ggsave("runs/run3/Workforce-Run3-PosteriorPredictiveCheck.png", ppc_plot, width = 8, height = 5, dpi = 300)
