# JAGS model based on the following models:
#
# Jags-Ydich-XmetMulti-Mlogistic.R 
# Jags-Ymet-XmetSsubj-MrobustHier.R 
# 
# Accompanies the book:
#   Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#   A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

# source("DBDA2E-utilities.R")

#===============================================================================

genMCMC2Heirarchy = function( data , xName="x", yName="y", sName="s", 
                    numSavedSteps=10000 , thinSteps=1 , saveName=NULL ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault ) { 
  require(runjags)
  
  stopifnot(length(sName) == 2)
  
  #-----------------------------------------------------------------------------
  # THE DATA.
  y = data[,yName]
  x = as.matrix(data[,xName],ncol=length(xName))
  s1 = as.numeric(data[,sName[1]])
  s2 = as.numeric(data[,sName[2]])
  
  # Do some checking that data make sense:
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
  if ( any( !is.finite(x) ) ) { stop("All x values must be finite.") }
  cat("\nCORRELATION MATRIX OF PREDICTORS:\n ")
  # show( round(cor(x),3) )
  cat("\n")
  flush.console()

  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    x = x ,
    y = y ,
    Nx = dim(x)[2] ,
    Ntotal = dim(x)[1],
    # Heirarchical info
    s1 = s1,
    s2 = s2, 
    Ngroup1 = max(s1),
    Ngroup2 = max(s2)
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  # Standardize the data:
  data {
    for ( j in 1:Nx ) {
      xm[j]  <- mean(x[,j])
      xsd[j] <-   sd(x[,j])
      for ( i in 1:Ntotal ) {
        zx[i,j] <- ( x[i,j] - xm[j] ) / xsd[j]
      }
    }
  }
  # Specify the model for standardized data:
  model {
    for ( i in 1:Ntotal ) {
      # In JAGS, ilogit is logistic:
      y[i] ~ dbern(ilogit(zbeta0[s2[i]] + sum( zbeta[1:Nx, s2[i]] * zx[i,1:Nx])))
    }

    # Distribution for betas based on group 2
    for (n in 1:Ngroup2) {
      for ( j in 1:Nx ) {
        zbeta[j, n] ~ dnorm(zbetagrp1mu[s1[n]], 1/(zbetagrp1sigma)^2)
      }
      zbeta0[n] ~ dnorm(zbeta0grp1mu[s1[n]], 1/(zbeta0grp1sigma)^2)
    }

    # Distribution for group 1
    for (l in 1:Ngroup1) {
      zbetagrp1mu[l] ~ dnorm(zbetamu, 1/(zbetasigma)^2)
      zbeta0grp1mu[l] ~ dnorm(zbeta0mu, 1/(zbeta0sigma)^2)
    }

    # Priors vague on standardized scale:
    zbeta0mu ~ dnorm(0, 1/2^2)
    zbetamu ~ dnorm(0, 1/2^2)
  
    # Variance prior for all levels of heirarchy
    zbeta0sigma ~ dunif( 1.0E-3 , 1.0E+3 )
    zbetasigma ~ dunif( 1.0E-3 , 1.0E+3 )
    zbeta0grp1sigma ~ dunif( 1.0E-3 , 1.0E+3 )
    zbetagrp1sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  
    # Transform to original scale:
    for (n in 1:Ngroup2) {
      beta[1:Nx, n] <- zbeta[1:Nx, n] / xsd[1:Nx]
      beta0[n] <- zbeta0[n] - sum( zbeta[1:Nx, n] * xm[1:Nx] / xsd[1:Nx] )
    }
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="TEMPmodel.txt" )
  
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Let JAGS do it...
  
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c("beta0", "beta", 'zbeta0mu', 'zbetamu', 'zbeta0grp1mu',  
                 'zbetagrp1mu', 'zbeta0sigma','zbetasigma',
                 'zbeta0grp1sigma', 'zbetagrp1sigma', "zbeta0" , "zbeta" )
  adaptSteps = 1000  # Number of steps to "tune" the samplers
  burnInSteps = 1500
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

genMCMC1Heirarchy = function( data , xName="x", yName="y", sName="s", 
                              numSavedSteps=10000 , thinSteps=1 , saveName=NULL ,
                              runjagsMethod=runjagsMethodDefault , 
                              nChains=nChainsDefault ) { 
  require(runjags)
  
  stopifnot(length(sName) == 1)
  
  #-----------------------------------------------------------------------------
  # THE DATA.
  y = data[,yName]
  x = as.matrix(data[,xName],ncol=length(xName))
  s2 = as.numeric(factor(data[,sName]))

  # Do some checking that data make sense:
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
  if ( any( !is.finite(x) ) ) { stop("All x values must be finite.") }
  cat("\nCORRELATION MATRIX OF PREDICTORS:\n ")
  # show( round(cor(x),3) )
  cat("\n")
  flush.console()
  
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    x = x ,
    y = y ,
    Nx = dim(x)[2] ,
    Ntotal = dim(x)[1],
    # Heirarchical info
    s2 = s2, 
    Ngroup2 = max(s2)
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  # Standardize the data:
  data {
  for ( j in 1:Nx ) {
  xm[j]  <- mean(x[,j])
  xsd[j] <-   sd(x[,j])
  for ( i in 1:Ntotal ) {
  zx[i,j] <- ( x[i,j] - xm[j] ) / xsd[j]
  }
  }
  }
  # Specify the model for standardized data:
  model {
  for ( i in 1:Ntotal ) {
  # In JAGS, ilogit is logistic:
  y[i] ~ dbern(ilogit(zbeta0[s2[i]] + sum( zbeta[1:Nx, s2[i]] * zx[i,1:Nx])))
  }
  
  # Distribution for betas based on group 2
  for (n in 1:Ngroup2) {
    for ( j in 1:Nx ) {
      zbeta[j, n] ~ dnorm(zbetamu, 1/(zbetasigma)^2)
    }
    zbeta0[n] ~ dnorm(zbeta0mu, 1/(zbeta0sigma)^2)
  }
  
  # Priors vague on standardized scale:
  zbeta0mu ~ dnorm(0, 1/2^2)
  zbetamu ~ dnorm(0, 1/2^2)
  
  # Variance prior for all levels of heirarchy
  zbeta0sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbetasigma ~ dunif( 1.0E-3 , 1.0E+3 )

  # Transform to original scale:
  for (n in 1:Ngroup2) {
    beta[1:Nx, n] <- zbeta[1:Nx, n] / xsd[1:Nx]
    beta0[n] <- zbeta0[n] - sum( zbeta[1:Nx, n] * xm[1:Nx] / xsd[1:Nx] )
    }
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="TEMPmodel.txt" )
  
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Let JAGS do it...
  
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c("beta0", "beta", 'zbeta0mu', 'zbetamu',  
                 'zbeta0sigma','zbetasigma', "zbeta0" , "zbeta" )
  adaptSteps = 1000  # Number of steps to "tune" the samplers
  burnInSteps = 1500
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
  mcmcMat = as.matrix(codaSamples)
  paramName = colnames(mcmcMat)
  for ( pName in paramName ) {
    summaryInfo = rbind( summaryInfo , summarizePost( mcmcMat[,pName] ) )
  }
  rownames(summaryInfo) = paramName
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================

plotMCMC = function(codaSamples, data, xName = "x", yName = "y", 
                              sName = "s",showCurve = FALSE, 
                              compVal = 0.5, rope = NULL, ropeDiff = NULL,
                              pairsPlot = FALSE, compValDiff = 0.0,
                              saveName = NULL, saveType = "jpg") {
  # showCurve is TRUE or FALSE and indicates whether the posterior should
  #   be displayed as a histogram (by default) or by an approximate curve.
  # pairsPlot is TRUE or FALSE and indicates whether scatterplots of pairs
  #   of parameters should be displayed.
  #-----------------------------------------------------------------------------
  y = data[,yName]
  x = data[,xName]
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  zbeta0 = mcmcMat[,grep("^zbeta0$|^zbeta0\\[",colnames(mcmcMat))]
  zbeta  = mcmcMat[,grep("^zbeta$|^zbeta\\[",colnames(mcmcMat))]
  if ( ncol(x)==1 ) { zbeta = matrix( zbeta , ncol=1 ) }
  beta0 = mcmcMat[,grep("^beta0$|^beta0\\[",colnames(mcmcMat))]
  beta  = mcmcMat[,grep("^beta$|^beta\\[",colnames(mcmcMat))]
  if ( ncol(x)==1 ) { beta = matrix( beta , ncol=1 ) }
  zbeta0mu <- mcmcMat[,"zbeta0mu"]
  zbetamu <- mcmcMat[, "zbetamu"]
  
  if (length(sName) > 1) {
    s1 = factor(data[,sName[1]])
    s2 = factor(data[,sName[2]])
    nSubj1 = length(unique(s1)) # should be same as max(s)
    nSubj2 = length(unique(s2)) # should be same as max(s)
    beta0mugrp1 <- mcmcMat[,grep("^zbeta0grp1mu$|^zbeta0grp1mu\\[",colnames(mcmcMat))]
    betamugrp1 <- mcmcMat[,grep("^zbetagrp1mu$|^zbetagrp1mu\\[",colnames(mcmcMat))]
  } else {
    s2 = factor(data[,sName])
    nSubj2 = length(unique(s2)) # should be same as max(s)
  }

  
  
  #-----------------------------------------------------------------------------
  if ( pairsPlot ) {
    
    createPairPlot <- function(mat, labels) {
      # Plot the parameters pairwise, to see correlations:
      openGraph()
      nPtToPlot = 1000
      plotIdx = floor(seq(1,chainLength,by=chainLength/nPtToPlot))
      panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...) {
        usr = par("usr"); on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        r = (cor(x, y))
        txt = format(c(r, 0.123456789), digits=digits)[1]
        txt = paste(prefix, txt, sep="")
        if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
        text(0.5, 0.5, txt, cex=1.5 ) # was cex=cex.cor*r
      }
      pairs(mat, labels, lower.panel=panel.cor , col="skyblue" )
      
      if ( !is.null(saveName) ) {
        saveGraph( file=paste(saveName,"PostPairs",sep=""), type=saveType)
      }  
      
      # Betas per heirarchy
      for (i in 1:Nsubj2) {
        betaNms <- grep(paste0('beta\\[\\d,', i, '\\]'), colnames(beta), value = TRUE)
        createPairPlot(cbind(beta0[, i], beta[, betaNms] )[plotIdx,],
                       labels=c(paste("Group:", levels(s2)[i], 'beta0'), 
                                paste("Group:", levels(s2)[i], betaNms, "\n",xName)))   
      }
      
      # Group Means
      createPairPlot(cbind(beta0mugrp1, betamugrp1)[plotIdx,],
                     labels=c(paste0("beta[", i, "]"), 
                              paste0(betaNms, "\n",xName)))   
      
      # Overall Amts
      createPairPlot(cbind(beta0mu, betamu)[plotIdx,],
                     labels=c(paste0("beta[", i, "]"), 
                              paste0(betaNms, "\n",xName)))   
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
  
  #....................................
  # Original scale:

  # Intercept
  Nidx = ncol(beta0)
  openGraph(width=2.5*Nidx,height=2.0*Nidx)
  layout(mat = matrix(c(rep(1, Nidx), 2:(Nidx^2 + 1)), 
                      Nidx + 1, Nidx, byrow = TRUE),
         heights = c(1, rep(3, Nidx)))
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, 'Original Scale Beta0: Intercept', cex = 2)
  for ( t1Idx in 1:Nidx ) {
    for ( t2Idx in 1:Nidx ) {
      parName1 = paste0('beta0[', t1Idx, "]")
      parName2 = paste0('beta0[', t2Idx, "]")
      grp1 <- levels(s2)[t1Idx]
      grp2 <- levels(s2)[t2Idx]
      if ( t1Idx > t2Idx) {  
        # plot.new() # empty plot, advance to next
        par( mar=c(3.5,3.5,1,1) , mgp=c(2.0,0.7,0) , pty="s" )
        nToPlot = 700
        ptIdx = round(seq(1,chainLength,length=nToPlot))
        plot(mcmcMat[ptIdx,parName2], mcmcMat[ptIdx,parName1], cex.lab=1.75,
               xlab=grp2 , ylab=grp1 , col="skyblue" )
        abline(0,1,lty="dotted")
      } else if (t1Idx == t2Idx) {
        par(mar=c(3.5,1,1,1) , mgp=c(2.0,0.7,0) , pty="m" )
        postInfo = plotPost(mcmcMat[,parName1], cex.lab = 1.75 , 
                             compVal=compVal, ROPE=rope , cex.main=1.5 ,
                             xlab=grp1, main="" )
        includeRows = ( s2 == levels(s2)[t1Idx] ) # rows of this subject in data
        dataPropor = sum(y[includeRows])/sum(includeRows) 
        points( dataPropor , 0 , pch="+" , col="red" , cex=3 )
      } else if ( t1Idx < t2Idx ) {
        par( mar=c(3.5,1,1,1) , mgp=c(2.0,0.7,0) , pty="m" )
        postInfo = plotPost( mcmcMat[,parName1]-mcmcMat[,parName2] , 
                             compVal=compValDiff , ROPE=ropeDiff , cex.main=1.5 ,
                             xlab=paste0(grp1,"-",grp2) , main="" , 
                             cex.lab=1.75 )
        includeRows1 = ( s2 == levels(s2)[t1Idx] ) # rows of this subject in data
        dataPropor1 = sum(y[includeRows1])/sum(includeRows1) 
        includeRows2 = ( s2 == levels(s2)[t2Idx] ) # rows of this subject in data
        dataPropor2 = sum(y[includeRows2])/sum(includeRows2) 
        points( dataPropor1-dataPropor2 , 0 , pch="+" , col="red" , cex=3 )
      }
    }
  }
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName, 'Standardized_Scale_Intercept',
                          sep=""), type=saveType)
  }
  
  # Betas
  for (pred in 1:length(xName)) {
    Nidx = nSubj2
    openGraph(width=2.5*Nidx,height=2.0*Nidx)
    layout(mat = matrix(c(rep(1, Nidx), 2:(Nidx^2 + 1)), 
                        Nidx + 1, Nidx, byrow = TRUE),
           heights = c(1, rep(3, Nidx)))
    par(mar = c(0,0,0,0))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, paste('Original Scale Beta:', xName[pred]), cex = 2)
    for ( t1Idx in 1:Nidx ) {
      for ( t2Idx in 1:Nidx ) {
        parName1 = paste0('beta[', pred, ",", t1Idx, "]")
        parName2 = paste0('beta[', pred, ",", t2Idx, "]")
        grp1 <- levels(s2)[t1Idx]
        grp2 <- levels(s2)[t2Idx]
        if ( t1Idx > t2Idx) {  
          # plot.new() # empty plot, advance to next
          par( mar=c(3.5,3.5,1,1) , mgp=c(2.0,0.7,0) , pty="s" )
          nToPlot = 700
          ptIdx = round(seq(1,chainLength,length=nToPlot))
          plot(mcmcMat[ptIdx,parName2], mcmcMat[ptIdx,parName1], cex.lab=1.75,
               xlab=grp2 , ylab=grp1 , col="skyblue" )
          abline(0,1,lty="dotted")
        } else if (t1Idx == t2Idx) {
          par(mar=c(3.5,1,1,1) , mgp=c(2.0,0.7,0) , pty="m" )
          postInfo = plotPost(mcmcMat[,parName1], cex.lab = 1.75 , 
                              compVal=compVal, ROPE=rope , cex.main=1.5 ,
                              xlab=grp1, main="" )
          includeRows = ( s2 == levels(s2)[t1Idx] ) # rows of this subject in data
          dataPropor = sum(y[includeRows])/sum(includeRows) 
          points( dataPropor , 0 , pch="+" , col="red" , cex=3 )
        } else if ( t1Idx < t2Idx ) {
          par( mar=c(3.5,1,1,1) , mgp=c(2.0,0.7,0) , pty="m" )
          postInfo = plotPost( mcmcMat[,parName1]-mcmcMat[,parName2] , 
                               compVal=compValDiff , ROPE=ropeDiff , cex.main=1.5 ,
                               xlab=paste0(grp1,"-",grp2) , main="" , 
                               cex.lab=1.75 )
          includeRows1 = ( s2 == levels(s2)[t1Idx] ) # rows of this subject in data
          dataPropor1 = sum(y[includeRows1])/sum(includeRows1) 
          includeRows2 = ( s2 == levels(s2)[t2Idx] ) # rows of this subject in data
          dataPropor2 = sum(y[includeRows2])/sum(includeRows2) 
          points( dataPropor1-dataPropor2 , 0 , pch="+" , col="red" , cex=3 )
        }
      }
    }
    if ( !is.null(saveName) ) {
      saveGraph( file=paste(saveName, 'Original_Scale_Beta', 
                            xName[pred],sep=""), type=saveType)
    }
  }
  
  #....................................
  # Standardized scale:
  Nidx = ncol(beta0)
  openGraph(width=2.5*Nidx,height=2.0*Nidx)
  layout(mat = matrix(c(rep(1, Nidx), 2:(Nidx^2 + 1)), 
                      Nidx + 1, Nidx, byrow = TRUE),
         heights = c(1, rep(3, Nidx)))
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, 'Standardized Scale Beta0: Intercept', cex = 2)
  for ( t1Idx in 1:Nidx ) {
    for ( t2Idx in 1:Nidx ) {
      parName1 = paste0('zbeta0[', t1Idx, "]")
      parName2 = paste0('zbeta0[', t2Idx, "]")
      grp1 <- levels(s2)[t1Idx]
      grp2 <- levels(s2)[t2Idx]
      if ( t1Idx > t2Idx) {  
        # plot.new() # empty plot, advance to next
        par( mar=c(3.5,3.5,1,1) , mgp=c(2.0,0.7,0) , pty="s" )
        nToPlot = 700
        ptIdx = round(seq(1,chainLength,length=nToPlot))
        plot(mcmcMat[ptIdx,parName2], mcmcMat[ptIdx,parName1], cex.lab=1.75,
             xlab=grp2 , ylab=grp1 , col="skyblue" )
        abline(0,1,lty="dotted")
      } else if (t1Idx == t2Idx) {
        par(mar=c(3.5,1,1,1) , mgp=c(2.0,0.7,0) , pty="m" )
        postInfo = plotPost(mcmcMat[,parName1], cex.lab = 1.75 , 
                            compVal=compVal, ROPE=rope , cex.main=1.5 ,
                            xlab=grp1, main="" )
        includeRows = ( s2 == levels(s2)[t1Idx] ) # rows of this subject in data
        dataPropor = sum(y[includeRows])/sum(includeRows) 
        points( dataPropor , 0 , pch="+" , col="red" , cex=3 )
      } else if ( t1Idx < t2Idx ) {
        par( mar=c(3.5,1,1,1) , mgp=c(2.0,0.7,0) , pty="m" )
        postInfo = plotPost( mcmcMat[,parName1]-mcmcMat[,parName2] , 
                             compVal=compValDiff , ROPE=ropeDiff , cex.main=1.5 ,
                             xlab=paste0(grp1,"-",grp2) , main="" , 
                             cex.lab=1.75 )
        includeRows1 = ( s2 == levels(s2)[t1Idx] ) # rows of this subject in data
        dataPropor1 = sum(y[includeRows1])/sum(includeRows1) 
        includeRows2 = ( s2 == levels(s2)[t2Idx] ) # rows of this subject in data
        dataPropor2 = sum(y[includeRows2])/sum(includeRows2) 
        points( dataPropor1-dataPropor2 , 0 , pch="+" , col="red" , cex=3 )
      }
    }
  }
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName, 'Standardized_Scale_Intercept', sep=""), 
               type=saveType)
  }
  
  
  # Betas
  for (pred in 1:length(xName)) {
    Nidx = nSubj2
    openGraph(width=2.5*Nidx,height=2.0*Nidx)
    layout(mat = matrix(c(rep(1, Nidx), 2:(Nidx^2 + 1)), 
                        Nidx + 1, Nidx, byrow = TRUE),
           heights = c(1, rep(3, Nidx)))
    par(mar = c(0,0,0,0))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, paste('Standardized Scale Beta:', xName[pred]), cex = 2)
    for ( t1Idx in 1:Nidx ) {
      for ( t2Idx in 1:Nidx ) {
        parName1 = paste0('zbeta[', pred, ",", t1Idx, "]")
        parName2 = paste0('zbeta[', pred, ",", t2Idx, "]")
        grp1 <- levels(s2)[t1Idx]
        grp2 <- levels(s2)[t2Idx]
        if ( t1Idx > t2Idx) {  
          # plot.new() # empty plot, advance to next
          par( mar=c(3.5,3.5,1,1) , mgp=c(2.0,0.7,0) , pty="s" )
          nToPlot = 700
          ptIdx = round(seq(1,chainLength,length=nToPlot))
          plot(mcmcMat[ptIdx,parName2], mcmcMat[ptIdx,parName1], cex.lab=1.75,
               xlab=grp2 , ylab=grp1 , col="skyblue" )
          abline(0,1,lty="dotted")
        } else if (t1Idx == t2Idx) {
          par(mar=c(3.5,1,1,1) , mgp=c(2.0,0.7,0) , pty="m" )
          postInfo = plotPost(mcmcMat[,parName1], cex.lab = 1.75 , 
                              compVal=compVal, ROPE=rope , cex.main=1.5 ,
                              xlab=grp1, main="" )
          includeRows = ( s2 == levels(s2)[t1Idx] ) # rows of this subject in data
          dataPropor = sum(y[includeRows])/sum(includeRows) 
          points( dataPropor , 0 , pch="+" , col="red" , cex=3 )
        } else if ( t1Idx < t2Idx ) {
          par( mar=c(3.5,1,1,1) , mgp=c(2.0,0.7,0) , pty="m" )
          postInfo = plotPost( mcmcMat[,parName1]-mcmcMat[,parName2] , 
                               compVal=compValDiff , ROPE=ropeDiff , cex.main=1.5 ,
                               xlab=paste0(grp1,"-",grp2) , main="" , 
                               cex.lab=1.75 )
          includeRows1 = ( s2 == levels(s2)[t1Idx] ) # rows of this subject in data
          dataPropor1 = sum(y[includeRows1])/sum(includeRows1) 
          includeRows2 = ( s2 == levels(s2)[t2Idx] ) # rows of this subject in data
          dataPropor2 = sum(y[includeRows2])/sum(includeRows2) 
          points( dataPropor1-dataPropor2 , 0 , pch="+" , col="red" , cex=3 )
        }
      }
    }
    if ( !is.null(saveName) ) {
      saveGraph( file=paste(saveName, 'Standardized_Scale_Beta', 
                            xName[pred],sep=""), type=saveType)
    }
  }
  
  if (length(sName) > 1) {
    #.................................................
    # Group comparisons
    Nidx = ncol(beta0mugrp1)
    openGraph(width=2.5*Nidx,height=2.0*Nidx)
    layout(mat = matrix(c(rep(1, Nidx), 2:(Nidx^2 + 1)), 
                        Nidx + 1, Nidx, byrow = TRUE),
           heights = c(1, rep(3, Nidx)))
    par(mar = c(0,0,0,0))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, 'Standardized Scale Group Intercept Mu', cex = 2)
    for ( t1Idx in 1:Nidx ) {
      for ( t2Idx in 1:Nidx ) {
        parName1 = paste0('zbeta0grp1mu[', t1Idx, "]")
        parName2 = paste0('zbeta0grp1mu[', t2Idx, "]")
        grp1 <- levels(s1)[t1Idx]
        grp2 <- levels(s1)[t2Idx]
        if ( t1Idx > t2Idx) {  
          # plot.new() # empty plot, advance to next
          par( mar=c(3.5,3.5,1,1) , mgp=c(2.0,0.7,0) , pty="s" )
          nToPlot = 700
          ptIdx = round(seq(1,chainLength,length=nToPlot))
          plot(mcmcMat[ptIdx,parName2], mcmcMat[ptIdx,parName1], cex.lab=1.75,
               xlab=grp2 , ylab=grp1 , col="skyblue" )
          abline(0,1,lty="dotted")
        } else if (t1Idx == t2Idx) {
          par(mar=c(3.5,1,1,1) , mgp=c(2.0,0.7,0) , pty="m" )
          postInfo = plotPost(mcmcMat[,parName1], cex.lab = 1.75 , 
                              compVal=compVal, ROPE=rope , cex.main=1.5 ,
                              xlab=grp1, main="" )
          includeRows = ( s1 == levels(s1)[t1Idx] ) # rows of this subject in data
          dataPropor = sum(y[includeRows])/sum(includeRows) 
          points( dataPropor , 0 , pch="+" , col="red" , cex=3 )
        } else if ( t1Idx < t2Idx ) {
          par( mar=c(3.5,1,1,1) , mgp=c(2.0,0.7,0) , pty="m" )
          postInfo = plotPost( mcmcMat[,parName1]-mcmcMat[,parName2] , 
                               compVal=compValDiff , ROPE=ropeDiff , cex.main=1.5 ,
                               xlab=paste0(grp1,"-",grp2) , main="" , 
                               cex.lab=1.75 )
          includeRows1 = ( s1 == levels(s1)[t1Idx] ) # rows of this subject in data
          dataPropor1 = sum(y[includeRows1])/sum(includeRows1) 
          includeRows2 = ( s1 == levels(s1)[t2Idx] ) # rows of this subject in data
          dataPropor2 = sum(y[includeRows2])/sum(includeRows2) 
          points( dataPropor1-dataPropor2 , 0 , pch="+" , col="red" , cex=3 )
        }
      }
    }
    if ( !is.null(saveName) ) {
      saveGraph( file=paste(saveName,"Group_Mu_Intercept",sep=""), type=saveType)
    }
    
    
    # Betas
    Nidx = nSubj1
    openGraph(width=2.5*Nidx,height=2.0*Nidx)
    layout(mat = matrix(c(rep(1, Nidx), 2:(Nidx^2 + 1)), 
                        Nidx + 1, Nidx, byrow = TRUE),
           heights = c(1, rep(3, Nidx)))
    par(mar = c(0,0,0,0))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, 'Standardized Scale Group Beta Mu', cex = 2)
    for ( t1Idx in 1:Nidx ) {
      for ( t2Idx in 1:Nidx ) {
        parName1 = paste0('zbetagrp1mu[', t1Idx, "]")
        parName2 = paste0('zbetagrp1mu[', t2Idx, "]")
        grp1 <- levels(s1)[t1Idx]
        grp2 <- levels(s1)[t2Idx]
        if ( t1Idx > t2Idx) {  
          # plot.new() # empty plot, advance to next
          par( mar=c(3.5,3.5,1,1) , mgp=c(2.0,0.7,0) , pty="s" )
          nToPlot = 700
          ptIdx = round(seq(1,chainLength,length=nToPlot))
          plot(mcmcMat[ptIdx,parName2], mcmcMat[ptIdx,parName1], cex.lab=1.75,
               xlab=grp2 , ylab=grp1 , col="skyblue" )
          abline(0,1,lty="dotted")
        } else if (t1Idx == t2Idx) {
          par(mar=c(3.5,1,1,1) , mgp=c(2.0,0.7,0) , pty="m" )
          postInfo = plotPost(mcmcMat[,parName1], cex.lab = 1.75 , 
                              compVal=compVal, ROPE=rope , cex.main=1.5 ,
                              xlab=grp1, main="" )
          includeRows = ( s1 == levels(s1)[t1Idx] ) # rows of this subject in data
          dataPropor = sum(y[includeRows])/sum(includeRows) 
          points( dataPropor , 0 , pch="+" , col="red" , cex=3 )
        } else if ( t1Idx < t2Idx ) {
          par( mar=c(3.5,1,1,1) , mgp=c(2.0,0.7,0) , pty="m" )
          postInfo = plotPost( mcmcMat[,parName1]-mcmcMat[,parName2] , 
                               compVal=compValDiff , ROPE=ropeDiff , cex.main=1.5 ,
                               xlab=paste0(grp1,"-",grp2) , main="" , 
                               cex.lab=1.75 )
          includeRows1 = ( s1 == levels(s1)[t1Idx] ) # rows of this subject in data
          dataPropor1 = sum(y[includeRows1])/sum(includeRows1) 
          includeRows2 = ( s1 == levels(s1)[t2Idx] ) # rows of this subject in data
          dataPropor2 = sum(y[includeRows2])/sum(includeRows2) 
          points( dataPropor1-dataPropor2 , 0 , pch="+" , col="red" , cex=3 )
        }
      }
    }
    if ( !is.null(saveName) ) {
      saveGraph( file=paste(saveName,"Group_Mu_Beta",sep=""), type=saveType)
    }
    #   openGraph(width=2.5*Nidx,height=2.0*Nidx)
    # layout(mat = matrix(c(rep(1, Nidx), 2:(Nidx^2 + 1)), 
    #                     Nidx + 1, Nidx, byrow = TRUE),
    #        heights = c(1, rep(3, Nidx)))
    # par(mar = c(0,0,0,0))
    # plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    # text(x = 0.5, y = 0.5, paste('Original Scale Group Mu:', xName[pred]), cex = 2)
    #   for ( t1Idx in 1:Nidx ) {
    #     for ( t2Idx in 1:Nidx ) {
    #       parName1 = paste0('zbetagrp1mu[', pred, ",", t1Idx, "]")
    #       parName2 = paste0('zbetagrp1mu[', pred, ",", t2Idx, "]")
    #       grp1 <- levels(s1)[t1Idx]
    #       grp2 <- levels(s1)[t2Idx]
    #       if ( t1Idx > t2Idx) {  
    #         # plot.new() # empty plot, advance to next
    #         par( mar=c(3.5,3.5,1,1) , mgp=c(2.0,0.7,0) , pty="s" )
    #         nToPlot = 700
    #         ptIdx = round(seq(1,chainLength,length=nToPlot))
    #         plot(mcmcMat[ptIdx,parName2], mcmcMat[ptIdx,parName1], cex.lab=1.75,
    #              xlab=grp2 , ylab=grp1 , col="skyblue" )
    #         abline(0,1,lty="dotted")
    #       } else if (t1Idx == t2Idx) {
    #         par(mar=c(3.5,1,1,1) , mgp=c(2.0,0.7,0) , pty="m" )
    #         postInfo = plotPost(mcmcMat[,parName1], cex.lab = 1.75 , 
    #                             compVal=compVal, ROPE=rope , cex.main=1.5 ,
    #                             xlab=grp1, main="" )
    #         includeRows = ( s1 == levels(s1)[t1Idx] ) # rows of this subject in data
    #         dataPropor = sum(y[includeRows])/sum(includeRows) 
    #         points( dataPropor , 0 , pch="+" , col="red" , cex=3 )
    #       } else if ( t1Idx < t2Idx ) {
    #         par( mar=c(3.5,1,1,1) , mgp=c(2.0,0.7,0) , pty="m" )
    #         postInfo = plotPost( mcmcMat[,parName1]-mcmcMat[,parName2] , 
    #                              compVal=compValDiff , ROPE=ropeDiff , cex.main=1.5 ,
    #                              xlab=paste0(grp1,"-",grp2) , main="" , 
    #                              cex.lab=1.75 )
    #         includeRows1 = ( s1 == levels(s1)[t1Idx] ) # rows of this subject in data
    #         dataPropor1 = sum(y[includeRows1])/sum(includeRows1) 
    #         includeRows2 = ( s1 == levels(s1)[t2Idx] ) # rows of this subject in data
    #         dataPropor2 = sum(y[includeRows2])/sum(includeRows2) 
    #         points( dataPropor1-dataPropor2 , 0 , pch="+" , col="red" , cex=3 )
    #       }
    #     }
    #   }
    # }
  }
  
  #.................................................
  # Pop comparisons
  panelCount = 1  
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"Pop_Info") )
    
  histInfo = plotPost( zbeta0mu , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(zmu[beta*0]) , main=paste("Intercept, Pop Level") )
  
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"Pop Info") )
  histInfo = plotPost( zbetamu , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(zmu[beta*1]) , main=paste("Slope, Pop Level") )

  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"Pop_Info",sep=""), type=saveType)
  }
  #-----------------------------------------------------------------------------
}
#===============================================================================
