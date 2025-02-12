runJAGS <- function(summaryData, nTrials, X, jagsData, jagsParameters, jagsInits, 
                    n.chains=4, n.burnin=250, n.iter=2000, n.thin=1, modelFile="./EZHBDDM.bug",
                    track_allParameters = track_allParameters){
  # Prepare data
  correct <- summaryData$A
  varRT   <- summaryData$Vrt
  meanRT  <- summaryData$Mrt
  n_trials <- nTrials
  n_Participants    <- nrow(summaryData)
  # Run model and get samples
  tic <- clock::date_now(zone="UTC")
  suppressMessages(samples <- jags(data=jagsData, 
                                   parameters.to.save=jagsParameters, 
                                   model=modelFile, 
                                   n.chains=n.chains, 
                                   n.iter=n.iter, 
                                   n.burnin=n.burnin, 
                                   n.thin=n.thin, 
                                   DIC=T, 
                                   inits=jagsInits))
  toc <- clock::date_now(zone="UTC")
  clock <- as.numeric(toc-tic, units="secs")  # Record time
  object <- samples$BUGSoutput$sims.array
  rhats <- apply(object,3,getRhat)
  estimates <- list()
  error <- list
  credInterval <- list()
  for(i in 1:length(jagsParameters)){
    posteriorParameters <- extractSamples(jagsParameters[i], samples)
    if(length(dim(posteriorParameters))==3){
      meanPost    <- apply(posteriorParameters,3,mean)
      sdPost      <- apply(posteriorParameters,3,sd)
      percentiles <- apply(posteriorParameters,3, quantile, probs=c(0.025,0.975))
    }else{
      meanPost    <- mean(posteriorParameters)
      sdPost      <- sd(posteriorParameters)
      percentiles <- quantile(posteriorParameters,probs = c(0.025,0.975))
    }
    estimates    <- c(estimates, list(meanPost))
    error        <- c(error, list(sdPost))
    credInterval <- c(credInterval, list(percentiles))
  }
  names(estimates)   <- jagsParameters
  names(credInterval) <- jagsParameters
  
  return(list("estimates" = estimates, "sd" = error, "credInterval" = credInterval,
              "rhats" = rhats, "clock" = clock))
}