runSims <- function(nParticipants, nTrials, nDatasets = 10, modelType = NA, criterion = NA, n.chains = 3, 
                    n.burnin=250, n.iter=2000, n.thin=1, forceSim = FALSE, rhatCheck=TRUE, modelFile = NA,
                    track_allParameters = FALSE, redo_if_bad_rhat=FALSE, output.folder="../results/"){
  grand_tic <- clock::date_now(zone="UTC")
  #################################
  # Initial checks
  #################################
  # Load necessary R libraries
  suppressMessages(library(R2jags))
  # Make sure modelType is valid
  if(is.na(modelType)){    modelType = "hierarchical"
  }else{  valid.models <- c("hierarchical", "metaregression", "ttest") 
        if(!(modelType %in% valid.models)){
          stop("Please specify a valid modelType: 'hierarchical' (default), 'metaregression' 'ttest'")
        }
  }
  # Identify output File
  outputFile <- nameOutput(nTrials, nParticipants, nDatasets, modelType, fromPrior=TRUE, output.folder = output.folder)
 
  # Check if we needToRun simulations again (overruled by 'forceSim')
  if(!forceSim){
    if(file.exists(outputFile)){    
              needToRun <- FALSE
              load(outputFile)
    }else{ needToRun <- TRUE}
  }else{  needToRun <- TRUE  }
  ###################################
  # Run simulation study (if needed)
  ###################################
  if(needToRun){
    ######################
    #      SET UP        #
    ######################
    settings <- list("nPart"= nParticipants, "nTrials"= nTrials,
                     "modelType" = modelType, "nDatasets" = nDatasets,
                     "prior" = default_priors(modelType))
    # If the model includes an effect (betaweight)
    if(modelType!="hierarchical"){   
      # Make sure we have a valid "criterion" (default to 'drift')
      if(is.na(criterion)){    criterion <- "drift"   }
      X <- 0:nParticipants   # Default predictor       
      if(modelType=="ttest"){   X <- X %% 2    # Dummy predictor
      }else{   X <- X/nParticipants          }        
      settings <- c(settings, list("X" = X, "criterion" = criterion))
    }else{    X <- NA    }

    # ~~~~~~~~~~~~~~~~ JAGS variables
    # Define parameters to be tracked on JAGS, according to the modelType
    if(track_allParameters){jagsParameters <- c("bound_mean", "drift_mean", "nondt_mean", "bound", "nondt",
                                                "drift_sdev", "nondt_sdev", "bound_sdev", "drift")
    }else{                  jagsParameters <- c("bound_mean", "drift_mean", "nondt_mean")                    }
    if(modelType!="hierarchical"){  jagsParameters <- c(jagsParameters, "betaweight")  }
    # Data to be passed to JAGS
    jagsData <- list("n_Participants", "n_trials", "correct", "meanRT", "varRT")
    # init values
    jagsInits <- default_inits(n.chains, nParticipants)  
    # ~~~~~~~~~~~~~~~~ Storing objects
    # Count number of parameters
    if(track_allParameters){  nParams <- (length(jagsParameters)-3) + (nParticipants*3)    
    }else{                    nParams <- (length(jagsParameters))                           }
    MatEstimates <- matrix(NA, nrow=nDatasets, ncol=nParams)
    MatTrueVal   <- matrix(NA, nrow=nDatasets, ncol=nParams)
    ArrayCredInt <- array(NA, dim=c(nDatasets,nParams,2))
    MatRhats     <- matrix(NA, nrow=nDatasets, ncol=(nParams+1))
    seed_id <- rep(NA,nDatasets)
    # Provide hierarchical mean parameters
    bound_mean <- rnorm(nDatasets, settings$prior$bound_mean_mean, 
                        settings$prior$bound_mean_sdev)
    drift_mean <- rnorm(nDatasets, settings$prior$drift_mean_mean, 
                        settings$prior$drift_mean_sdev)
    nondt_mean <- rnorm(nDatasets, settings$prior$nondt_mean_mean, 
                        settings$prior$nondt_mean_sdev)
    if(modelType!="hierarchical"){
        if(criterion=="nondt"){   betaweight <- runif(nDatasets,  0, 1)
        }else{                    betaweight <- runif(nDatasets, -1, 1)       }
    }
    ######################
    #   Run iterations   #
    ######################
    repetition_counts <- 0
    for(k in 1:nDatasets){
      rhat_not_verified <- TRUE
      seed <- k
      set.seed(seed)
      cat("============>> Dataset", k, "of", nDatasets,"\n")
      # Get data
      #########
      indiv_parameters <- getIndPar(nParticipants, drift_mean[k], bound_mean[k], nondt_mean[k])
      summData = simSumStats(indiv_parameters, nTrials)
      parameter_set <- list("bound" = indiv_parameters$bound,
                            "drift" = indiv_parameters$drift,
                            "nondt" = indiv_parameters$nondt,
                            "bound_mean" = bound_mean[k],
                            "drift_mean" = drift_mean[k],
                            "nondt_mean" = nondt_mean[k])
      # Obtain estimates  
      ###################
      while(rhat_not_verified){
        if(k>100){
          runJags <- try(runJAGS(summaryData = summData, nTrials = nTrials, X = X, 
                                 jagsData = jagsData, jagsParameters = jagsParameters, 
                                 jagsInits = jagsInits, n.chains = n.chains, 
                                 modelFile = modelFile, track_allParameters = track_allParameters))
          while(inherits(runJags, "try-error")){
            repetition_counts <- repetition_counts+1
            seed <- seed+10000
            set.seed(seed)
            cat("============>> Dataset", k, "of", nDatasets,"+",repetition_counts,"\n")
            bound.mean <- rnorm(1, settings$prior$bound_mean_mean, 
                                settings$prior$bound_mean_sdev)
            drift.mean <- rnorm(1, settings$prior$drift_mean_mean, 
                                settings$prior$drift_mean_sdev)
            nondt.mean <- rnorm(1, settings$prior$nondt_mean_mean, 
                                settings$prior$nondt_mean_sdev)
            if(modelType!="hierarchical"){
              if(criterion=="nondt"){   beta <- runif(1,  0, 1)
              }else{                    beta <- runif(1, -1, 1)       }
            }
            indiv_parameters <- getIndPar(n_participants = nParticipants, dmean = drift.mean, 
                                          bmean = bound.mean, nmean = nondt.mean)
            summData = simSumStats(indiv_pars = indiv_parameters, n_trials =  nTrials)
            parameter_set <- list("bound" = indiv_parameters$bound,
                                  "drift" = indiv_parameters$drift,
                                  "nondt" = indiv_parameters$nondt,
                                  "bound_mean" = bound.mean,
                                  "drift_mean" = drift.mean,
                                  "nondt_mean" = nondt.mean)
            runJags <- try(runJAGS(summaryData = summData, nTrials = nTrials, X = X, 
                                   jagsData = jagsData, jagsParameters = jagsParameters, 
                                   jagsInits = jagsInits, n.chains = n.chains, 
                                   modelFile = modelFile, track_allParameters = track_allParameters))
          }
        }else{
          runJags <- runJAGS(summaryData = summData, nTrials = nTrials, X = X, 
                             jagsData = jagsData, jagsParameters = jagsParameters, 
                             jagsInits = jagsInits, n.chains = n.chains, 
                             modelFile = modelFile, track_allParameters = track_allParameters)
        }
        
        count_bad_rhats <- sum(runJags$rhats[jagsParameters]>1.05)
        if((!redo_if_bad_rhat)|(count_bad_rhats==0)){ rhat_not_verified <-  FALSE}
        seed <- seed+10000
      }
      MatRhats[k,] <- runJags$rhats
      c <- 0; d <- 0
      for(j in 1:length(runJags$estimates)){
        this <-  names(runJags$estimates[j])
        m <- length(runJags$estimates[[this]])
        w <- length(parameter_set[[this]])
        MatEstimates[k,(c+1):(c+m)] <- runJags$estimates[[this]]
        MatTrueVal[k,(d+1):(d+w)]   <- parameter_set[[this]]
        if(is.vector(runJags$credInterval[[j]])){
          ArrayCredInt[k,(c+1):(c+m),1] <- runJags$credInterval[[this]][1]
          ArrayCredInt[k,(c+1):(c+m),2] <- runJags$credInterval[[this]][2]
        }else{
          ArrayCredInt[k,(c+1):(c+m),1] <- runJags$credInterval[[this]][1,]
          ArrayCredInt[k,(c+1):(c+m),2] <- runJags$credInterval[[this]][2,]
        }
        c <- c+m; d <- d+w
      }
      
      seed_id[k] <- seed
    }
    
    paramNames <- NA
    paramNames2 <- NA
    for(j in 1:length(runJags$estimates)){
      this <-  names(runJags$estimates[j])
      if(is.vector(runJags$credInterval[[this]])){
        paramNames <- c(paramNames, names(runJags$credInterval[this]))
      }else{
        paramNames <- c(paramNames, colnames(runJags$credInterval[[this]]))
      }
      if(length(parameter_set[[this]])==1){
        paramNames2 <- c(paramNames2, names(parameter_set[this]))
      }else{
        labels <- paste(names(parameter_set[this]), "[",1:length(parameter_set[[this]]),"]",sep="")
        paramNames2 <- c(paramNames2, labels)
      }
    }
    paramNames <- paramNames[-1]
    paramNames2 <- paramNames2[-1]
    colnames(MatEstimates) <- paramNames
    colnames(ArrayCredInt) <- paramNames
    colnames(MatTrueVal)   <- paramNames2
    colnames(MatRhats) <- names(runJags$rhats)
    
    if(rhatCheck){check_Rhat(MatRhats)}
    
    grand_toc <- clock::date_now(zone="UTC")
    total_time <- difftime(grand_toc, grand_tic, units="mins")
    output <- list("rhats"  = MatRhats, "estimates" = MatEstimates, "credIntervals" = ArrayCredInt,
                   "trueValues" = MatTrueVal, "settings" = settings, "n.chains" = n.chains,
                   "totalTime" = total_time, "seed_id" = seed_id)
    save(output, file=outputFile)
    
    cat("Running this simulation study took ", total_time, "minutes.\n")
  }else{  cat("This simulation had been run before.\nLoading stored results: COMPLETE!\n",
              "Running this simulation study took ", output$totalTime, "minutes.\n")
    if(rhatCheck){   check_Rhat(output$rhats)      }
  }
  return(output)
}