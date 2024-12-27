###############################################################################
# Load all necessary functions to run a simulation study with summary statistics
###############################################################################
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

online_source <- "https://raw.githubusercontent.com/Adrifelcha/EZ-project/refs/heads/main/code/functions/" 
######################
# Simulation set up
######################
# A function to name output file
source(paste(online_source, "namingFunctions.R", sep=""))
# Load default priors for the hierarchical mean parameters
source("./default_priors.R")

################
# JAGS set up 
################
# Write JAGS model
source(paste(online_source, "write_JAGSmodel.R", sep=""))
# Function to pass data to JAGS
source(paste(online_source, "data_toJAGS.R", sep=""))
# Function to load default init values
source(paste(online_source, "default_inits.R",sep=""))
source("./runJAGS.R")
##################
# MAnage samples
##################
# Function to extract samples
source(paste(online_source, "extractSamples.R", sep=""))
# Function to calculate Rhat
source("https://raw.githubusercontent.com/Adrifelcha/EZ-project/refs/heads/main/code/functions/rhat.R")
# Function to double check all rhats < 1.05
source(paste(online_source, "check_Rhat.R", sep=""))

# Function to run simulations
source("./main_runSims.R")

#####################
# Show results
#####################
# Show results of a single simulation study
source("https://raw.githubusercontent.com/Adrifelcha/EZ-project/refs/heads/main/code/functions/plot_recoveryPanel.R")
source("https://raw.githubusercontent.com/Adrifelcha/EZ-project/refs/heads/main/code/functions/plot_singleSim.R")