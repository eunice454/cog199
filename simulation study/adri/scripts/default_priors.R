default_priors <- function(modelType=NA){
  prior <- data.frame("bound_mean_mean" = 1.50,  "bound_mean_sdev" = 0.20, 
                      "drift_mean_mean" = 0.00,  "drift_mean_sdev" = 0.30, 
                      "nondt_mean_mean" = 0.40,  "nondt_mean_sdev" = 0.05)
  if(modelType!="hierarchical"){
    prior$betaweight_mean = 0
    prior$betaweight_sdev = 1
  }
  return(prior)
}