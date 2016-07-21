##***********************************************************************
## Functions in this script were tested under the following environment
##
## > sessionInfo()
## R version 3.3.0 (2016-05-03)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: OS X 10.11.5 (El Capitan)
##
## locale:
##  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
##
## attached base packages:
##   [1] stats     graphics  grDevices utils     datasets  methods   base
##
## other attached packages:
##   [1] brms_0.9.1       rstan_2.9.0-3    ggplot2_2.1.0    data.table_1.9.6
##
##************************************************************************
#
#### Load libraries ####
library(data.table)
library("brms")
library("plyr")

#### Extract summaries from BRMS modelresults ####
summarize_multinomial_brms <- function(modelresults) {
  #'
  #' @param modelresults from the BRMS model
  #' @example
  #' modelresults <- brms_modelresults
  #'
  #' TODO: add gaussian summaries and generalize for logistic regression
  if (modelresults$family$family == "categorical" & modelresults$family$link == "logit") {
    #' get all the betas from the model. Burn-in (warm-up) has already been accounted for.
    model_betas <- posterior_samples(x = modelresults, add_chain = T)
    #' get betas for each contrast
    contrast_10 <- as.matrix(model_betas[, grep("b_trait1:", colnames(model_betas), value = T)])
    contrast_20 <- as.matrix(model_betas[, grep("b_trait2:", colnames(model_betas), value = T)])
    colnames(contrast_10) <- gsub("b_trait1:(.*)", "\\1", colnames(contrast_10))
    colnames(contrast_20) <- gsub("b_trait2:(.*)", "\\1", colnames(contrast_20))
    contrast_21 <- contrast_20[, colnames(contrast_10)] - contrast_10
    betalist <- list(contrast_10 = contrast_10,
                     contrast_20 = contrast_20,
                     contrast_21 = contrast_21)
    allstats <- sapply(names(betalist), function(bname) {
      #' @example
      #' bname <- names(betalist)[1]
      #'
      betavals <- betalist[[bname]]
      betaci <- sapply(colnames(betavals), function(bcol) {
        #' @example
        #' bcol <- 4
        #' include probability mass at 0.
        #' include data after discarding burn-in.
        #'
        cival <- c(quantile(x = betavals[, bcol],
                            probs = c(0.025, 0.975)),
                   medbeta = median(betavals[, bcol]))
        names(cival) <- c("lowerci", "upperci", "medbeta")
        significant <- ifelse(cival["lowerci"] * cival["upperci"] > 0, "yes", "no")
        names(significant) <- NULL
        out <- c(cival, significant = significant)
        return(out)
      })
      outstats <- data.frame(t(betaci))
      outstats$variable <- rownames(outstats)
      return(data.table(outstats))
    }, simplify = F)
    return(allstats)
  }
}
# TODO: Gaussian, Binary

#### Compute error rate on model predictions ####
compute_error_rate <- function(modelresults){
  #'
  #' @param modelresults from the BRMS model
  #' @example
  #'
  #' modelresults <- brms_modelresults
  #'
  modelpredictions <- predict(object = modelresults)
  predictions <- ldply(sapply(1:nrow(modelpredictions), function(brow){
    #' brow <- 1
    #'
    out <- modelpredictions[brow, ]/sum(modelpredictions[1,])
    response <- modelresults$data$response[brow]
    # is the class with the highest posterior probability same as the true class?
    classification_error <- which.max(out) != response + 1
    likelihood_error <- 1 - out[response + 1]
    # collect for perplexity | -log(predicted_probability of the observed response)
    log_prob <- -log(out[response + 1])
    errors <- data.frame(classification_error = classification_error,
                         likelihood_error = likelihood_error,
                         logprob = log_prob, stringsAsFactors = F)
    return(errors)
  }, simplify = F))
  err <- apply(predictions, 2, mean)
  return(err)
}

#### Compute R^2 of the regression model ####
get_rsquared <- function(modelresults){
  #'
  #' @param modelresults from the BRMS model
  #' modelresults <- brms_modelresults
  #'
  model_predictions <- predict(object = modelresults)
  # Squared correlation between response and the predicted value
  rsquared <- cor(modelresults$data$response_recoded, model_predictions[ ,"Estimate"]) ^ 2
  # L2 distance between prediction and response
  ss_residual <- sum((modelresults$data$response_recoded - model_predictions[ ,"Estimate"]) ^ 2)
  # variance in the response times n
  ss_total  <- sum((modelresults$data$response_recoded - mean(modelresults$data$response_recoded)) ^ 2)
  # totol response variance explained by covariates
  rsquared_ss <- 1 - (ss_residual/ss_total)
  return(c(cor_obspred_sq = rsquared, var_explained = rsquared_ss))
}
