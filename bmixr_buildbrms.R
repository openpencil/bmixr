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

#### Regression modeling code ####
run_brms <- function(datasub,
                     samples_of_interest = NULL,
                     response_of_interest,
                     covariates,
                     family_of_regression,
                     randomeffects_term,
                     numcores,
                     numiter,
                     name_of_output_file,
                     include_randomslopes = F,
                     interaction_variable1_pattern = NULL,
                     interaction_variable2 = NULL,
                     type_of_prior = NULL,
                     treedepth = 10,
                     adaptdelta = 0.99) {
  #' Estimate regression model with random effects
  #'
  #' @param datasub dataset with response variable and covariates of interest
  #' @param samples_of_interest subset of response values of interest
  #' @param response_of_interest column which contains the response of interest
  #' @param covariates to include in the model (include ALL covariates!!)
  #' @param family_of_regression "multinomial", "binary" or "gaussian"
  #' @param include_randomslopes boolean, if T, individual slopes estimated using random effects
  #' @param interaction_variable1_pattern for covariate with interactions of interest
  #' @param interaction_variable2  covariates for building interaction terms with variable1
  #' @param randomeffects_term variable with the random effects in the model
  #' @param numcores number of cores to use
  #' @param numiter number of MCMC iterations
  #' # TODO: what is a good prior for gaussian models?
  #' @param type_of_prior array e.g. c(set_prior("normal (0, 8)")) |default NULL
  #' @param treedepth depth of tree evaluated at each MCMC iteration (STAN parameter), default 10
  #' @param adaptdelta larger values slow down MCMC, increases accuracy of posteriors, default 0.8
  #' @param name_of_output_file name of file in which the model findings will be saved
  #'
  #' @examples
  #' datasub <- some_dataframe
  #' samples_of_interest <- c(0, 1, 2)
  #' response_of_interest <- "some_response"
  #' covariates <- c("vector", "of", "covariates")
  #' family_of_regression  <- "multinomial"
  #' include_randomslopes <- F
  #' interaction_variable1_pattern <- "pattern_for_the_variable_with_interactions"
  #' interaction_variable2 <- "other_variables_for_pairwise_interaction_terms"
  #' randomeffects_variable <- "randomeffects_variable"
  #' numcores <- 4
  #' numiter <- 10
  #' type_of_prior <- c(set_prior("normal (0, 8)"))
  #' treedepth <- 10
  #' adaptdelta = 0.99
  #'
  #' @output brms regression output will be saved in the working directory
  #'
  if (family_of_regression != "gaussian") {
    regdata <- datasub[get(response_of_interest) %in% samples_of_interest]
    indices_of_interest <- seq_along(samples_of_interest) - 1
    names(indices_of_interest) <- as.character(samples_of_interest)
    regdata[, response_recoded := indices_of_interest[as.character(get(response_of_interest))]]
  } else {
    regdata <- datasub
    regdata[, response_recoded := get(response_of_interest)]
  }
  # remove mathematical operators from colnames so they are not parsed with the regression formula.
  setnames(x = regdata, old = colnames(regdata),
           new = gsub("-", "_", gsub("\\((.*)\\)", "_\\1", colnames(regdata))))
  # formula specifics
  # No global intercept will be estimated: 0
  # trait has c - 1 levels, where c is the number of categories of a response variable,
  # estimates separate intercepts  for each response category
  # trait:covariate estimates estimates the regression coefficients for each level of response.
  # Interactions Y/N
  covs_of_interest <- setdiff(covariates, c(response_of_interest, randomeffects_term))
  if (length(interaction_variable1_pattern) == 0 & length(interaction_variable2) == 0) {
    covariate_vector <- paste(covs_of_interest, collapse = "+")
  } else {
    interaction_vars <- grep(interaction_variable1_pattern, covs_of_interest, value = T)
    non_interaction_vars <- setdiff(covs_of_interest, interaction_vars)
    covariate_vector <- paste(c(non_interaction_vars, interaction_vars,
                                paste(interaction_vars, interaction_variable2, sep = ":")),
                              collapse = "+")
  }
  if (family_of_regression == "multinomial") {
    formula_1 <- sprintf("response_recoded ~ 0 + trait + trait:(%s)", covariate_vector)
  } else {
    formula_1 <- sprintf("response_recoded ~ %s", covariate_vector)
  }
  # Random Slopes Y/N
  if (include_randomslopes == F) {
    cat("no_randomslopes", randomeffects_term, "\n")
    formula <- paste(formula_1, sprintf("(1|%s)", randomeffects_term), sep = " + ")
  } else {
    # NB: Randomslopes will not work for more than 33 variables. Don't ask why.
    cat("including_randomslopes", randomeffects_term, "\n")
    formula <- paste(formula_1, sprintf("(1+%s|%s)", covariate_vector, randomeffects_term), sep = " + ")
  }
  brms_family <- ifelse(family_of_regression == "multinomial" &
                          length(unique(regdata[["response_recoded"]]) > 2), "categorical",
                        ifelse(family_of_regression == "binary",  "bernoulli", "gaussian"))
  mixed_model <- brm(formula = as.formula(formula), data = regdata, family = brms_family,
                     cores = numcores, chains = numcores, iter = numiter, warmup = 0.1*numiter,
                     prior = type_of_prior,
                     control = list(max_treedepth = treedepth, adapt_delta = adaptdelta))
  datetoday <- format(x = Sys.time(), format = "%m%d%Y")
  saveRDS(object = mixed_model, file = sprintf("./%s_%s.RDS", name_of_output_file, datetoday))
}

## Check for MCMC convergence ##
## load the saved model
## plot(name_of_the_brms_model)
