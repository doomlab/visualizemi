#' Bootstrapping Measurement Invariance
#'
#' A function to help bootstrap expected rates of
#' "replication" for measurement invariance as compared
#' to a random assignment of groups. This function is
#' mostly designed for estimation for registered reports
#' but can also be used after a model is completed to
#' estimate rates of expected differences in the
#' future.
#'
#' @param saved_configural A saved \code{lavaan} model
#' of the "configural" model that includes the grouping
#' variable but no other equality constraints.
#' @param data The dataframe for the estimation
#' @param model The original model lavaan syntax
#' @param group The grouping variable column as a character
#' @param nboot The number of bootstraps you would like
#' to calculate. Please note: large models with many parameters
#' will run slowly depending on your computer.
#' @param invariance_index The name of the fix index you
#' want to use for invariance testing. For example,
#' "cfi", or "rmsea" - must be a character vector and
#' part of the \code{fitmeasures} provided by
#' \code{lavaan}.
#' @param invariance_rule The difference between your
#' previous model fit index and the new equality model
#' index that you would accept as invariant. Must be a
#' numeric value.
#' @param group.equal The equality constraints you would
#' like to impose as a vector, in order. This argument
#' is the same as \code{lavaan} syntax. For example,
#' c("loadings", "intercepts", "residuals").
#'
#' @return A dataframe of the bootstrapped results.
#' This set includes the proportion of invariant tests
#' for bootstrapping and random results and the
#' effect size h of this difference.
#'
#' @keywords multigroup cfa, sem, lavaan
#' @import lavaan dplyr
#' @importFrom tidyr pivot_longer
#' @include globals.R
#'
#' @examples
#' HS.model <- ' visual  =~ x1 + x2 + x3
#' textual =~ x4 + x5 + x6
#' speed   =~ x7 + x8 + x9 '
#'
#' library(lavaan)
#'
#' data("HolzingerSwineford1939")
#'
#' saved_model <- cfa(HS.model, data = HolzingerSwineford1939,
#'  meanstructure = TRUE,
#'  group = "sex",
#'  group.equal = c("loadings"))
#'
#' saved_configural <- cfa(HS.model, data = HolzingerSwineford1939,
#'  meanstructure = TRUE,
#'  group = "sex")
#'
#' # saved_boot <- bootstrap_rr(
#' #  saved_configural = saved_configural,
#' #  data = HolzingerSwineford1939,
#' #  model = HS.model,
#' #  group = "sex",
#' #  nboot = 100, # don't do this - this is to make it run fast
#' #  invariance_index = "cfi",
#' #  invariance_rule = .01,
#' #  group.equal = c("loadings", "intercepts", "residuals"))
#'
#' #  saved_boot
#'
#' @rdname bootstrap_rr
#' @export

bootstrap_rr <- function(saved_configural,
                         data,
                         model,
                         group,
                         nboot = 1000,
                         invariance_index,
                         invariance_rule,
                         group.equal){


  tol = 1e-5
  # Deal with missing information  ------------------------------------------
  if(missing(saved_configural)){stop("You must include the saved model.")}
  if(saved_configural@Data@data.type != "full"){stop("You must have full data for this function.")}

  # get the data
  DF <- data

  # try catch
  # test the model
  test_model_configural <- function(temp.fit, temp.DF, group, model) {
    tryCatch(
      {
          temp.partials <- update(temp.fit,
                                  model = model,
                                  data = temp.DF,
                                  group = group)
          return(temp.partials)

      }, warning = function(x){
        # just move on
        return(temp.partials <- NULL)
      }, error = function(x){
        # just move on
        return(temp.partials <- NULL)
      }
    )
  }

  test_model <- function(temp.fit, temp.DF, group.equal, group, model) {
    tryCatch(
      {

          temp.partials <- update(temp.fit,
                                  data = temp.DF,
                                  model,
                                  group.equal = group.equal,
                                  group = group)
          return(temp.partials)

      }, warning = function(x){
        # just move on
        return(temp.partials <- NULL)
      }, error = function(x){
        # just move on
        return(temp.partials <- NULL)
      }
    )
  }

  boot_results <- list()
  random_boot_results <- list()

  for (i in 1:nboot){
    # create bootstrapped data and random data
    # you have to have these values
    total_break <- 0
    while(total_break < 2){
      total_break <- 0

      random_group_results <- sample(DF[ , group], size = nrow(DF),
                                     replace = TRUE)

      temp.DF <- DF %>%
        slice_sample(n = nrow(DF), replace = TRUE) %>%
        mutate(random_group = random_group_results)

      # update the model
      temp.fit <- test_model_configural(saved_configural, temp.DF, group, model)
      random.fit <- test_model_configural(saved_configural, temp.DF, "random_group", model)

      # get rule comparison
      if(!is.null(temp.fit)){
        comparison <- fitmeasures(temp.fit, invariance_index)
        names(comparison) <- "configural"
        total_break <- total_break + 1
      }

      if(!is.null(random.fit)){
        random_comparison <- fitmeasures(random.fit, invariance_index)
        names(random_comparison) <- "random_configural"
        total_break <- total_break + 1
      }
    }

    # test if "invariant"
    # do steps they are interested in
    invariance.models <- list()
    random.invariance.models <- list()
    for (p in 1:length(group.equal)){

      # make it NULL
      temp.partials <- NULL
      random.partials <- NULL

      temp.partials <- test_model(temp.fit, temp.DF, group.equal[1:p], group, model)
      random.partials <- test_model(random.fit, temp.DF, group.equal[1:p], "random_group", model)

      # test if not null
      if(!is.null(temp.partials)){
        invariance.models[[p]] <- fitmeasures(temp.partials, invariance_index)
        names(invariance.models[[p]]) <- group.equal[p]
      } else {
        invariance.models[[p]] <- NA
        names(invariance.models[[p]]) <- group.equal[p]
      }


    # test if not null
    if(!is.null(random.partials)){
      random.invariance.models[[p]] <- fitmeasures(random.partials, invariance_index)
      names(random.invariance.models[[p]]) <- group.equal[p]
    } else {
      random.invariance.models[[p]] <- NA
      names(random.invariance.models[[p]]) <- group.equal[p]
    }

    }

    boot_results[[i]] <- c(comparison, unlist(invariance.models))
    random_boot_results[[i]] <- c(random_comparison, unlist(random.invariance.models))

    }

  boot_DF <- bind_rows(boot_results) %>%
    mutate(model_number = 1:nboot) %>%
    pivot_longer(cols = -model_number,
                 names_to = "model",
                 values_to = "fit_index") %>%
    group_by(model_number) %>%
    mutate(model_difference = lag(fit_index) - fit_index) %>%
    filter(!is.na(model_difference)) %>%
    mutate(invariant = model_difference <= (invariance_rule+tol)) %>%
    filter(!invariant) %>%
    slice_head() %>%
    group_by(model) %>%
    summarize(non_invariant = n()/nboot) %>%
    full_join(
      bind_rows(random_boot_results) %>%
        mutate(model_number = 1:nboot) %>%
        pivot_longer(cols = -model_number,
                     names_to = "model",
                     values_to = "fit_index") %>%
        group_by(model_number) %>%
        mutate(model_difference = lag(fit_index) - fit_index) %>%
        filter(!is.na(model_difference)) %>%
        mutate(invariant = model_difference <= (invariance_rule+tol)) %>%
        filter(!invariant) %>%
        slice_head() %>%
        group_by(model) %>%
        summarize(random_non_invariant = n()/nboot),
      by = "model"
    ) %>%
    mutate(non_invariant = ifelse(is.na(non_invariant), 0, non_invariant),
           random_non_invariant = ifelse(is.na(random_non_invariant), 0, random_non_invariant),
           h_nmi = 2*(asin(sqrt(non_invariant))-asin(sqrt(random_non_invariant))),
           h_mi = 2*(asin(sqrt(1-non_invariant))-asin(sqrt(1-random_non_invariant))),
           h_nmi_p = h_nmi / pi,
           h_mi_p = h_mi / pi)

  # print out results/suggestions
  return(boot_DF)
}

