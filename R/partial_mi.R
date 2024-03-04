#' Calculation of Partial Invariance by Parameter
#'
#' This function iterates over a constrained model
#' and examines potential areas for partial invariance.
#' The function returns a list of models and a dataframe
#' so that the user can decide which items may be relaxed
#' for partial invariance.
#'
#' @param saved_model A \code{lavaan} model with
#' at least one level of constraints across groups.
#' @param data Dataframe from the original saved model.
#' @param model Model from the original saved model.
#' @param group Grouping variable column from the
#' original saved model. Use a character vector.
#' @param group.equal A vector of names of the
#' constraints to test for multigroup analysis, in
#' the order that they should be tested. See \code{lavaan}
#' for possible options.
#' @param partial_step The level of partial invariance
#' you would like to test. You can use "loadings",
#' "regressions", "intercepts", "residuals", or "thresholds".
#' Note that for some models this may estimate more than you
#' want (i.e., the syntax for factor covariances is \code{~~}
#' as well as residuals) but these can be excluded from the
#' from the final dataframe.
#'
#' @return Models with each constraint relaxed individually
#' and a summary table of fit indices from each model.
#'
#' \item{models}{Saved \code{lavaan} models for each
#' parameter relaxed.}
#' \item{fit_table}{A dataframe of fit indices for
#' all models to investigate places for partial invariance.}
#'
#' @keywords multigroup cfa, sem, lavaan
#' @import lavaan dplyr ggplot2
#'
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
#' saved_mgcfa <- mgcfa(model = HS.model,
#'  data = HolzingerSwineford1939,
#'  group = "sex",
#'  group.equal = c("loadings", "intercepts", "residuals"),
#'  meanstructure = TRUE)
#'
#' saved_mi <- partial_mi(saved_model =
#'   saved_mgcfa$invariance_models$model.loadings,
#'   data = HolzingerSwineford1939,
#'   model = HS.model,
#'   group = "sex",
#'   group.equal = "loadings",
#'   partial_step = "loadings")
#'
#'   saved_mi$fit_table
#'   saved_mi$models
#'
#' @rdname partial_mi
#' @export

partial_mi <- function(saved_model,
                       data,
                       model,
                       group,
                       group.equal,
                       partial_step){


  # Deal with missing information  ------------------------------------------

  # intercepts: the intercepts of the observed variables ~1
  # means: the intercepts/means of the latent variables ~1
  # residuals: the residual variances of the observed variables ~~
  # residual.covariances: the residual covariances of the observed variables ~~
  # lv.variances: the (residual) variances of the latent variables ~~
  # lv.covariances: the (residual) covariances of the latent varibles ~~
  # regressions: all regression coefficients in the model ~
  # loadings: loadings on factors =~
  # thresholds: when using ordered models |

  if(missing(saved_model) | missing(partial_step) |
     missing(group.equal) | missing(data) | missing(model) |
     missing(group)) {
    stop("You must define all parameters!")
  }

  # deal with tibbles
  data <- as.data.frame(data)

  # get partial_step
  op_filter <- ifelse(
    partial_step == "loadings", "=~", ifelse(
      partial_step == "regressions", "~", ifelse(
        partial_step == "intercepts", "~1", ifelse(
          partial_step == "residuals", "~~", ifelse(
            partial_step == "thresholds", "|", NA
          )
        )
      )
    )
  )

  # first get the data to change
  partial.values <- tidy(saved_model) %>%
    filter(op == op_filter) %>%
    pull(term) %>%
    unique()

  # run the analysis
  partial.save <- list()
  fit.save <- list()

  for (i in 1:length(partial.values)){
    partial.save[[partial.values[i]]] <- lavaan::update(object = saved_model,
                                                        data = data,
                                                        model = model,
                                                        group = group,
                                                        group.equal = group.equal,
                                                        group.partial = partial.values[i])
    fit.save[[partial.values[i]]] <- fitmeasures(partial.save[[partial.values[i]]])
  }

  fit_table <- bind_rows(fit.save) %>%
    mutate(free.parameter = partial.values) %>%
    relocate(free.parameter)

  # print out results/suggestions
  return(list("models" = partial.save,
         "fit_table" = fit_table))
}

