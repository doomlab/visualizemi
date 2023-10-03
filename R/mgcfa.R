#' Multigroup Confirmatory Factor Analysis
#'
#' A function to calculate multigroup confirmatory
#' factor analysis given a list of model constraints.
#' Returns all models using \code{lavaan} package
#' and tidy dataframes of relevant information.
#'
#' @param model A character string representing the
#' overall model formatted in \code{lavaan} syntax
#' @param data A dataframe that includes a grouping
#' variable for multigroup analysis. Define \code{data}
#' and \code{group} OR \code{sample.cov},
#' \code{sample.mean}, and \code{sample.nobs}.
#' @param sample.cov A sample variance-covariance
#' matrix that includes row and/or column names that
#' match the observed variables in the model. For
#' multigroup analyses, you need to format this data
#' as a list of variance-covariance matrices with
#' names for each group.
#' @param sample.mean The mean vector for observed
#' variables. For multigroup analyses, you need to format
#' this data as a list with names for each group.
#' It is assumed they are in the same order as
#' \code{sample}.cov.
#' @param sample.nobs Number of observations for each
#' group. For multigroup analyses, you need to format
#' this data as a list with names for each group.
#' It is assumed they are in the same order as
#' \code{sample.cov}.
#' @param group A character name of the group column
#' included in the \code{data} argument.
#' @param group.equal A vector of names of the
#' constraints to test for multigroup analysis, in
#' the order that they should be tested. See \code{lavaan}
#' for possible options.
#' @param group.partial A vector of constraints to relax
#' for partial group models.
#' @param conf.level The confidence interval for model
#' coefficients.
#' @param ... Other arguments to be included to configure
#' the \code{cfa} function from \code{lavaan}. For
#' example, you can include arguments for ordered
#' models, clustering, sampling.weights, or estimators.
#'
#' @return A list of outputs from the multigroup
#' procedure.
#'
#' \item{model_coef}{A tidy dataframe of model
#' coefficients}
#' \item{model_fit}{A tidy dataframe of model fit indices}
#' \item{model.overall}{If a dataframe is provided,
#' the overall CFA model with all groups combined is included}
#' \item{model.GROUPS}{The individual group models}
#' \item{model.configural}{The configural model that
#' puts together both groups in the same model with
#' no constraints}
#' \item{model.CONSTRAINTS}{The model constraints provided
#' in group.equal are converted into models with those
#' increasing equality constraints}
#'
#' @keywords multigroup cfa, sem, lavaan
#' @import lavaan dplyr
#' @importFrom broom tidy glance
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
#' @rdname mgcfa
#' @export

mgcfa <- function(model,
                  data = NULL,
                  sample.cov = NULL,
                  sample.mean = NULL,
                  sample.nobs = NULL,
                  group = NULL,
                  group.equal,
                  group.partial = NULL,
                  conf.level = conf.level,
                  ...){

  # Deal with missing information -------------------------------------------

  total_data <- sum(c(is.null(data), is.null(group)))
  total_not_data <- sum(c(is.null(sample.cov),
                          is.null(sample.mean),
                          is.null(sample.nobs)))
  if (total_data != 2){
    if (total_not_data != 3){
      stop("You must define the data and group OR the sample.cov, sample.mean, and sample.nobs.")
    }
  }


  # Do this if data ---------------------------------------------------------

  if(!is.null(data)){

    # group information
    group_names <- unique(data[ , group])
    data$group <- data[ , group]

    # do the overall models
    model.overall <- cfa(model = model,
                         data = data,
                         ...)

    group_models <- list()
    for (i in 1:length(group_names)){

      group_models[[paste0("model.", group_names[i])]] <-
        cfa(model = model,
            data = subset(data, group == group_names[i]),
            ...)
    }

    model.configural <- cfa(model = model,
                            data = data,
                            group = group,
                            ...)

    equal_models <- list()
    # do steps they are interested in
    for (i in 1:length(group.equal)){

      equal_models[[paste0("model.", group.equal[i])]] <-
        cfa(model = model,
            data = data,
            group = group,
            group.equal = group.equal[1:i],
            group.partial = group.partial,
            ...)
    }

  } else {

    group_names <- names(sample.cov)
    names(sample.nobs) <- names(sample.mean) <- group_names

    # do the overall models
    model.overall <- NULL

    group_models <- list()
    for (i in 1:length(group_names)){

      group_models[[paste0("model.", group_names[i])]] <-
        cfa(model = model,
            sample.nobs = sample.nobs[[group_names[i]]],
            sample.mean = sample.mean[[group_names[i]]],
            sample.cov = sample.cov[[group_names[i]]],
            ...)
    }

    model.configural <- cfa(model = model,
                            data = data,
                            sample.nobs = sample.nobs,
                            sample.mean = sample.mean,
                            sample.cov = sample.cov,
                            ...)

    equal_models <- list()
    # do steps they are interested in
    for (i in 1:length(group.equal)){

      equal_models[[paste0("model.", group.equal[i])]] <-
        cfa(model = model,
            sample.nobs = sample.nobs,
            sample.mean = sample.mean,
            sample.cov = sample.cov,
            group.equal = group.equal[1:i],
            group.partial = group.partial,
            ...)
    }

  }

  # Put together models -----------------------------------------------------
  model_coef <- bind_rows(
    tidy(model.overall, conf.level = conf.level) %>% mutate(model = "Overall"))

  model_fit <- bind_rows(
    glance(model.overall) %>% mutate(model = "Overall"))

  for (i in 1:length(group_models)){
    model_coef <- bind_rows(
      model_coef,
      tidy(group_models[[i]], conf.level = conf.level) %>%
        mutate(model = paste0("Group ", group_names[i]))
    )

    model_fit <- bind_rows(
      model_fit,
      glance(group_models[[i]]) %>%
        mutate(model = paste0("Group ", group_names[i]))
    )
  }

  model_coef <- bind_rows(
    model_coef,
    tidy(model.configural, conf.level = conf.level) %>% mutate(model = "Configural")
  )

  model_fit <- bind_rows(
    model_fit,
    glance(model.configural) %>% mutate(model = "Configural")
  )

  for (i in 1:length(equal_models)){
    model_coef <- bind_rows(
      model_coef,
      tidy(equal_models[[i]], conf.level = conf.level) %>%
        mutate(model = group.equal[i])
    )

    model_fit <- bind_rows(
      model_fit,
      glance(equal_models[[i]]) %>%
        mutate(model = group.equal[i])
    )
  }

  return(list(
    "model_coef" = model_coef,
    "model_fit" = model_fit,
    "model_overall" = model.overall,
    "group_models" = group_models,
    "model_configural" = model.configural,
    "invariance_models" = equal_models
  ))

}
