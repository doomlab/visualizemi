## ---- include = FALSE----------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------
library(visualizemi)
library(lavaan)
library(knitr)
library(ggplot2)
library(introdataviz)
library(ggridges)

## ------------------------------------------------
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

                    # cfa model
saved_mgcfa <- mgcfa(model = HS.model,
                     # dataset in data frame 
                     data = HolzingerSwineford1939, 
                     # grouping variable column 
                     group = "sex",
                     # lavaan syntax for group constraints
                     group.equal = c("loadings", "intercepts", "residuals"), 
                     # any other lavaan cfa arguments
                     meanstructure = TRUE)

# note you can also include sample.nobs, sample.cov, and sample.mean if 
# you only have the correlation or covariance matrices

## ------------------------------------------------
kable(head(saved_mgcfa$model_coef))

kable(saved_mgcfa$model_fit)

# overall
saved_mgcfa$model_overall

# groups
# saved_mgcfa$group_models$model.1
# saved_mgcfa$group_models$model.2

# configural
# saved_mgcfa$model_configural

# constraints
# saved_mgcfa$invariance_models$model.loadings
# saved_mgcfa$invariance_models$model.intercepts
# saved_mgcfa$invariance_models$model.residuals

## ------------------------------------------------
                            # a saved model from mgcfa or any lavaan model
saved_partial <- partial_mi(saved_model = saved_mgcfa$invariance_models$model.loadings,
                            # dataframe of the original data
                            data = HolzingerSwineford1939, 
                            # model syntax from lavaan
                            model = HS.model, 
                            # group variable column
                            group = "sex",
                            # the equality constraints you have in the model
                            group.equal = c("loadings"),
                            # which step to test the partial invariance on
                            partial_step = "loadings")
# note that the group.equal and partial_step are not always the same 

kable(saved_partial$fit_table)

# saved_partial$models$`visual =~ x1`
# saved_partial$models$`visual =~ x2`
# and so on

## ------------------------------------------------
saved_mgcfa.partial <- mgcfa(model = HS.model,
                     # dataset in data frame 
                     data = HolzingerSwineford1939, 
                     # grouping variable column 
                     group = "sex",
                     # lavaan syntax for group constraints
                     group.equal = c("loadings", "intercepts", "residuals"), 
                     # any other lavaan cfa arguments
                     meanstructure = TRUE,
                     group.partial = c("speed =~ x9")
)

kable(saved_mgcfa.partial$model_fit)

## ------------------------------------------------
                          # a table of tidy coefficients, use broom::tidy 
                          # and create "model" column if you don't use mgcfa
                          # be sure to use the partial model or one with out 
                          # constraints or this graph will be boring
saved_mi_plots <- plot_mi(data_coef = saved_mgcfa.partial$model_coef,
                          # which model do you want to plot from model column
                          model_step = "loadings", 
                          # name of observed item
                          item_name = "x9", 
                          # LV limits to graph
                          x_limits = c(-1,1), 
                          # Y min and max in data
                          y_limits = c(min(HolzingerSwineford1939$x9),
                                       max(HolzingerSwineford1939$x9)), 
                          # what ci do you want
                          conf.level = .95, 
                          # what model results do you want
                          model_results = saved_mgcfa.partial$invariance_models$model.loadings, 
                          # which latent is the observed variable on
                          # important for cross-loaded variables
                          lv_name = "speed", 
                          # if you have more than two groups, which two do you want
                          plot_groups = NULL)

saved_mi_plots$complete

## ------------------------------------------------
saved_boot_model <- bootstrap_rr(
  # saved configural model to start at
  saved_configural = saved_mgcfa$model_configural,
  # dataset for the analysis
  data = HolzingerSwineford1939,
  # model lavaan syntax
  model = HS.model,
  # group variable in the dataset
  group = "sex", 
  # number of bootstraps
  # this is set to a low number to compile quickly for cran
  nboot = 50,
  # name of the fit measure you want to use, make sure it's lavaan
  invariance_index = "cfi",
  # rule for the difference in fit indices
  invariance_rule = .01,
  # what order of steps do you want to test? 
  group.equal = c("loadings", "intercepts", "residuals")
)

kable(saved_boot_model)

## ------------------------------------------------
saved_boot_partial <- bootstrap_partial(
  # the model you want to test 
  # use the model before your invariant one 
  # similar set up to partial_mi
  saved_model = saved_mgcfa.partial$model_configural,
  # the dataframe
  data = HolzingerSwineford1939,
  # the original model syntax
  model = HS.model,
  # the grouping variable column
  group = "sex",
  # run more, but this package vignette needs to knit fast
  nboot = 50,
  # what index are you using for invariance?
  # match this to lavaan's name under fitmeasures()
  invariance_index = "cfi",
  # what rule are you using?
  invariance_rule = .01,
  # what are we comparing against? 
  invariance_compare = unname(fitmeasures(saved_mgcfa.partial$model_configural, "cfi")),
  # which step you want to estimate effect size for
  partial_step = c("loadings"), 
  # which parameters you want to hold constrained
  group.equal = c("loadings")
  )

kable(head(saved_boot_partial$boot_DF))

kable(saved_boot_partial$boot_summary)

kable(saved_boot_partial$boot_effects)

saved_boot_partial$invariance_plot

saved_boot_partial$effect_invariance_plot

saved_boot_partial$density_plot

