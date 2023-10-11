#' Bootstrapping Partial Invariance
#'
#' A function to bootstrap a partially invariant
#' multigroup CFA model. This function bootstraps
#' a model with a focus on freeing one parameter
#' at a time for a specific partial invariance model.
#' For example, you can use this function to view all
#' differences in intercepts for observed variables
#' between two groups. This function returns the
#' bootstrapped values, a summary table, and two
#' invariance plots comparing bootstrapped data to
#' randomized group assignment.
#'
#' @param saved_model A saved lavaan model of the
#' step you would like to estimate partial
#' invariance for.
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
#' previous model fit index and the new partial model
#' index that you would accept as invariant. Must be a
#' numeric value.
#' @param invariance_compare The comparison fit index to
#' judge if a model is invariant. This value will be compared
#' to the calculated partial invariance models and the
#' difference will be calculated to compare to \code{invariance_rule}.
#' @param partial_step Which parameter you would like to
#' bootstrap partial invariance for. You can use "loadings",
#' "regressions", "intercepts", "residuals", or "thresholds".
#' Note that for some models this may estimate more than you
#' want (i.e., the syntax for factor covariances is \code{~~}
#' as well as residuals) but these can be excluded from the
#' from the final dataframe.
#' @param group.equal The equality constraints for this
#' bootstrapping - use this parameter in case you want to
#' estimate effect size for a specific partial step
#' but continue to hold all other things constrained
#'
#' @return A set of graphs and dataframes.
#'
#' \item{invariance_plot}{A ggplot2 object that visualizes
#' the bootstrapped and random estimates difference score
#' between groups. These are separated by "invariant" and
#' "non-invariant" models.}
#' \item{effect_invariance_plot}{A ggplot2 object that
#' visualizes the effect sizes of the differences of the
#' invariance plot in a forest plot style. }
#' \item{density_plot}{A ggplot2 object that
#' visualizes the raw data of the boot_DF comparing
#' groups and invariance results.}
#' \item{boot_DF}{A dataframe of the bootstrapped and
#' randomized results.}
#' \item{boot_summary}{A dataframe of the summary of the
#' bootstrapped and random results with effect sizes. Note:
#' these last two dataframes can be used to recreate the
#' visualizations in your own style. }
#' \item{boot_effects}{A dataframe of the summary of the
#' bootstrapped and random results with effect sizes for
#' replication rates. }
#'
#' @keywords multigroup cfa, sem, lavaan
#' @import lavaan dplyr ggplot2 broom ggridges
#' @importFrom tidyr pivot_wider
#' @importFrom broom tidy glance
#' @include globals.R
#'
#' @examples
#' library(lavaan)
#' HS.model <- ' visual  =~ x1 + x2 + x3
#' textual =~ x4 + x5 + x6
#' speed   =~ x7 + x8 + x9 '
#'
#' saved_model <- cfa(HS.model, data = HolzingerSwineford1939,
#'  meanstructure = TRUE,
#'  group = "sex",
#'  group.equal = c("loadings"))
#'
#' # not run to save load time
#' # saved_boot <- bootstrapped_partial(
#' #  saved_model = saved_model,
#' #  data = HolzingerSwineford1939,
#' #  model = HS.model,
#' #  group = "sex",
#' #  nboot = 100,
#' #  invariance_index = "cfi",
#' #  invariance_rule = .01,
#' #  invariance_compare = fitmeasures(saved_model, "cfi"),
#' #  partial_step = c("loadings"),
#' #  group.equal = c("loadings"))
#'
#' # saved_boot$boot_DF
#' # saved_boot$boot_summary
#' # saved_boot$invariance_plot
#' # saved_boot$effect_invariance_plot

#' @rdname bootstrapped_partial
#' @export

bootstrapped_partial <- function(saved_model,
                                 data,
                                 model,
                                 group,
                                 nboot = 1000,
                                 invariance_index,
                                 invariance_rule,
                                 invariance_compare,
                                 partial_step,
                                 group.equal){


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

  tol = 1e-5

  if(missing(saved_model)){stop("You must include a saved lavaan model.")}
  if(missing(invariance_index)){stop("You must includedan invariance_index")}
  if(missing(invariance_rule)){stop("You must include an invariance_rule.")}
  if(missing(invariance_compare)){stop("You must include an invariance_compare.")}
  if(missing(data)){stop("You must include a dataframe to analyze.")}
  if(missing(model)){stop("You must include the original model syntax.")}
  if(missing(group)){stop("You must include a group column as a character vector.")}
  if(missing(group.equal)){stop("You must include the equality constraints.")}

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

  # deal with busted models
  test_model <- function(saved_model, temp.DF, i, partial.values, group,
                         model, group.equal) {
    tryCatch(
      {
          temp.model <- lavaan::update(saved_model,
                                         data = temp.DF,
                                         group = group,
                                         group.partial = partial.values[i],
                                         model = model,
                                         group.equal = group.equal)
          return(temp.model)

      }, warning = function(x){
        # just move on
        return(temp.model <- NULL)
      }, error = function(x){
        # just move on
        return(temp.model <- NULL)
      }
    )
  }

  boot_results <- list()
  for (p in 1:nboot){

    # bootstrap the data
    DF <- data
    random_group_results <- sample(DF[ , group], size = nrow(DF),
                                   replace = TRUE)
    temp.DF <- DF %>%
      slice_sample(n = nrow(DF), replace = TRUE) %>%
      mutate(random_group = random_group_results)

    # loop over model and update with relaxed parameter
    temp.parameters <- list()
    for (i in 1:length(partial.values)){

      # make it NULL
      temp.model <- NULL
      random.model <- NULL

      temp.model <- test_model(saved_model,
                               temp.DF, i,
                               partial.values,
                               group,
                               model, group.equal)
      random.model <- test_model(saved_model,
                                 temp.DF, i,
                                 partial.values,
                                 "random_group",
                                 model, group.equal)

      # get estimates of relaxed parameters
      if(!is.null(temp.model) & !is.null(random.model)){

        temp.parameters[[i]] <- tidy(temp.model) %>%
          filter(op == op_filter) %>%
          filter(term == partial.values[i]) %>%
          select(term, std.all, group) %>%
          mutate(type = "boot") %>%
          bind_rows(
            tidy(random.model) %>%
              filter(op == op_filter) %>%
              filter(term == partial.values[i]) %>%
              select(term, std.all, group) %>%
              mutate(type = "random")
          ) %>%
          pivot_wider(id_cols = term,
                      names_from = c("type", "group"),
                      values_from = "std.all") %>%
          mutate(boot_fit = unname(fitmeasures(temp.model, invariance_index)),
                 random_fit = unname(fitmeasures(random.model, invariance_index)))
      }


    }

    boot_results[[p]] <- bind_rows(temp.parameters)

  }

  # put together results
  boot_results <- boot_results[lapply(boot_results, length)>0]
  boot_DF <- as.data.frame(bind_rows(boot_results))
  boot_DF$boot_difference <- unname(boot_DF[ , 2] - boot_DF[ , 3])
  boot_DF$random_difference <- unname(boot_DF[ , 4] - boot_DF[ , 5])
  boot_DF$boot_index_difference <- (invariance_compare - boot_DF$boot_fit) <= (invariance_rule+tol)
  boot_DF$random_index_difference <- (invariance_compare - boot_DF$random_fit) <= (invariance_rule+tol)

  boot_long <- bind_rows(
    boot_DF %>%
      select(term, boot_difference, boot_index_difference) %>%
      rename(difference = boot_difference,
             index_difference = boot_index_difference) %>%
      mutate(type = "Bootstrapped"),
    boot_DF %>%
      select(term, random_difference, random_index_difference) %>%
      rename(difference = random_difference,
             index_difference = random_index_difference) %>%
      mutate(type = "Random")
  )

  label_graph <- c(
      `TRUE` = "Invariant",
      `FALSE` = "Non-Invariant"
    )

  invariance_plot <-
    ggplot(boot_long, aes(term, difference, color = type)) +
    stat_summary(fun = mean,
                 geom = "point",
                 position = position_dodge(width = 0.90)) +
    stat_summary(fun.data = mean_cl_normal,
                 geom = "errorbar",
                 position = position_dodge(width = 0.90),
                 width = .2) +
    theme_classic() +
    xlab("Estimated Item") +
    ylab("Difference Score Between Groups") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "bottom") +
    facet_wrap(~index_difference, labeller = as_labeller(label_graph)) +
    scale_color_discrete(name = "Type of Estimate") +
    coord_flip()

  boot_summary <- boot_DF %>%
    group_by(term, boot_index_difference) %>%
    summarize(mean_boot_1 = mean(boot_1, na.rm = T),
              mean_boot_2 = mean(boot_2, na.rm = T),
              sd_boot_1 = sd(boot_1, na.rm = T),
              sd_boot_2 = sd(boot_2, na.rm = T),
              mean_boot_difference = mean(boot_difference, na.rm = T),
              sd_boot_difference = sd(boot_difference, na.rm = T),
              mean_boot_fit = mean(boot_fit, na.rm = T),
              sd_boot_fit = sd(boot_fit, na.rm = T),
              n_boot = n(),
              .groups = "keep") %>%
    rename(invariant = boot_index_difference) %>%
    left_join(
      boot_DF %>%
        group_by(term, random_index_difference) %>%
        summarize(mean_random_1 = mean(random_1, na.rm = T),
                  mean_random_2 = mean(random_2, na.rm = T),
                  sd_random_1 = sd(random_1, na.rm = T),
                  sd_random_2 = sd(random_2, na.rm = T),
                  mean_random_difference = mean(random_difference, na.rm = T),
                  sd_random_difference = sd(random_difference, na.rm = T),
                  mean_random_fit = mean(random_fit, na.rm = T),
                  sd_random_fit = sd(random_fit, na.rm = T),
                  n_random = n(),
                  .groups = "keep") %>%
        rename(invariant = random_index_difference),
      by = c("term", "invariant")
    )

  boot_effects <- boot_summary %>%
    ungroup() %>%
    select(term, invariant, n_boot, n_random) %>%
    filter(invariant == FALSE) %>%
    select(-invariant) %>%
    rename(non_invariant = n_boot ,
           random_non_invariant = n_random) %>%
    mutate(non_invariant = non_invariant / nboot,
           random_non_invariant = random_non_invariant / nboot,
           non_invariant = ifelse(is.na(non_invariant), 0, non_invariant),
           random_non_invariant = ifelse(is.na(random_non_invariant), 0, random_non_invariant),
           h_nmi = 2*(asin(sqrt(non_invariant))-asin(sqrt(random_non_invariant))),
           h_mi = 2*(asin(sqrt(1-non_invariant))-asin(sqrt(1-random_non_invariant))))


  for (i in 1:nrow(boot_summary)){

    if(!is.na(boot_summary$n_boot[i]) &
       boot_summary$n_boot[i] >= nboot*.10){
      temp <- calculate_d(m1 = boot_summary$mean_boot_1[i],
                          m2 = boot_summary$mean_boot_2[i],
                          sd1 = boot_summary$sd_boot_1[i],
                          sd2 = boot_summary$sd_boot_2[i],
                          n1 = boot_summary$n_boot[i],
                          n2 = boot_summary$n_boot[i])

      boot_summary$d_boot_low[i] <- temp$dlow
      boot_summary$d_boot[i] <- temp$d
      boot_summary$d_boot_high[i] <- temp$dhigh
    } else{
      boot_summary$d_boot_low[i] <- NA
      boot_summary$d_boot[i] <- NA
      boot_summary$d_boot_high[i] <- NA
    }

    if(!is.na(boot_summary$n_random[i]) &
       boot_summary$n_random[i] >= nboot*.10){
      temp <- calculate_d(m1 = boot_summary$mean_random_1[i],
                          m2 = boot_summary$mean_random_2[i],
                          sd1 = boot_summary$sd_random_1[i],
                          sd2 = boot_summary$sd_random_2[i],
                          n1 = boot_summary$n_random[i],
                          n2 = boot_summary$n_random[i])

      boot_summary$d_random_low[i] <- temp$dlow
      boot_summary$d_random[i] <- temp$d
      boot_summary$d_random_high[i] <- temp$dhigh
    } else{
      boot_summary$d_random_low[i] <- NA
      boot_summary$d_random[i] <- NA
      boot_summary$d_random_high[i] <- NA
    }

    }

  boot_summary_long <- bind_rows(
    boot_summary %>%
      select(term, invariant, n_boot,
             starts_with("d_boot")) %>%
      rename(n = n_boot,
             d_low = d_boot_low,
             d = d_boot,
             d_high = d_boot_high) %>%
      mutate(type = "Bootstrapped"),
    boot_summary %>%
      select(term, invariant, n_random,
             starts_with("d_random")) %>%
      rename(n = n_random,
             d_low = d_random_low,
             d = d_random,
             d_high = d_random_high) %>%
      mutate(type = "Random")
  )

  effect_invariance_plot <-
    ggplot(boot_summary_long, aes(term, d, color = type)) +
    theme_classic() +
    facet_wrap(~invariant, labeller = as_labeller(label_graph)) +
    scale_color_discrete(name = "Type of Estimate") +
    geom_point(aes(size = n)) +
    geom_errorbar(aes(ymin = d_low, ymax = d_high),
                  width = .2) +
    xlab("Estimated Item") +
    ylab("Group Effect Size") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "bottom") +
    geom_hline(yintercept = 0) +
    coord_flip()

    DF_long <-
    bind_rows(
      boot_DF %>%
        select(term, boot_index_difference, boot_1, boot_2) %>%
        pivot_longer(cols = c(boot_1, boot_2),
                     names_to = "group",
                     values_to = "estimate") %>%
        rename(invariant = boot_index_difference) %>%
        mutate(type = "Bootstrapped"),
      boot_DF %>%
        select(term, random_index_difference, random_1, random_2) %>%
        pivot_longer(cols = c(random_1, random_2),
                     names_to = "group",
                     values_to = "estimate") %>%
        rename(invariant = random_index_difference) %>%
        mutate(type = "Random")

    ) %>%
      mutate(group = gsub("boot_1|random_1", "Group 1", group),
             group = gsub("boot_2|random_2", "Group 2", group),
             invariant = factor(invariant,
                                levels = c("TRUE", "FALSE"),
                                labels = c("Invariant", "Non-Invariant")))

  density_plot <- ggplot(DF_long, aes(x = estimate, y = term, fill = group)) +
    geom_density_ridges(
      # jittered_points = TRUE, # useful if you want to see the dots
      alpha = 0.7, scale = 0.9
    ) +
    # or facet_grid(~invariant*type) for every parameter on an entire row
    facet_grid(invariant~type) +
    theme_classic() +
    xlab("Estimated Score") +
    ylab("Parameter") +
    scale_fill_discrete(name = "Group")

  # print out results/suggestions
  return(list(
    "invariance_plot" = invariance_plot,
    "effect_invariance_plot" = effect_invariance_plot,
    "density_plot" = density_plot,
    "boot_DF" = boot_DF,
    "boot_summary" = boot_summary,
    "boot_effects" = boot_effects
  ))
}

