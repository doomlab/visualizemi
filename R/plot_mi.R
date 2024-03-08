#' Visualization for Measurement Invariance
#'
#' A function to create visualizations for
#' measurement invariance comparison between
#' two groups. Please note: this function graphs
#' two groups at a time - use the function several
#' times if you have more than two groups.
#'
#' @param data_coef A tidy dataframe of the coefficients
#' from \code{mgcfa} function in this package.
#' @param model_step The model step you would like to
#' plot using this function. Check the "Model" column
#' in the coefficients to pick. Include this as a
#' character vector.
#' @param item_name The name of the item/variable
#' you want to plot for comparison. Check out the
#' "term" column to know your options. Usually
#' observed variables from your model.
#' @param x_limits The limits of the latent variable
#' to plot on on the latent mean graph.
#' @param y_limits The limits of the observed variable
#' that you included in \code{item_name}.
#' @param conf.level The confidence interval level you
#' would like to graph.
#' @param model_results The saved summary from \code{cfa}
#' or \code{mgcfa} of the same model name you listed
#' in \code{model_step}.
#' @param lv_name The name of the latent variable you
#' would like to plot. Ensure this name is the same
#' as what is in the \code{data_coef}.
#' @param plot_groups An optional variable to
#' denote which two groups you would like to plot. Used
#' when more than two groups are included in the output.
#'
#' @return A list of graphs created to visualize
#' measurement invariance.
#'
#' \item{complete}{A ggplot2 object that includes all
#' visualization stacked together for publication.}
#' \item{intercept}{A ggplot2 object that includes the
#' intercept section of the overall graph.}
#' \item{mean}{A ggplot2 object that includes the
#' latent means section of the overall graph.}
#' \item{variance}{A ggplot2 object that includes the
#' variance/residuals section of the overall graph.}
#'
#' @keywords multigroup cfa, sem, lavaan
#' @import lavaan dplyr ggplot2 introdataviz
#' @importFrom cowplot plot_grid get_legend
#' @include globals.R
#'
#' @examples
#'
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
#'  saved_mgcfa <- mgcfa(model = HS.model,
#'   data = HolzingerSwineford1939,
#'   group = "sex",
#'   group.equal = c("loadings", "intercepts", "residuals"),
#'   meanstructure = TRUE)
#'
#' # saved_mi_plots <- plot_mi(data_coef = saved_mgcfa$model_coef,
#' #  model_step = "Configural", # which model
#' #  item_name = "x1", # name of observed item
#' #  x_limits = c(-1,1), # LV limits to graph
#' #  y_limits = c(min(HolzingerSwineford1939$x1),
#' #             max(HolzingerSwineford1939$x1)), # Y min and max in data
#' #  conf.level = .95, # what ci do you want
#' #  model_results = saved_mgcfa$model_configural, # what model results do you want
#' #  lv_name = "visual", # which latent is the observed variable on
#' #  plot_groups = NULL)
#' @include globals.R
#'
#' @rdname plot_mi
#' @export

plot_mi <- function(data_coef, # output from model_coef
                       model_step, # which model
                       item_name, # name of observed item
                       x_limits = c(-1,1), # LV limits to graph
                       y_limits, # Y min and max in data
                       conf.level = .95, # what ci do you want
                       model_results, # what model results do you want
                       lv_name, # which latent is the observed variable on
                       plot_groups = NULL
){


  # Deal with missing information  ------------------------------------------
  if(missing(data_coef) |
     missing(model_step) |
     missing(item_name) |
     missing(y_limits) |
     missing(model_results) |
     missing(lv_name)) {
    stop("You must define data coefficients, model step, item name,
         y limits, model results, and the name of the latent variable.")
  }

  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("introdataviz", quietly = TRUE)

  # calculate cutoff
  cutoff <- qt(p = (1-conf.level)/2,
               df = sum(unlist(model_results@Data@nobs)),
               lower.tail = F)

  # get group variable
  group_var <- model_results@Data@group
  group_labels <- model_results@Data@group.label

  # pick the first two if it's null
  if(is.null(plot_groups)){
    plot_groups <- group_labels[1:2]
  }

  # first get the data
  graph.data <- data_coef %>% # put in tidy coefficients
    filter(model == model_step) %>% # pick a model
    filter(grepl(paste0(item_name, " |", item_name, "$"), term)) %>%  # pick a question
    mutate(group = factor(group, levels = names(table(data_coef$group)),
                          labels = group_labels)) %>%
    filter(group %in% plot_groups)

  # make ribbon data y = slope*x + intercept for ci for slopes
  ribbondata <- bind_rows(
    data.frame(
      x = seq(from = x_limits[1] - 1,
              to = x_limits[2] + 1,
              by = .05),
      group = unique(graph.data$group)[1]
    ) %>%
      mutate(ymin = (graph.data %>% filter(op == "=~") %>%
                       slice_head() %>% pull(estimate) * x) -
               (cutoff*graph.data %>% filter(op == "=~") %>%
                  slice_head() %>% pull(std.error)) +
               graph.data %>% filter(op == "~1") %>%
               slice_head() %>% pull(estimate),
             ymax = (graph.data %>% filter(op == "=~") %>%
                       slice_head() %>% pull(estimate) * x) +
               (cutoff*graph.data %>% filter(op == "=~") %>%
                  slice_head() %>% pull(std.error)) +
               graph.data %>% filter(op == "~1") %>%
               slice_head() %>% pull(estimate)),
    data.frame(
      x = seq(from = x_limits[1] - 1,
              to = x_limits[2] + 1,
              by = .05),
      group = unique(graph.data$group)[2]
    ) %>%
      mutate(ymin = (graph.data %>% filter(op == "=~") %>%
                       slice_tail() %>% pull(estimate) * x) -
               (cutoff*graph.data %>% filter(op == "=~") %>%
                  slice_tail() %>% pull(std.error)) +
               graph.data %>% filter(op == "~1") %>%
               slice_tail() %>% pull(estimate),
             ymax = (graph.data %>% filter(op == "=~") %>%
                       slice_tail() %>% pull(estimate) * x) +
               (cutoff*graph.data %>% filter(op == "=~") %>%
                  slice_tail() %>% pull(std.error)) +
               graph.data %>% filter(op == "~1") %>%
               slice_tail() %>% pull(estimate))
  )

  # make point data to draw on the intercepts
  pointdata <- data.frame(
    x = c(0,0),
    y = graph.data %>% filter(op == "~1") %>% pull(estimate),
    group = graph.data %>% filter(op == "~1") %>% pull(group),
    ymin = graph.data %>% filter(op == "~1") %>% pull(estimate) -
      cutoff * graph.data %>% filter(op == "~1") %>% pull(std.error),
    ymax = graph.data %>% filter(op == "~1") %>% pull(estimate) +
      cutoff * graph.data %>% filter(op == "~1") %>% pull(std.error)
  )

  # make the line data to draw on the slopes
  linedata <- data.frame(
    slope = graph.data %>% filter(op == "=~") %>% pull(estimate),
    intercept = graph.data %>% filter(op == "~1") %>% pull(estimate),
    group = graph.data %>% filter(op == "~1") %>% pull(group)
  )

  # make the distributions for the residuals
  violindata <- data.frame(
    y = c(rnorm(n = 1000,
                mean = graph.data %>% filter(op == "~~") %>%
                  slice_head() %>% pull(estimate),
                sd = graph.data %>% filter(op == "~~") %>%
                  slice_head() %>% pull(std.error)),
          rnorm(n = 1000,
                mean = graph.data %>% filter(op == "~~") %>%
                  slice_tail() %>% pull(estimate),
                sd = graph.data %>% filter(op == "~~") %>%
                  slice_tail() %>% pull(std.error))),
    group = c(rep(graph.data %>% filter(op == "~~") %>%
                    slice_head() %>% pull(group), 1000),
              rep(graph.data %>% filter(op == "~~") %>%
                    slice_tail() %>% pull(group), 1000)),
    x = 1
  )

  # make the latent mean data for right panel
  latent_means <- lavPredict(model_results,
                             type = "lv",
                             label = TRUE,
                             assemble = TRUE,
                             append.data = TRUE)

  latent_means$lv <- latent_means[ , lv_name]
  latent_means$group <- latent_means[ , group_var]

  # make a plot of the variance
  variance_plot <-
    ggplot(violindata, aes(x = 1, y = y, color = group, fill = group)) +
    geom_split_violin() +
    theme_void() +
    theme(legend.position = "none") +
    stat_summary(fun = "mean",
                 geom = "crossbar",
                 width = 0.5,
                 colour = "black")

  # make the plot with intercepts and slopes
  intercept_plot <-
    ggplot() +
    # basic set up
    theme_classic() +
    xlab("Latent Variable") +
    ylab("Observed Variable") +
    coord_cartesian(xlim = x_limits, ylim = y_limits) +
    # plot the intercepts
    geom_point(data = pointdata,
               aes(x = x, y = y, color = group),
               inherit.aes = FALSE) +
    geom_errorbar(data = pointdata,
                  aes(x = x, ymin = ymin, ymax = ymax, color = group),
                  inherit.aes = FALSE, width = .10) +
    # plot the slopes
    geom_abline(data = linedata,
                aes(slope = slope, intercept = intercept, color = group)) +
    geom_ribbon(data = ribbondata,
                aes(x = x, ymin = ymin, ymax = ymax, fill = group),
                inherit.aes = FALSE, alpha = .2) +
    scale_color_discrete(name = "Group") +
    scale_fill_discrete(name = "Group") +
    geom_vline(xintercept = 0) +
    theme(axis.line.y = element_blank())

  # make the latent means plot
  mean_plot <- ggplot(latent_means, aes(x = lv, fill = group)) +
    geom_density(alpha = .2) +
    theme_classic() +
    xlab("Latent Variable") +
    ylab("Density") +
    geom_vline(data = latent_means %>% group_by(group) %>% summarize(mean = mean(lv)),
               aes(xintercept = mean, color = group)) +
    theme(legend.position = "none") +
    coord_cartesian(xlim = x_limits)

  y_range = abs(y_limits[2] - y_limits[1])

  # line up the two plots
  prow <- suppressWarnings(plot_grid(
    intercept_plot +
      ggtitle("Item Invariance") +
      theme(legend.position = "none") +
      annotation_custom(ggplotGrob(variance_plot),
                        xmin = .25, xmax = 1,
                        ymin = y_limits[1], ymax = y_limits[2]-y_range/1.8),
    mean_plot +
      ggtitle("Latent Mean Distribution") +
      theme(legend.position = "none"),
    align = 'vh',
    hjust = -1,
    nrow = 1
  ))

  # march 8, 2024 get_legend from cowplot
  # ggplot2 update broke the legends
  get_legend <- function(plot, legend = NULL) {

    gt <- ggplotGrob(plot)

    pattern <- "guide-box"
    if (!is.null(legend)) {
      pattern <- paste0(pattern, "-", legend)
    }

    indices <- grep(pattern, gt$layout$name)

    not_empty <- !vapply(
      gt$grobs[indices],
      inherits, what = "zeroGrob",
      FUN.VALUE = logical(1)
    )
    indices <- indices[not_empty]

    if (length(indices) > 0) {
      return(gt$grobs[[indices[1]]])
    }
    return(NULL)
  }

  # get the legend
  legend_b <- get_legend(
    intercept_plot +
      guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )

  # send out the plot
  together <- suppressWarnings(plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1)))

  return(list("complete" = together,
              "intercept" = intercept_plot,
              "mean" = mean_plot,
              "variance" = variance_plot))

}
