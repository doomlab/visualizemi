#' \eqn{d_s} for Between Subjects with Pooled SD Denominator
#'
#' This function displays d for two between subjects groups
#' and gives the central and non-central confidence interval
#' using the pooled standard deviation as the denominator. Taken
#' from MOTE for this package. Use MOTE for full options.
#'
#' To calculate \eqn{d_s}, mean two is subtracted from mean one and divided
#' by the pooled standard deviation.
#' \deqn{d_s = \frac{M_1 - M_2}{S_{pooled}}}
#'
#'
#' @param m1 mean group one
#' @param m2 mean group two
#' @param sd1 standard deviation group one
#' @param sd2 standard deviation group two
#' @param n1 sample size group one
#' @param n2 sample size group two
#' @param a significance level
#'
#' @return Provides the effect size (Cohen's *d*) with associated
#' central and non-central confidence intervals,
#' the *t*-statistic, the confidence intervals associated with the means of
#' each group, as well as the standard deviations and standard errors
#' of the means for each group. The one-tailed confidence interval
#' is also included for sensitivity analyses.
#'
#' \item{d}{effect size}
#' \item{dlow}{lower level confidence interval of d value}
#' \item{dhigh}{upper level confidence interval of d value}
#'
#' @keywords effect size, independent t, between-subjects, pooled
#' standard deviation, pooled sd
#' @import methods
#' @importFrom stats nlm optimize pt qt rnorm sd
#'
#' @examples
#' calculate_d(m1 = 4, m2 = 3,
#' sd1 = 1, sd2 = 1,
#' n1 = 100, n2 = 100,
#' a = .05)
#'
#' @export

calculate_d <- function(m1 = NULL, m2 = NULL,
                         sd1 = NULL, sd2 = NULL,
                         n1 = NULL, n2 = NULL,
                         a = .05) {

  if (a < 0 || a > 1) {stop("Alpha should be between 0 and 1.")}

    # deal with missing values
    if (missing(m1)){stop("Be sure to include m1 for the first mean.")}
    if (missing(m2)){stop("Be sure to include m2 for the second mean.")}
    if (missing(sd1)){stop("Be sure to include sd1 for the first mean.")}
    if (missing(sd2)){stop("Be sure to include sd2 for the second mean.")}
    if (missing(n1)){stop("Be sure to include the sample size n1 for the first group.")}
    if (missing(n2)){stop("Be sure to include the sample size n2 for the second group.")}

    # calculate d
    spooled <- sqrt( ((n1 - 1) * sd1 ^ 2 + (n2 - 1) * sd2 ^ 2) / (n1 + n2 - 2))
    d <- (m1 - m2) / spooled

    # calculate t
    se1 <- sd1 / sqrt(n1)
    se2 <- sd2 / sqrt(n2)
    sepooled <- sqrt((spooled ^ 2 / n1 + spooled ^ 2 / n2))
    t <- (m1 - m2) / sepooled

  # calculate noncentral ci
  ncpboth <- noncentral_t(t, (n1 - 1 + n2 - 1), conf.level = (1 - a), sup.int.warns = TRUE)
  dlow <- ncpboth$Lower.Limit / sqrt(((n1 * n2) / (n1 + n2)))
  dhigh <- ncpboth$Upper.Limit / sqrt(((n1 * n2) / (n1 + n2)))

  output = list("d" = d, #d stats
                "dlow" = dlow,
                "dhigh" = dhigh)

  return(output)
}

#' @rdname calculate_d
#' @export
