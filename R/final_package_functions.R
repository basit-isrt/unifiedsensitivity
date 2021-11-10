#' Obtain the point estimate intervals of counterfactual
#' means
#'
#' Obtain the partially identified point estimate intervals of counterfactual
#' means under the Odds Ratio Framework
#'
#' @param A treatment indicator
#' @param Y outcome
#' @param gamma sensitivity parameter (log odds ratio)
#' @param fitted.probs fitted (generalized) propensity scores
#'
#' @return The partially identified interval of the IPW estimator.
#'
#' @export
#'
get.extrema.OR <- function(A, Y, gamma = 0, fitted.probs) {
  fitted.logit <- qlogis(fitted.probs)

  eg <- exp(-fitted.logit)
  Y <- Y[A == 1]
  eg <- eg[A == 1]
  eg <- eg[order(-Y)]
  Y <- Y[order(-Y)]

  ## maximization
  num.each.low <- Y * (1 + exp(-gamma) * eg)
  num.each.up <- Y * (1 + exp(gamma) * eg)
  num <-
    c(0, cumsum(num.each.up)) + c(rev(cumsum(rev(num.each.low))), 0)

  den.each.low <- (1 + exp(-gamma) * eg)
  den.each.up <- (1 + exp(gamma) * eg)
  den <-
    c(0, cumsum(den.each.up)) + c(rev(cumsum(rev(den.each.low))), 0)

  maximum <- max(num / den)
  ## print(den[which.max(num/den)] / n)

  ## minimization
  num <-
    c(0, cumsum(num.each.low)) + c(rev(cumsum(rev(num.each.up))), 0)
  den <-
    c(0, cumsum(den.each.low)) + c(rev(cumsum(rev(den.each.up))), 0)
  minimum <- min(num / den)
  ## print(den[which.min(num/den)] / n)

  point_estimate <- c(minimum, maximum)
  names(point_estimate) <- c("lower", "upper")
  return(point_estimate)
}


#' @describeIn get.extrema.OR Obtain the partially identified point estimate
#' intervals of counterfactual means under the Risk Ratio Framework
#'
#' @inheritParams get.extrema.OR
#'
#' @return The partially identified interval of the IPW estimator.
#'
#' @export
#'
get.extrema.RR <- function(A, Y, gamma = 0, fitted.probs) {
  fitted.log <- log(fitted.probs)
  gamma1 <- fitted.log
  gamma1[fitted.log < -gamma] <- -gamma

  eg <- exp(-fitted.log)
  Y <- Y[A == 1]
  eg <- eg[A == 1]
  gamma1 <- gamma1[A == 1]
  eg <- eg[order(-Y)]
  Y <- Y[order(-Y)]
  gamma1 <- gamma1[order(-Y)]

  ## maximization
  num.each.low <- Y * (exp(gamma1) * eg)
  num.each.up <- Y * (exp(gamma) * eg)
  num <-
    c(0, cumsum(num.each.up)) + c(rev(cumsum(rev(num.each.low))), 0)

  den.each.low <- (exp(gamma1) * eg)
  den.each.up <- (exp(gamma) * eg)
  den <-
    c(0, cumsum(den.each.up)) + c(rev(cumsum(rev(den.each.low))), 0)

  maximum <- max(num / den)

  ## minimization
  num <-
    c(0, cumsum(num.each.low)) + c(rev(cumsum(rev(num.each.up))), 0)
  den <-
    c(0, cumsum(den.each.low)) + c(rev(cumsum(rev(den.each.up))), 0)

  minimum <- min(num / den)

  point_estimate <- c(minimum, maximum)
  names(point_estimate) <- c("lower", "upper")
  return(point_estimate)
}

#' Estimate the generalized propensity scores
#'
#' Estimate the generalized propensity scores using multinomial logistic
#' regression
#'
#' @param data A \code{\link[base:data.frame]{data.frame}} containing all variables required for the analysis
#' @param A_name the name of the treatment variable in the data frame
#' @param gps.formula an object of class \code{\link[stats:formula]{formula}} that specifies the regression model for GPS
#'
#' @return A matrix of the estimated GPS
#'
#' @export
#'

gps.estimate <- function(data, A_name, gps.formula) {
  if (typeof(gps.formula) != "character") {
    gps.formula = Reduce(paste0, deparse(gps.formula))
  }
  gps.formula = gsub(" ", "", gps.formula)

  fitted.probs <- nnet::multinom(formula = gps.formula,
                                 data = data,
                                 trace = F)$fitted.values

  ## Creating the propensity score matrix in the case of binary treatment
  if (dim(fitted.probs)[2] == 1) {
    fitted.probs <- cbind(1 - fitted.probs, fitted.probs)
  }


  return(fitted.probs)
}




#' Fit the outcome regression model
#'
#' @inheritParams gps.estimate
#' @param outreg.formula an object of class \code{\link[stats:formula]{formula}} that specifies the regression model for the outcome regression
#' @param outreg.family a description of the error distribution and link function to be used in the model. (See \code{\link[stats:family]{family}} for details of family functions.)
#' @return The fitted outcome regression model object
#'
#' @export
#'

outreg <-
  function(data, outreg.formula, outreg.family = "gaussian") {
    if (typeof(outreg.formula) != "character") {
      outreg.formula = Reduce(paste0, deparse(outreg.formula))
    }
    outreg.formula = gsub(" ", "", outreg.formula)

    outreg.model = glm(formula = outreg.formula,
                       family = outreg.family,
                       data = data)

    return(outreg.model)
  }


extrema.cf.mean <-
  function(data_dummies, A, Y, gamma = 0, fitted.probs, outreg.formula,
           outreg.family = "gaussian", AIPW = FALSE, method = "RR") {
    data_dummies[A == 1,]
    if (AIPW) {
      outreg.model <- outreg(data_dummies[A == 1,], outreg.formula,
                             outreg.family)
      Y.fitted <-
        predict(outreg.model, data_dummies, type = "response")
    } else {
      Y.fitted <- rep(0, length(Y))
    }

    if (method == "RR")
      out <-
        get.extrema.RR(A, Y - Y.fitted, gamma, fitted.probs) + mean(Y.fitted)
    else if (method == "OR")
      out <-
        get.extrema.OR(A, Y - Y.fitted, gamma, fitted.probs) + mean(Y.fitted)

    return(out)
  }

#' Sensitivity analysis for observational studies
#'
#' @inheritParams gps.estimate
#' @inheritParams outreg
#' @param Y_name the name of the outcome variable in the data frame
#' @param contrast the linear contrast to define the target estimand
#' @param gamma sensitivity parameter (log odds ratio or log risk ratio)
#' @param AIPW should the doubly robust AIPW estimator be implemented? (TRUE or FALSE)
#' @param method The sensitivity analysis framework to be used ("OR" or "RR")
#'
#' @import stats
#'
#' @return The partially identified interval of the target estimand
#'
#' @export
#'

extrema.os.unified <-
  function(data,
           A_name,
           Y_name,
           gps.formula,
           contrast,
           gamma = 0,
           AIPW = FALSE,
           method = "RR",
           outreg.formula = NULL,
           outreg.family = "gaussian") {
    Y <- data[, Y_name]
    A <- factor(data[, A_name])

    fitted.probs <- gps.estimate(data, A_name, gps.formula)

    treatment <- c(0, 0)
    control <- c(0, 0)

    dummies <- fastDummies::dummy_cols(data[, A_name])[, -1]
    colnames(dummies) <- paste0("A", seq_along(contrast))

    data_dummies <- cbind(data, dummies)

    for (i in seq_along(contrast)) {
      if (contrast[i] > 0) {
        cf.mean <-
          extrema.cf.mean(
            data_dummies,
            data_dummies[, paste0("A", i)],
            Y,
            gamma,
            fitted.probs[, i],
            outreg.formula,
            outreg.family,
            AIPW,
            method
          )

        treatment <- treatment + contrast[i] * cf.mean
        # get.extrema.mul(A = data_dummies[ , paste0("A", i)], Y = Y,
        #                 gamma = gamma, fitted.prob = fitted.probs[ , i])
      }

      else if (contrast[i] < 0) {
        cf.mean <-
          extrema.cf.mean(
            data_dummies,
            data_dummies[, paste0("A", i)],
            Y,
            gamma,
            fitted.probs[, i],
            outreg.formula,
            outreg.family,
            AIPW,
            method
          )
        control <- control + contrast[i] * cf.mean

      }
    }

    return(treatment + rev(control))
  }

#' @describeIn extrema.os.unified Obtain a (1 - alpha) confidence interval under
#' the specified sensitivity analysis framework
#'
#' @inheritParams extrema.os.unified
#' @param alpha Significance level
#' @param parallel Should parallel computing be used? (TRUE or FALSE)
#' @param B Number of Bootstrap resamples.
#'
#' @return A (1 - alpha) confidence interval
#'
#' @import parallel
#' @export
#'
#' @examples
#' # sensitivity analysis using the mtcars dataset
#' # Let mpg be the outcome and cyl be the treatment with three levels
#'
#' # sensitivity analysis for pairwise ATE between cylinder type 4 and 8
#'
#' # gamma = 0 (No unmeasured confounding)
#' ## point estimate interval
#' extrema.os.unified(data = mtcars, A_name = "cyl", Y_name = "mpg",
#'                    gps.formula = cyl ~ disp + hp + drat + wt, gamma = 0,
#'                    contrast = c(1, 0, -1), method = "RR")
#'
#' ## 90% Confidence interval
#' set.seed(100)
#' bootsens.os.unified(data = mtcars, A_name = "cyl", Y_name = "mpg",
#'                     gps.formula = cyl ~ disp + hp + drat + wt, gamma = 0,
#'                     contrast = c(1, 0, -1), parallel = TRUE, alpha = 0.1)
#'
#' # gamma = 0.5
#' ## point estimate interval
#' extrema.os.unified(data = mtcars, A_name = "cyl", Y_name = "mpg",
#'                    gps.formula = cyl ~ disp + hp + drat + wt, gamma = 0.5,
#'                    contrast = c(1, 0, -1), method = "RR")
#'
#' ## 90% Confidence interval
#' set.seed(100)
#' bootsens.os.unified(data = mtcars, A_name = "cyl", Y_name = "mpg",
#'                     gps.formula = cyl ~ disp + hp + drat + wt, gamma = 0.5,
#'                     contrast = c(1, 0, -1), parallel = TRUE, alpha = 0.1)
bootsens.os.unified <-
  function(data,
           A_name,
           Y_name,
           gps.formula,
           contrast,
           gamma = 0,
           AIPW = FALSE,
           method = "RR",
           outreg.formula = NULL,
           outreg.family = "gaussian",
           alpha = 0.1,
           parallel = FALSE,
           B = 1000) {
    no.cores <- ifelse(parallel, parallel::detectCores(), 1)

    n <- dim(data)[1]

    out <- parallel::mclapply(1:B, function(iter) {
      s <- sample(1:n, n, TRUE)

      res <-
        tryCatch(
          extrema.os.unified(
            data = data[s, ],
            A_name,
            Y_name,
            gps.formula,
            contrast,
            gamma,
            AIPW,
            method,
            outreg.formula,
            outreg.family
          ),
          error = function(e) {
            print(e)
          }
        )

      res
    },
    mc.cores = no.cores)
    out <- do.call(rbind, out)

    c(quantile(out[, 1], alpha / 2, na.rm = TRUE),
      quantile(out[, 2], 1 - alpha / 2, na.rm = TRUE))

  }
