#' @import purrr
#' @import stats
#' @import parallel
#' @import furrr
#' @import future
#' @import readr
#' @importFrom magrittr %>%
#' @importFrom utils capture.output
#' @details
#' Linear Regression, parallel Linear Regression, parallel General Linear regression with Little Bag of Bootstraps modified based on the original blblm package
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))



blbcoef <- function(fit) {
  coef(fit)
}


blbsigma <- function(fit) {
  sqrt(sum((fit$residuals^2)) / fit$df.residual)
}


split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}



lm_each_subsample <- function(formula, data, n, B) {
  replicate(B, lm_each_boot(formula, data, n), simplify = FALSE)
}



lm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  lm1(formula, data, freqs)
}



lm1 <- function(formula, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  fit <- lm(formula, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}


map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}


#### blblm parallel part ####
#' @title fit a blblm_parallel class
#' @description fit bootstrap LM models and utilize parallel computing
#' @param formula formula indicating y and x variables
#'
#' @param data dataframe
#' @param m number of splits
#' @param B number of bootstraps
#' @param cluster number of clusters to use for parallel computing. when = 1 will be running on single thread.
#' @return blblm_parallel object
#' @export
#' @examples
#' # proj = blblm::blblm_parallel(disp~mpg, mtcars, m = 3, B = 100, cluster = 4)
#'
blblm_parallel <- function(formula, data, m = 10, B = 5000, cluster = 1){
  data_list <- split_data(data, m)
  availCores <- parallel::detectCores()
  stopifnot(cluster %in% 1:availCores == TRUE)
  if (cluster >1){
    #cl <- makeCluster(cluster)
    suppressWarnings(plan(multiprocess, workers = cluster))
    estimates = future_map(data_list, function(z){
      lm_each_subsample(formula = formula, data = z, n = nrow(data), B = B)
    })
    #stopCluster(cl)

    res <- list(estimates = estimates, formula = formula)
    #closeAllConnections()
  }
  else{
    estimates = lapply(data_list, function(z){
      lm_each_subsample(formula = formula, data = z, n = nrow(data), B = B)
    })
    #stopCluster(cl)
    res <- list(estimates = estimates, formula = formula)
  }

  class(res) <- "blblm_parallel"
  invisible(res)
}

#' @title coefficients of blblm_parallel project
#' @description coefficients method of blblm_parallel class returning coefficients
#' @param object a fitted blblm_parallel class
#'
#' @param ... others
#'
#' @export
#' @method coef blblm_parallel
#' @examples
#' # proj = blblm::blblm_parallel(disp~mpg, mtcars, m = 3, B = 100, cluster = 4)
#' # coef(proj)
coef.blblm_parallel<- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
  #closeAllConnections()
}

#' @title MSE of blblm_parallel object
#' @description coefficients method of blblm_parallel class returning coefficients
#' @param object a fitted blblm_parallel class
#' @param confidence bool indicating if you want a confidence interval of sigma
#' @param level confidence level
#' @param cluster number of clusters to use for parallel computing. when = 1 will be running on single thread.
#' @param ... others
#'
#' @export
#' @method sigma blblm_parallel
#' @examples
#' # proj = blblm::blblm_parallel(disp~mpg, mtcars, m = 3, B = 100, cluster = 4) # not run
#' # sigma(proj, confidence = TRUE, level = 0.99, cluster = 1) # not run
sigma.blblm_parallel <- function(object, confidence = FALSE, level = 0.95, cluster=1, ...) {
  est <- object$estimates

  stopifnot(cluster %in% 1:parallel::detectCores())
  if (cluster == 1){
    sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
    if (confidence) {
      alpha <- 1 - level
      limits <- est %>%
        map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
        set_names(NULL)
      return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
    }
    else {
      return(sigma)
    }
  }
  else {
    suppressWarnings(future::plan(future::multiprocess, workers = cluster))
    sigma <- mean(future_map_dbl(est, ~ mean(map_dbl(., "sigma"))))
    if (confidence) {
      alpha <- 1 - level
      limits <- est %>%
        future_map_mean(~ quantile(future_map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
        set_names(NULL)
      return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
    }
    else {
      return(sigma)
    }
  }
  #closeAllConnections()
}

#' @title predict blblm_parallel
#' @description prediction using blblm models and parallel computing
#' @param object a blblm_parallel object
#' @param new_data a data frame
#' @param confidence boolean specify if you need confidence interval
#' @param level confidence level
#' @param cluster integer > 1 will uses future for parallel computing
#' @param ... others
#'
#' @export
#' @method predict blblm_parallel
#'
#' @examples
#' # proj1 = blblm::blblm_parallel(disp~mpg, mtcars, m = 3, B = 100, cluster = 4)
#' # predict(proj1, data.frame('mpg' = c(24,21)), confidence = TRUE, cluster = 4)
#'
predict.blblm_parallel <- function(object, new_data, confidence = FALSE, level = 0.95, cluster = 1, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms.formula(object$formula,data = new_data), "term.labels")), new_data)
  stopifnot(cluster %in% 1:parallel::detectCores())
  if (cluster > 1){
    suppressWarnings(future::plan(future::multiprocess, workers = cluster))
    if (confidence) {
      future_map_mean(est, ~ future_map_cbind(., ~ X %*% .$coef) %>%
                        apply(1, mean_lwr_upr, level = level) %>%
                        t())
    }
    else {
      future_map_mean(est, ~ future_map_cbind(., ~ X %*% .$coef) %>% rowMeans())
    }
    #closeAllConnections()
    future:::ClusterRegistry("stop")
  }
  else {
    if (confidence) {
      map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
                 apply(1, mean_lwr_upr, level = level) %>%
                 t())
    } else {
      map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
    }
  }

}

#' @title confidence intervals of blblm_parallel coefficients
#' @description coefficients method of blblm_parallel class returning coefficients
#' @param object a fitted blblm_parallel class
#' @param parm a vector of parameter names or NULL for all parameters in the model
#' @param level confidence level
#' @param cluster number of clusters to use for parallel computing. when = 1 will be running on single thread.
#' @param ... others
#'
#' @export
#' @method confint blblm_parallel
#' @examples
#' # proj = blblm::blblm_parallel(disp~mpg, mtcars, m = 3, B = 100, cluster = 2)
#' # confint(proj, parm = NULL, level = 0.99, cluster = 1)
confint.blblm_parallel <- function(object, parm = NULL, level = 0.95, cluster = 1, ...) {

  stopifnot(cluster %in% 1:parallel::detectCores())

  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates

  if (cluster >1) {
    suppressWarnings(plan(multiprocess, workers = cluster))
    out <- future_map_rbind(parm, function(p) {
      future_map_mean(est, ~ future_map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
    })
    future:::ClusterRegistry("stop")
  }
  else {
    #suppressWarnings(plan(multiprocess, workers = cluster))
    out <- map_rbind(parm, function(p) {
      map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
    })
  }
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  #closeAllConnections()
  dimnames(out)[[1]] <- parm
  out
}

#' @title print method of blblm_parallel
#' @description print blblm_parallel
#' @param x a blblm_parallel fitted class
#'
#' @param ... others
#'
#' @export
#' @method print blblm_parallel
#' @examples
#' # proj = blblm::blblm_parallel(disp~mpg, mtcars, m = 3, B = 100, cluster = 4)
#' # print(proj)
print.blblm_parallel <- function(x, ...) {
  cat("blblm_parallel model:", capture.output(x$formula))
  cat("\n")
}



########################################################################################


#### blbglm parallel part


#' @title fit a blbglm_parallel class
#' @description fit bootstrap GLM models and utilize parallel computing
#' @param formula formula indicating y and x variables
#'
#' @param data dataframe
#' @param m number of splits
#' @param B number of bootstraps
#' @param cluster number of clusters to use for parallel computing. when = 1 will be running on single thread.
#' @param family family of glm, logit model as default.
#' @export
#' @examples
#' #proj = blblm::blbglm_parallel(vs~mpg+disp+gear, mtcars, m = 3, B = 100, cluster = 2)
#'
blbglm_parallel <- function(formula, data, family = binomial(link = "logit"), m = 10, B = 5000, cluster = 1){
  data_list <- split_data(data, m)
  availCores <- parallel::detectCores()
  family1 = family
  stopifnot(cluster %in% 1:availCores == TRUE)
  if (cluster >1){
    suppressWarnings(plan(multiprocess, workers = cluster))
    estimates = future_map(data_list, function(z){
      glm_each_subsample(formula = formula, data = z, n = nrow(data), B = B, family=family1)
    })
  }
  else {
    estimates = lapply(data_list, function(z){
      glm_each_subsample(formula = formula, data = z, n = nrow(data), B = B, family=family1)
    })
  }
  #closeAllConnections()
  future:::ClusterRegistry("stop")
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blbglm_parallel"
  invisible(res)
}




# compute the estimates
glm_each_subsample <- function(formula, data, n, B, family) {
  replicate(B, glm_each_boot(formula, data, n, family), simplify = FALSE)
}


# compute the regression estimates for a blb dataset
glm_each_boot <- function(formula, data, n, family) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  glm1(formula, data, freqs, family)
}


# estimate the regression estimates based on given the number of repetitions
glm1 <- function(formula, data, freqs, family) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wrong variable from the global scope.
  environment(formula) <- environment()
  fit <- glm(formula, data, weights = freqs, family = family)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}

#' @title print a blbglm_parallel project
#' @description print a blbglm_parallel project
#' @param x a blbglm_parallel project
#' @param ... other arguments
#' @export
#' @method print blbglm_parallel
#' @examples
#' #proj = blblm::blbglm_parallel(vs~mpg+disp, mtcars, m = 3, B = 100, cluster = 2)
#' #print(proj)
print.blbglm_parallel <- function(x, ...) {
  cat("blbglm_parallel model:", capture.output(x$formula))
  cat("\n")
}


#' @title coefficients of blbglm_parallel project
#' @description coefficients method of blbglm_parallel class returning coefficients
#' @param object a fitted blbglm_parallel class
#'
#' @param ... others
#'
#' @export
#' @method coef blbglm_parallel
#' @examples
#' #proj = blblm::blbglm_parallel(vs~mpg+disp, mtcars, m = 3, B = 100, cluster = 2)
#' #coef(proj)
coef.blbglm_parallel<- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}

#' @title confint method of blbglm_parallel
#' @description confidence interval of blbglm_parallel
#' @param object a blblm_parallel fitted class
#'
#' @param parm specify parameter, otherwise all X variables in the formula
#' @param level confidence level
#' @param cluster for parallel computing if cluster > 1
#' @param ... others
#'
#' @export
#' @method confint blbglm_parallel
#' @examples
#' #proj = blblm::blbglm_parallel(vs~mpg+disp, mtcars, m = 3, B = 100, cluster = 2)
#' #confint(proj, level = 0.95, parm = NULL, cluster = 2)
confint.blbglm_parallel <- function(object, parm = NULL, level = 0.95, cluster = 1, ...) {

  stopifnot(cluster %in% 1:parallel::detectCores())

  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates

  if (cluster >1) {
    suppressWarnings(plan(multiprocess, workers = cluster))
    out <- future_map_rbind(parm, function(p) {
      future_map_mean(est, ~ future_map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
    })
    #closeAllConnections()
    future:::ClusterRegistry("stop")
  }
  else {
    #suppressWarnings(plan(multiprocess, workers = cluster))
    out <- map_rbind(parm, function(p) {
      map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
    })
  }
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' @title MSE of a blbglm_parallel object
#' @description coefficients method of blbglm_parallel class returning coefficients
#' @param object a fitted blblm_parallel class
#' @param confidence bool indicating if you want a confidence interval of sigma
#' @param level confidence level
#' @param cluster number of clusters to use for parallel computing. when = 1 will be running on single thread.
#' @param ... others
#'
#' @export
#' @method sigma blbglm_parallel
#' @examples
#' #proj = blblm::blbglm_parallel(vs~mpg+disp, mtcars, m = 3, B = 100, cluster = 2)
#' #sigma(proj, confidence = TRUE, level = 0.99, cluster = 2)
sigma.blbglm_parallel <- function(object, confidence = FALSE, level = 0.95, cluster=1, ...) {
  est <- object$estimates

  stopifnot(cluster %in% 1:parallel::detectCores())
  if (cluster == 1){
    sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
    if (confidence) {
      alpha <- 1 - level
      limits <- est %>%
        map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
        set_names(NULL)
      return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
    }
    else {
      return(sigma)
    }
  }
  else{
    suppressWarnings(future::plan(future::multiprocess, workers = cluster))
    sigma <- mean(future_map_dbl(est, ~ mean(map_dbl(., "sigma"))))
    if (confidence) {
      alpha <- 1 - level
      limits <- est %>%
        future_map_mean(~ quantile(future_map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
        set_names(NULL)
      #closeAllConnections()
      future:::ClusterRegistry("stop")
      return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
    }
    else {
      #closeAllConnections()
      return(sigma)
    }
  }
}


#### common parallel mapping functions

future_map_mean <- function(.x, .f, ...) {
  (furrr::future_map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

future_map_cbind <- function(.x, .f, ...) {
  furrr::future_map(.x, .f, ...) %>% reduce(cbind)
}

future_map_rbind <- function(.x, .f, ...) {
  furrr::future_map(.x, .f, ...) %>% reduce(rbind)
}

########################################################################
#' @title blblm_parallel on a list of csv files with same structure
#' @description fit bootstrap LM models and utilize parallel computing on a list of files. notice the files should be csv and all column names should exactly match.
#' @param formula formula indicating y and x variables
#'
#' @param file_names dataframe
#' @param B number of bootstraps
#' @param cluster number of clusters to use for parallel computing. when = 1 will be running on single thread. more than 2 cores is not compatible in packages by default
#'
#' @export
#' @examples
#' # listnames = c('~/blblm/files/data01.csv', '~/blblm/files/data02.csv',
#' # '~/blblm/files/data03.csv', '~/blblm/files/data04.csv')
#' # projec = blblm_parallel_list(mpg ~ cyl + hp, listnames, B = 100, cluster = 1)
blblm_parallel_list <- function(formula, file_names, B = 5000, cluster = 1){
  availCores <- parallel::detectCores()
  stopifnot(cluster %in% 1:availCores == TRUE)
  if (cluster >1){
    #cl <- makeCluster(cluster)
    suppressWarnings(plan(multiprocess, workers = cluster))
    estimates = future_map(file_names, function(z){
      suppressMessages({file01 =readr::read_csv(z)})
      lm_each_subsample(formula = formula, data = file01, n = nrow(file01), B = B)
      #close(file01)
    })
    #stopCluster(cl)
    #closeAllConnections()
    future:::ClusterRegistry("stop")
  }
  else{
    estimates = lapply(file_names, function(z){
    suppressMessages({file01 = readr::read_csv(z)})
    lm_each_subsample(formula = formula, data = file01, n = nrow(file01), B = B)
    })
  }
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm_parallel"
  invisible(res)
}