
#------------------------------------------------
# if NULL then replace with chosen value, otherwise keep original value
#' @noRd
define_default <- function(x, default) {
  if (is.null(x)) {
    return(default)
  } else {
    return(x)
  }
}

#------------------------------------------------
# if a single value is provided then expand to a vector of length n
#' @noRd
force_vector <- function(x, n) {
  if (length(x) == 1) {
    return(rep(x,n))
  } else {
    return(x)
  }
}

#------------------------------------------------
# calculate midpoints of a vector
#' @noRd
midpoints <- function(x) {
  return((x[-1] + x[-length(x)])/2)
}

# -----------------------------------
# takes matrix as input, converts to list format for use within Rcpp code
#' @noRd
matrix_to_rcpp <- function(x) {
  return(split(x, f = 1:nrow(x)))
}

# -----------------------------------
# takes list format returned from Rcpp and converts to matrix
#' @noRd
rcpp_to_matrix <- function(x) {
  ret <- matrix(unlist(x), nrow = length(x), byrow = TRUE)
  return(ret)
}

# -----------------------------------
# takes list format returned from Rcpp and converts to three-dimensional array.
# Array indexing is in the same order as the underlying list, for example
# x[i,j,k] is equivalent to l[[i]][[j]][[k]]
#' @noRd
rcpp_to_array <- function(x) {
  ret <- array(unlist(x), dim = c(length(x[[1]][[1]]), length(x[[1]]), length(x)))
  ret <- aperm(ret, perm = c(3,2,1))
  return(ret)
}

#------------------------------------------------
# return 95% quantile
#' @importFrom stats quantile
#' @noRd
quantile_95 <- function(x) {
  ret <- quantile(x, probs = c(0.025, 0.5, 0.975))
  names(ret) <- c("Q2.5", "Q50", "Q97.5")
  return(ret)
}

#------------------------------------------------
# deal with possible NAs in the object x, which can be any object type that
# accommodates is.na() calls. For example, vectors and matrices are fine, but
# for lists it will only evaluate at the top level.
# Includes the following return options:
#   1 = stop() if all NA, otherwise return FALSE
#   2 = return TRUE if all NA, otherwise return FALSE
#   3 = stop() if any NA, otherwise return FALSE
#   4 = return TRUE if any NA, otherwise return FALSE
#' @noRd
process_NAs <- function(x, na_option = 3) {
  
  if (all(is.na(x))) {
    if (na_option == 1) {
      stop("x cannot contain all NAs")
    } else if (na_option == 2) {
      return(TRUE)
    }
  }
  if (any(is.na(x))) {
    if (na_option == 3) {
      stop("x cannot contain any NAs")
    } else if (na_option == 4) {
      return(TRUE)
    }
  }
  return(FALSE)
}

#------------------------------------------------
# sum values in log space, i.e. do log(sum(exp(x))), while avoiding some
# underflow/overflow issues. Robust to large +ve and large -ve values, including
# allowing -Inf values as long as there is at least one finite value in x.
# Includes options for dealing with NA values - see process_NAs().
#' @noRd
log_sum <- function(x, na_option = 3) {
  
  # deal with NAs
  if (process_NAs(x, na_option = na_option)) {
    return(NA)
  }
  
  # strip NAs if present
  x <- x[!is.na(x)]
  
  # deal with case that x a single value
  if (length(x) == 1) {
    return(x)
  }
  
  # deal with infinite edge cases
  if (any(!is.finite(x) & (x > 0))) {
    return(Inf)
  }
  if (all(!is.finite(x) & (x < 0))) {
    return(-Inf)
  }
  
  # sum in log space while allowing large +ve and -ve values
  x_max <- max(x)
  ret <- x_max + log(sum(exp(x - x_max)))
  return(ret)
}

#------------------------------------------------
# alternative version of log_sum() that allows for small values of x, along with
# large -ve values used in combination with other values. For example, the
# following values of x are dealt with more accurately using log_sum2():
# x <- c(0, -40)  # requires negative_x_correction
# x <- c(1e-20, -46)  # requires small_x_correction and negative_x_correction
#' @noRd
log_sum2 <- function(x, na_option = 3) {
  
  # deal with NAs
  if (process_NAs(x, na_option = na_option)) {
    return(NA)
  }
  
  # strip NAs if present
  x <- x[!is.na(x)]
  
  # deal with case that x a single value
  if (length(x) == 1) {
    return(x)
  }
  
  # deal with infinite edge cases
  if (any(!is.finite(x) & (x > 0))) {
    return(Inf)
  }
  if (all(!is.finite(x) & (x < 0))) {
    return(-Inf)
  }
  
  # if all x are very close to 0 then chain together approximations:
  # e^x = 1 + x
  # log(n + x) = log(n) + x / n
  # to move these values outside the log function
  if (all((x >= -1e-11) & (x <= 1e-11))) {
    ret <- log(length(x)) + mean(x)
    return(ret)
  }
  
  # if any x are close to 0 then use approximation above for these values only.
  # Results in a modified version of x and a correction factor to be added to
  # final result
  small_x_correction <- 0
  if (any((x >= -1e-11) & (x <= 1e-11))) {
    w_small <- which((x >= -1e-11) & (x <= 1e-11))
    small_x_correction <- sum(x[w_small]) / (length(w_small) + sum(exp(x[-w_small])))
    x <- c(log(length(w_small)), x[-w_small])
  }
  
  # rescale x to be relative to largest value
  m <- max(x)
  x <- x - m
  
  # if any x are large -ve values then use following approximation:
  # log(a + sum(exp(x))) = log(a) + sum(exp(x)) / a
  # Results in a modified version of x and a correction factor to be added to
  # final result
  negative_x_correction <- 0
  if (any(x < -23)) {
    w_negative <- which(x < -23)
    negative_x_correction <-  sum(exp(x[w_negative])) / sum(exp(x[-w_negative]))
    x <- x[-w_negative]
  }
  
  # use naive formula on remaining terms and add in corrections
  ret <- m + log(sum(exp(x))) + negative_x_correction + small_x_correction
  return(ret)
}

#------------------------------------------------
# take the mean of values in log space, i.e. do log(mean(exp(x))), while
# avoiding some underflow/overflow issues. Robust to large +ve and large -ve
# values, including allowing -Inf values as long as there is at least one finite
# value in x. Includes options for dealing with NA values - see process_NAs().
#' @noRd
log_mean <- function(x, na_option = 3) {
  
  # deal with NAs
  if (process_NAs(x, na_option = na_option)) {
    return(NA)
  }
  
  # strip NAs if present
  x <- x[!is.na(x)]
  
  # deal with case that x a single value
  if (length(x) == 1) {
    return(x)
  }
  
  # deal with infinite edge cases
  if (any(!is.finite(x) & (x > 0))) {
    return(Inf)
  }
  if (all(!is.finite(x) & (x < 0))) {
    return(-Inf)
  }
  
  # take mean in log space while allowing large +ve and -ve values
  x_max <- max(x)
  ret <- x_max + log(sum(exp(x - x_max))) - log(length(x))
  return(ret)
}

#------------------------------------------------
# alternative version of log_mean() that allows for small values of x, along with
# large -ve values used in combination with other values. For example, the
# following values of x are dealt with more accurately using log_mean2():
# x <- rep(1e-20, 3)
#' @noRd
log_mean2 <- function(x, na_option = 3) {
  
  # deal with NAs
  if (process_NAs(x, na_option = na_option)) {
    return(NA)
  }
  
  # strip NAs if present
  x <- x[!is.na(x)]
  
  # deal with case that x a single value
  n_x <- length(x)
  if (n_x == 1) {
    return(x)
  }
  
  # deal with infinite edge cases
  if (any(!is.finite(x) & (x > 0))) {
    return(Inf)
  }
  if (all(!is.finite(x) & (x < 0))) {
    return(-Inf)
  }
  
  # if all x are very close to 0 then chain together approximations:
  # e^x = 1 + x
  # log(n + x) = log(n) + x / n
  # leading to simple result of mean of x
  if (all((x >= -1e-11) & (x <= 1e-11))) {
    return(mean(x))
  }
  
  # rescale x to be relative to largest value
  m <- max(x)
  x <- x - m
  
  # use naive formula on remaining terms and add in corrections
  ret <- m + log(sum(exp(x))) - log(n_x)
  return(ret)
}

#------------------------------------------------
# evaluates log(1 - exp(x)) for a single -ve value of x, while avoiding some
# underflow/overflow issues. Robust to reasonable large -ve values of x, up to
# about -700, and robust to small -ve values of x. Handles -Inf correctly.
# Includes the following options for dealing with NA values:
#   1 = stop() if x is NA
#   2 = return NA if x is NA
#' @noRd
log_one_minus <- function(x, na_option = 1) {
  
  # check inputs
  assert_single_pos(-x, zero_allowed = TRUE,
                    message = "x must be a single negative value or zero")
  
  # deal with NAs
  if (is.na(x)) {
    if (na_option == 1) {
      stop("x cannot be NA")
    } else {
      return(NA)
    }
  }
  
  # for x large -ve value, use approximation:
  # log(1 - exp(x)) = -exp(x)
  if (x < -23) {
    return(-exp(x))
  }
  
  # for x small -ve value, use approximation:
  # exp(x) = 1 + x
  # log(1 - exp(x)) = log(-x)
  if (x > -1e-10) {
    return(log(-x))
  }
  
  return(log(1 - exp(x)))
}

#------------------------------------------------
# geweke_pvalue
# return p-value of Geweke's diagnostic convergence statistic, estimated from
# package coda
#' @importFrom stats pnorm
#' @importFrom coda geweke.diag
#' @noRd
geweke_pvalue <- function(x) {
  ret <- 2*pnorm(abs(coda::geweke.diag(x)$z), lower.tail=FALSE)
  return(ret)
}

#------------------------------------------------
# check that geweke p-value non-significant at alpha significance level on
# values x[1:n]
#' @importFrom coda mcmc
#' @noRd
test_convergence <- function(x, n, alpha = 0.01) {
  # fail if n = 1
  if (n == 1) {
    return(FALSE)
  }
  
  # fail if ESS too small
  ESS <- try(coda::effectiveSize(x[1:n]), silent = TRUE)
  if (class(ESS) == "try-error") {
    return(FALSE)
  }
  if (ESS < 10) {
    return(FALSE)
  }
  
  # fail if geweke p-value < threshold
  g <- geweke_pvalue(mcmc(x[1:n]))
  ret <- (g > alpha)
  if (is.na(ret)) {
    ret <- FALSE;
  }
  
  # return
  return(ret)
}

#------------------------------------------------
# update progress bar
# pb_list = list of progress bar objects
# name = name of this progress bar
# i = new value of bar
# max_i = max value of bar (close when reach this value)
# close = whether to close when reach end
#' @importFrom utils setTxtProgressBar
#' @noRd
update_progress <- function(pb_list, name, i, max_i, close = TRUE) {
  setTxtProgressBar(pb_list[[name]], i)
  if (i == max_i & close) {
    close(pb_list[[name]])
  }
}

# -----------------------------------
# ask user a yes/no question. Return TRUE/FALSE.
#' @noRd
user_yes_no <- function(x = "continue? (Y/N): ") {
  user_choice <- NA
  while (!user_choice %in% c("Y", "y" ,"N", "n")) {
    user_choice <- readline(x)
  }
  return(user_choice %in% c("Y", "y"))
}

#------------------------------------------------
# recursive function for converting nested list (of any depth) to single string.
# Each list is enclosed in "{}" each element within a list is separated with
# ",", and values in a vector are separated with " ".
#' @noRd
list_to_text <- function(x, s = NULL) {
  s <- paste0(s, "{")
  for (i in 1:length(x)) {
    if (i > 1) {
      s <- paste0(s, ",")
    }
    if (is.list(x[[i]])) {
      s <- list_to_text(x[[i]], s)
    } else {
      s <- paste0(s, paste(x[[i]], collapse = " "))
    }
  }
  s <- paste0(s, "}")
  return(s)
}

#------------------------------------------------
# write a list x to file. See list_to_text() for string format
#' @noRd
write_text_list <- function(x, file_path) {
  
  # convert list to single string
  s <- list_to_text(x)
  
  # write to file
  writeLines(s, file_path)
  
}
