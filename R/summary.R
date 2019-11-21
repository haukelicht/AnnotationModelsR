# summary helpers ----

#' \code{\link{summary.em.fit}} helper: Rebuild call
#'
#' @param object an object of class "em.fit"
#'
#' @return a list
rebuild_call <- function(object) {
  
  stopifnot(is.recursive(object))
  stopifnot(!is.null(object$call))
  
  out <- list()
  out$args <- list()
  out$args$old <- as.list(object$call)
  arg_names <- names(args <- formals(as.character(out$args$old[[1]])))
  
  out$args$new <- unlist(c(
    out$args$old[which(names(out$args$old) %in% arg_names)]
    , args[which(!arg_names %in% names(out$args$old))]
  ))[arg_names]
  
  out$call <- sprintf("%s(%s)", as.character(out$args$old[[1]]), paste(arg_names, out$args$new, sep = " = ", collapse = ", "))
  
  out
}

#' \code{\link{summary.em.fit}} helper: Get data info
#'
#' @param object an object of class "em.fit"
#'
#' @return a scalar character vector 
get_data_info <- function(object) {
  with(
    object$info
    , sprintf(
      paste(
        "Annotations data:"
        , "  %d label classes"
        , "  %d items, No. of annotations per item: median = %d, mean (SD) = %1.2f (%1.2f)"
        , "  %d annotators, No. of annotations per annotator: median = %d, mean (SD) = %1.2f (%1.2f)"
        , sep = "\n"
      )
      , labels$n
      , items$n, items$median_annotations, items$mean_annotations, items$sd_annotations
      , annotators$n, annotators$median_annotations, annotators$mean_annotations, annotators$sd_annotations
    )
  )
}

#' \code{\link{summary.em.fit}} helper: Get iterations info
#'
#' @param object an object of class "em.fit"
#'
#' @return a scalar character vector 
get_iter_info <- function(object) {
  iter_info <- with(object
       , list(
         total = length(iterations)
         , start_ll = iterations[[1]]$log_likelihood
         , final_ll = iterations[[length(iterations)]]$log_likelihood
       )
  )
  
  sprintf(
    "Optimization: total %1.0f iterations, final log-likelihood = %+1.3f"
    , iter_info$total
    , iter_info$final_ll
  )
}

#' \code{\link{summary.em.fit}} helper: Get label summary statistics
#'
#' @param lab an numeric vector of class label estimates
#'
#' @importFrom e1071 skewness
#' @return a list
get_lab_sum <- function(lab) {
  list(
    mean = mean(lab)
    , sd = sd(lab)
    , skewness = e1071::skewness(lab)
    , quartiles = quantile(lab, c(.25, .5, .75))
  )
}

#' \code{\link{summary.em.fit}} helper: Printer for label summary statistics
#'
#' @param nm label class name
#' @param vals label class summary statistics
#'
#' @return a scalar character vector 
sprintf_lab_sum <- function(nm, vals) {
  sprintf(
    "%s:   %0.4f   %0.4f  %+0.4f   %0.4f  %0.4f  %0.4f"
    , nm
    , vals$mean
    , vals$sd
    , vals$skewness
    , vals$quartiles[1]
    , vals$quartiles[2]
    , vals$quartiles[3]
  )
}


#' Pad white space
#'
#' @param x character vector
#' @param n total No. characters
#' @param left padding to the left (default: \code{TRUE}). otherwise padding to the right
#'
#' @return a character vector
pad_ws <- Vectorize(function(x, n, left = TRUE) {
  if ((l <- nchar(x)) == n) 
    return(x)
  
  pad <- n-l
  stopifnot(pad >= 0L)
  pad <- paste(rep(" ", times = pad), collapse = "")
  
  if (left)
    return(paste0(pad, x))
  else 
    return(paste0(x, pad))
})


#' \code{\link{summary.em.fit}} helper: build summary table
#'
#' @param x a \code{data.frame} or \code{\link[tibble]{tibble}} with numeric columns, one for each label class 
#' @param indent integer. No. white-space indentations added to the left (default: 4)
#'
#' @return a scalar character vector 
build_sum_tab <- function(x, indent = 4L) {
  
  lab_sums <- lapply(x, get_lab_sum)
  
  nms <- sQuote(names(lab_sums))
  
  pad <- max(max(nchar(nms), 5))
  
  lab_sums_char <- mapply(
    sprintf_lab_sum
    , nm = pad_ws(nms, n = pad)
    , vals = lab_sums
  )
  
  ws_indent <- paste(rep(" ", indent), collapse = "")
  
  out <- paste(
    sprintf(
      "%s%s"
      , paste0(ws_indent, paste(rep(" ", pad+39), collapse = ""))
      , "Quantiles"
    )
    , sprintf(
      "    %sMean     S.D.    Skewness    25%%     50%%     75%%  "
      , pad_ws(x = "Label", n = pad+5L, left = FALSE)
    )
    , paste0(ws_indent, paste(rep("-", pad+54), collapse = ""))
    , paste(ws_indent, lab_sums_char, collapse = "\n")
    , sep = "\n"
  )
  
  out
}

# summary.em.fit ----
#' S3 method for class "summary.em.fit"
#' 
#' @description summary method for class "em.fit"
#'
#' @method summary em.fit
#' @param object an object of class "em.fit", usually a result of \code{\link{em.bin}}
#' @return an object of class "summary.em.fit"
#' @export
summary.em.fit <- function(object) {
  
  if (!inherits(object, "em.fit"))
    return(summary.default(object))
  
  # 'Call' section
  call <- rebuild_call(fit_cf)
  
  # 'Anntotations data' section
  data_info <- get_data_info(object)
  
  # 'Optimization' section
  iter_info <- get_iter_info(object)
  
  # 'Estimates' section
  nms <- sQuote(object$est_class_prevl[[1]])
  est_class_prevl <- sprintf(
    "  Label class prevalences:\n    %s"
    , paste(
      sprintf(
        "%s = %0.4f"
        , pad_ws(nms, max(nchar(nms)))
        , object$est_class_prevl[[2]]
      )
      , collapse = "\n    "
    )
  )
    
  # item estimates 
  i_ests <- object$est_class_probs[,-c(1,ncol(object$est_class_probs))]
  
  # annotator estimates
  idx <- which(object$est_annotator_params[[2]] == object$est_annotator_params[[3]])
  a_ests <- tidyr::spread(object$est_annotator_params[idx, -2], labeled, est_prob)[,-1]
  
  # combine estimtes output
  ests <- paste(
    "Estimates:", est_class_prevl
    , "", "  Label class probabilites:", build_sum_tab(i_ests)
    , "", "  Annotator accuracies:", build_sum_tab(a_ests)
    , sep = "\n"
  )
  
  # complete summary 
  sum <- paste0(
    "Expectation-Maximization Annotations Model Fit\n\n"
    , "Call:\n", call$call, "\n\n", data_info, "\n\n", iter_info, "\n\n", ests
  )
  
  class(sum) <- "summary.em.fit"
  
  return(sum)
}

#' S3 method for class "summary.em.fit"
#'
#' @method print summary.em.fit
#' @param x an object of class "summary.em.fit"
#' @param ... additional arguments are ignored
#' @return invisible
#' @export
print.summary.em.fit <- function(x, ...) {
  cat(x)
  invisible(NULL)
}