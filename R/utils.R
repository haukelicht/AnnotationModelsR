#' Get labeling
#'
#' @description Given a model object, function induces a labeling from class category probability estimates.
#'
#' @param object the model object
#'
#' @return a \code{\link[tibble]{tibble}} mapping items to class category estimates
#'
#' @export
get_labeling <- function(object) {

  if (!inherits(object, "em.fit")) {
    obj_class <- class(object)
    if (length(obj_class) > 1L)
      obj_class <- paste0("[", paste0('"', obj_class, '"', collapse = ", "), "]")
    stop("`get_labeling` not implemented for object of class ", obj_class)
  }

  k <- ncol(object$est_class_probs)
  cats <- names(object$est_class_probs)[-c(1, k)]
  object$est_class_probs$labeling <- cats[apply(object$est_class_probs[-c(1, k)], 1, which.max)]

  return(object$est_class_probs)
}
