
#' EM-fit plot: Item class probability estimate distributions
#'
#' @param em.fit an \code{em.fit} object es returned by \code{em}
#'
#' @return a \code{em.fit.plot} object inheriting from \code{\link[ggplot2]{ggplot}}
#'
#' @import ggplot2
#'
#' @export
emplot_item_class_distributions <- function(em.fit) {

  stopifnot(inherits(em.fit, "em.fit"))

  temp <- em.fit$est_class_probs
  temp[[paste0("_label_", names(temp)[2], "_")]]  <- temp[[2]] > temp[[3]]
  temp[["_prob_"]] <- ifelse(temp[[5]], temp[[2]], temp[[3]])
  temp[[5]] <- ifelse(temp[[5]], names(temp)[2], names(temp)[3])
  p <- ggplot(temp, aes(x = `_prob_`)) +
    geom_density(color = NA, fill = "black", alpha = .5) +
    facet_grid(cols = vars(`_label_no_`)) +
    labs(
      title = "Distribution of posterior label class probability estimates  by posterior labels"
      , subtitle = sprintf(
        paste(
          "Left panel: all items labeled '%s',"
          , "right panel all items labeled '%s'"
          , "based on .5 threshold"
        )
        , names(temp)[2], names(temp)[3]
      )
      , x = "Posterior class probability"
      , y = "Density"
    )

  p <- structure(p, class = c("em.fit.plot", class(p)))

  invisible(p)
}
#' EM-fit plot: Annotator ability estimate distributions
#'
#' @param em.fit an \code{em.fit} object es returned by \code{em}
#'
#' @return a \code{em.fit.plot} object inheriting from \code{\link[ggplot2]{ggplot}}
#'
#' @import ggplot2
#' @importFrom stats as.formula
#'
#' @export
emplot_annotator_ability_distributions <- function(em.fit) {
  stopifnot(inherits(em.fit, "em.fit"))

  temp <- unique(em.fit$est_annotator_params[[3]])

  p <- ggplot(
    em.fit$est_annotator_params[em.fit$est_annotator_params[[3]] == temp[1], ]
    , aes(est_prob)
  ) +
    geom_histogram(bins = nrow(em.fit$est_annotator_params)/20, alpha = .5) +
    facet_grid(as.formula(paste("~", names(em.fit$est_annotator_params)[2]))) +
    labs(
      title = "Distribution of posterior annotator ability parameter estimates by item class label"
      , subtitle = "Left and right panels depict label-specific estimates (inverses of error-rates)"
      , x = "Accuracy"
      , y = "Count"
    )

  p <- structure(p, class = c("em.fit.plot", class(p)))

  invisible(p)
}

#' EM-fit plot: estimation convergence
#'
#' @param em.fit an \code{em.fit} object es returned by \code{em}
#'
#' @return a \code{em.fit.plot} object inheriting from \code{\link[ggplot2]{ggplot}}
#'
#' @import ggplot2
#'
#' @export
emplot_convergence <- function(em.fit) {
  stopifnot(inherits(em.fit, "em.fit"))

  p <- ggplot(
    data.frame(
      ll = unlist(lapply(em.fit$iterations, `[[`, "log_likelihood"))
      , iter = 1:length(em.fit$iterations)
    )
    , aes(x = iter, y = ll)
  ) +
    geom_line() +
    labs(
      title = "Expectation-Maximization log-likelihood trace plot"
      , x = "Iteration"
      , y = "log-Likelihood"
    )

  p <- structure(p, class = c("em.fit.plot", class(p)))

  invisible(p)
}
