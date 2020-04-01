# Estimator helpers ----

#' Get annotations info
#'
#' @param x a \code{data.frame} with columns '_item', '_annotator' and '_label'
#'     holding corresponding integer indexes

#' @return a list
#'
#' @importFrom stats median sd setNames
get_annotations_info <- function(x) {

  info <- list()

  temp <- table(x[["_item"]])
  info$items <- list()
  info$items$n <- length(temp)
  info$items$mean_annotations <- mean(temp)
  info$items$sd_annotations <- tryCatch(sd(temp), error = function(err) 0.0)
  info$items$median_annotations <- median(temp)

  temp <- table(x[["_annotator"]])
  info$annotators <- list()
  info$annotators$n <- length(temp)
  info$annotators$mean_annotations <- mean(temp)
  info$annotators$sd_annotations <- tryCatch(sd(temp), error = function(err) 0.0)
  info$annotators$median_annotations <- median(temp)

  temp <- table(x[["_label"]])
  info$labels <- list()
  info$labels$n <- length(temp)
  info$labels$props <- setNames(as.vector(prop.table(temp)), names(temp))

  return(info)
}

#' Create EM indexes for binar model
#'
#' @param x a \code{data.frame} with columns '_item', '_annotator' and '_label'
#'     holding corresponding integer indexes
#' @param .echo a function used to print messages
#'
#' @return a list with elements
#'    \itemize{
#'      \item 'ii': item index values
#'      \item 'I': number of items
#'      \item 'jj': annotator index values
#'      \item 'J': number of annotators
#'      \item 'y': label index values
#'      \item 'K': number of label classes
#'      \item 'N': number of annotations
#'    }
create_em_indixes <- function(x, .echo) {

  out <- list()

  # item vector and number of items
  out$ii <- x[["_item"]]
  out$I <- length(unique(out$ii))

  # annotator vector and number of annotators
  out$jj <- x[["_annotator"]]
  out$J <- length(unique(out$jj))

  # label vector and number of categories (unique labels)
  out$y <- x[["_label"]]
  out$K <- length(unique(out$y))

  # number of annotations
  out$N <- length(out$ii)

  # test if complete panel design (i.e., all annotators labeled all items)
  complete_panel_design <- TRUE
  for (i in 1:out$I) {
    if (nrow(x[x[["_item"]] == i,]) != out$J)
      complete_panel_design <- FALSE
  }

  # report data sizes (if verbose)
  .echo(paste(
    "Data input:"
    , out$I, "items,"
    , out$J, "annotators, and"
    , out$K, "label categories"
  ))
  .echo(sprintf("Note: annotations made in %scomplete panel design", ifelse(complete_panel_design, "a ", "an in")), "\n")
  .echo("\n")

  return(out)
}

#' Initialize EM-parameters for binar model
#'
#' @param x a list as returned by \code{\link{create_em_indixes}}
#'
#' @param .ability.prior number in interval (0, 1), specifying the annotator ability prior used to initialize ability parameters for first iteration
#'      Defaults to 0.7
#'
#' @param .prevalence.prior vector of numbers in interval (0, 1),
#'      specifying prior probabilities on item label classes \eqn{\hat{\boldsymbol{pi}}}.
#'      Most have length equal to number of label classes (i.e., number of unique item labels in \code{data}).
#'      Defaults to \code{NULL}, in which case prior class label probabilities are set all equal to \eqn{1/K} (where \eqn{K} is the number of item label classes)
#'
#' @return a list with elements
#'     \itemize{
#'       \item 'theta_hat': annotator parameter matrices (a 3d-tensor with extension \eqn{K* K* J})
#'       \item 'pi_hat': prevalence estimates (a \eqn{K}-length numeric vector)
#'       \item 'E_z': \eqn{I* K} matrix of item label class probability estimates
#'     }
init_em_params <- function(x, .ability.prior, .prevalence.prior){

  # check input argument values
  stopifnot(is.list(x))
  stopifnot(all(c("y", "ii", "jj", "K", "J") %in% names(x)))

  stopifnot(length(.ability.prior) == length(x$K))
  stopifnot(length(.prevalence.prior) == length(x$K))

  # capture labels, item and worker IDs
  item_labels_ <- names(rev(sort(table(x$y))))
  item_ids_ <- unique(x$ii)
  worker_ids_ <- unique(x$jj)

  # create array with J K-by-K matrices defining annotators' confusion matrices
  theta_hat <- with(x, array(NA, c(K, K, J)))
  dimnames(theta_hat) <- list(item_labels_, item_labels_, worker_ids_)

  # set annotator accuracies for first iteration step
  for (j in 1:x$J) {
    theta_hat[,,j] <- (1-.ability.prior)/(x$K-1) # set off-diagonal (error-rates)
    diag(theta_hat[,,j]) <- .ability.prior # and diagonal
  }
  # NOTE: theta_hat[k,k',j] gives the probability that annotator j assigns the
  #       label k' to an item whose true category is k

  # set initial class prevalences/probabilites
  # pi_hat <- with(x, array(1/K, K))
  pi_hat <- .prevalence.prior
  dimnames(pi_hat) <- list(item_labels_)

  # create items-by-categories array to store estimated class probabilities of items
  E_z <- with(x, array(1/K, c(I, K)))
  dimnames(E_z) <- list(item_ids_, item_labels_)
  # NOTE: E_z[i,k] gives the probability that item i's true category is k

  return(
    list(
      theta_hat = theta_hat
      , pi_hat = pi_hat
      , E_z = E_z
    )
  )
}


#' Expectation-Maximization iterator for binary model
#'
#' @param p a list of initialized parameters as returned by \code{\link{init_em_params}}
#' @param max.iters positive integer, specifying the maximum number of EM iterations to perform.
#'      (CAUTION: low values may cause stopping before convergence)
#'      Defaults to 100.
#' @param beta.prior number in interval (0, 1), specifying the beta prior used to perform smoothing on \eqn{\hat{pi}} (label class prevalence estimates)
#' @param alpha.prior number in interval (0, 1), specifying the alpha prior used to perform smoothing on items' label class probability estimates
#' @param min.relative.diff number in interval (0, 1), specifying the minimum relative difference between log-likelihoods of two subsequent iteration steps.
#'      Governs when convergence is reached (the lower, the more iterations).
#'      Defaults to 1e-8
#' @describeIn create_em_indixes
#'
#' @return a list with elements
#'     \itemize{
#'       \item 'theta_hat': annotator parameter matrices (a 3d-tensor with extension \eqn{K* K* J})
#'       \item 'pi_hat': prevalence estimates (a \eqn{K}-length numeric vector)
#'       \item 'E_z': \eqn{I* K} matrix of item label class probability estimates
#'       \item 'iter_log': a list with as many elements as iterations until convergence.
#'           Each element is a list with elements
#'           'idx' (the iteration counter, starts at 1),
#'           'log_likelihood' (the log-likelihood value of the current iteration),
#'           'relative_diff' (the difference in log-likelihood values of the current relative to the previous iteration)
#'     }
em_iterate <- function(x, p, max.iters, beta.prior, alpha.prior, min.relative.diff, .echo) {

  # setup for first iteration step
  last_log_posterior = -Inf

  p$iter_log <- list()
  for (iter_step in 1:max.iters) {

    # --- Expectation Step --- #
    # set estimated class probabilities to current estimates of category prevalences
    for (i in 1:x$I)
      p$E_z[i,] <- p$pi_hat

    # for each annotation ...
    for (n in 1:x$N) {
      # update estimated category probabilities by taking into account annotator qualities
      p$E_z[x$ii[n], ] <- p$E_z[x$ii[n], ] * p$theta_hat[, x$y[n], x$jj[n]]
    }

    # rescale so that class probabilities sum to one
    p$E_z <- p$E_z / rowSums(p$E_z)

    # --- Maximization Step --- #
    # add beta smoothing on pi_hat
    nms <- dimnames(p$pi_hat)
    p$pi_hat <- as.array(rep(beta.prior, x$K))
    dimnames(p$pi_hat) <- nms
    for (i in 1:x$I)
      p$pi_hat <- p$pi_hat + p$E_z[i,]

    # ensure that probabilities sum to one
    p$pi_hat <- p$pi_hat / sum(p$pi_hat)

    # add alpha smoothing for theta_hat
    cnt <- array(alpha.prior, c(x$K, x$K, x$J))
    dimnames(cnt) <- dimnames(p$theta_hat)
    # for each annotation ...
    for (n in 1:x$N) {
      # ... add alpha prior to estimated class probabilities
      cnt[, x$y[n], x$jj[n]] <- cnt[, x$y[n], x$jj[n]] + p$E_z[x$ii[n], ]
    }

    # for each annotator ...
    for (j in 1:x$J) {
      # ... and for each category ...
      for (k in 1:x$K)
        # ensure that individual error-rates sum to one
        p$theta_hat[k,,j] <- cnt[k,,j] / sum(cnt[k,,j])
    }

    # create array for individual class posterior probabilities
    pp <- matrix(rep(p$pi_hat, x$I), ncol = length(p$pi_hat), byrow = TRUE)
    dimnames(pp) <- dimnames(p$E_z)

    # for each annotation ...
    for (n in 1:x$N) {
      # ... and for each category ...
      for (k in 1:x$K) {
        # ... update
        pp[x$ii[n], k] <- pp[x$ii[n], k] * p$theta_hat[k, x$y[n], x$jj[n]]
      }
    }

    # compute log of posterior probability of the data given current parameter estimates
    log_posterior <- sum(log(rowSums(pp)))

    # decide wether to progress
    if (iter_step == 1) {
      p$iter_log[[iter_step]] <- list(step = iter_step, log_likelihood = log_posterior, relative_diff = NA_real_)
      .echo("Estimation results:")
      .echo(
        sprintf(
          "  iteration nr. %3.0f: log posterior = %.7f"
          , iter_step, round(log_posterior, 7)
        )
      )
    } else {
      diff <- log_posterior - last_log_posterior
      relative_diff <- abs(diff / last_log_posterior)
      p$iter_log[[iter_step]] <- list(step = iter_step, log_likelihood = log_posterior, relative_diff = relative_diff)
      .echo(
        sprintf(
          "\r  iteration nr. %3.0f: log posterior = %.7f, relative difference = %.9f"
          , iter_step, round(log_posterior, 7), round(relative_diff, 9)
        )
        , appendLF = FALSE
      )
      if (relative_diff < min.relative.diff) {
        .echo("\nConverged!")
        break
      }
    }
    last_log_posterior <- log_posterior
  }

  return(p)
}


# Estimator ----

#' Estimate Dawid-Skene Expectation Maximization model for categorical annotation data
#'
#' @description Function fits an Dawid-Skene annotation model to categorical annotation data
#'     using the Expectation Maximization algorithm.
#'
#' @note Code adapted from Bob Carpenter's implementation
#'    of Dawid and Skene (1979) "Maximum Likelihood Estimation of Observer Error-Rates Using the EM Algorithm"
#'    (see \url{https://github.com/bob-carpenter/anno/blob/master/R/em-dawid-skene.R})
#'
#' @param data a data frame with rows corresponding to annotator \eqn{j}'s annotation \eqn{y_{ij}} of item \eqn{i}.
#'
#' @param item.col (unquoted name of) column containing item identifiers
#'
#' @param annotator.col (unquoted name of) column containing annotator identifiers
#'
#' @param label.col (unquoted name of) column containing label identifiers
#'
#' @param max.iters positive integer, specifying the maximum number of EM iterations to perform.
#'      (CAUTION: low values may cause stopping before convergence)
#'      Defaults to 100.
#' @param .prevalence.prior vector of numbers in interval (0, 1),
#'      specifying prior probabilities on item label classes \eqn{\hat{\boldsymbol{pi}}}.
#'      Most have length equal to number of label classes (i.e., number of unique item labels in \code{data}).
#'      Defaults to \code{NULL}, in which case prior class label probabilities are set all equal to \eqn{1/K} (where \eqn{K} is the number of item label classes)
#' @param .ability.prior number in interval (0, 1), specifying the annotator ability prior used to initialize ability parameters for first iteration
#'      Defaults to 0.7
#' @param .beta.prior number in interval (0, 1), specifying the beta prior used to perform smoothing on \eqn{\hat{\boldsymbol{pi}}} (label class prevalence estimates)
#'      Defaults to 0.01
#' @param .alpha.prior number in interval (0, 1), specifying the alpha prior used to perform smoothing on items' label class probability estimates
#'      Defaults to 0.01
#' @param .min.relative.diff number in interval (0, 1), specifying the minimum relative difference between log-likelihoods of two subsequent iteration steps.
#'      Governs when convergence is reached (the lower, the more iterations).
#'      Defaults to 1e-8
#' @param verbose logical. print-out data description and iteration history?
#'
#' @return a \code{em.fit} object, which is a list with elements
#'     \enumerate{
#'       \item{'call': the \code{call} object associated with the model
#'       \item{'info': a list object with elements
#'         \enumerate{
#'           \item 'items': a list with elements
#'             'n' (number of items),
#'             'mean_annotations' (average number of annotations per item),
#'             'sd_annotations' (std. dev. of number of annotations per item), and
#'             'median_annotations' (median number of annotations per item)
#'           \item 'annotators': a list with elements
#'             'n' (number of annotators),
#'             'mean_annotations' (average number of annotations per annotator),
#'             'sd_annotations' (std. dev. of number of annotations per annotator), and
#'             'median_annotations' (median number of annotations per annotator)
#'           \item 'labels': a list with elements
#'             'n' (number of label categories/classes), and
#'             'props' (proportions of items with label classes)
#'         }
#'       }
#'       \item{'est_class_probs': a \code{\link[tibble]{tibble}} with columns
#'         \itemize{
#'           \item '<item.col>' (original item identifier column name),
#'           \item one for each label class occuring in \code{data$<label.col>}
#'               with values holding label class posterior probability estimates \eqn{\hat{z}_{ij}},
#'           \item[] and
#'           \item 'majority_vote' (majority-winner label w/ random tie-breaing where required)
#'         }
#'       }
#'       \item{'est_class_prevl': a \code{\link[tibble]{tibble}} with as many rows as label classes and columns
#'         \enumerate{
#'           \item '<label.col>' (the original label identifier column name),
#'           \item 'est_prob' (the posterior prevalence estimate),
#'           \item 'prop_labels' (the label proportions, where labels are induced by applying a .5 threshold to posterior label class probability estimates),
#'           \item[] and
#'           \item 'prop_mv' (the majority-winner label proportions)
#'         }
#'       }
#'       \item{'est_annotator_params': a \code{\link[tibble]{tibble}} with \eqn{J* K* K} rows and columns
#'         \enumerate{
#'           \item '<annotator.col>' (the original annotator identifier column name),
#'           \item 'label' (the hypothetical 'true' label \eqn{\tilde{y}_i}),
#'           \item 'labeled' (\eqn{y_{ij}}, the label annotator \eqn{j} assigns to item \eqn{i}),
#'           \item[] and
#'           \item 'est_prob' (the probability annotator \eqn{j} assigns \eqn{y_{ij}} to item \eqn{i} with true label \eqn{\tilde{y}_i})
#'         }
#'       }
#'       \item{ 'iter_log': a list with as many elements as iterations until convergence.
#'         Each element is a list with elements
#'         \enumerate{
#'           \item 'idx' (the iteration counter, starts at 1),
#'           \item 'log_likelihood' (the log-likelihood value of the current iteration),
#'           \item[] and
#'           \item 'relative_diff' (the difference in log-likelihood values of the current relative to the previous iteration)
#'         }
#'       }
#'     }
#' @import rlang
#' @import dplyr
#' @importFrom tidyr spread
#' @importFrom tibble enframe
#' @importFrom stats runif
#'
#' @export
em <- function(
  data,
  item.col,
  annotator.col,
  label.col,
  max.iters = 100,
  .ability.prior = 0.7,
  .prevalence.prior = NULL,
  .beta.prior = 0.01,
  .alpha.prior = 0.01,
  .min.relative.diff = 1e-8,
  verbose = TRUE
){


  # define print function
  if (verbose)
    echo <- match.fun("message")
  else
    echo <- function(...) {}

  # test input
  stopifnot(is.data.frame(data))

  args <- as.list(match.call())

  req_colns <- c(args$item.col, args$annotator.col, args$label.col)

  if (!all(req_colns %in% colnames(data)))
    stop("input to argument `data' must be a data.frame object with columns: ", paste0("'", req_colns, "'", collapse = ", "))

  if (any((nm <- c("_item", "_annotator", "_label")) %in% colnames(data)))
    stop("`data' cannot have reserved columns ", paste0("'", nm, "'", collapse = ", "))

  # Setup data:

  # assign unique item, annotator and label indixes
  dm <- data %>%
    ungroup() %>%
    mutate(
      `_item` = group_indices(., !!enquo(item.col)) ,
      `_annotator` = group_indices(., !!enquo(annotator.col)),
      `_label` = group_indices(., !!enquo(label.col))
    )

  # create item, annotator and label indexes
  dat <- create_em_indixes(dm, echo)

  # Estimation:

  # if no prevalence prior is provided, set equal prior item label class probabilities
  if (is.null(.prevalence.prior))
      .prevalence.prior <- with(dat, array(1/K, K))

  # initialize parameters
  prms <- init_em_params(dat, .ability.prior, .prevalence.prior)

  # erform expectation-maximization iterations
  p <- em_iterate(
    x = dat
    , p = prms
    , max.iters = max.iters
    , beta.prior = .beta.prior
    , alpha.prior = .alpha.prior
    , min.relative.diff = .min.relative.diff
    , .echo = echo
  )

  # Prep output:

  # compute proportion of labels
  prop_labels <- rep(0, dat$K)
  for (k in 1:dat$K)
    prop_labels[k] <- sum(dat$y == k) / dat$N

  # compute majority votes
  majority_votes <- dm %>%
    group_by(`_item`, !!enquo(item.col), !!enquo(label.col)) %>%
    summarize(votes = n_distinct(!!enquo(annotator.col))) %>%
    group_by(`_item`, !!enquo(item.col), !!enquo(label.col), votes) %>%
    mutate(
      tie_breaker = runif(1),
      majority_vote = !!enquo(label.col)
    ) %>%
    group_by(!!enquo(item.col)) %>%
    filter(votes == max(votes)) %>%
    slice(which.max(tie_breaker)) %>%
    select(!!enquo(item.col), majority_vote) %>%
    ungroup()

  # item-specific label class probabilities
  z_out <- tibble(
    this_item = rep(1:dat$I, each=dat$K),
    cat_label = rep(1:dat$K, times=dat$I),
    est_prob = as.vector(t(p$E_z))
  ) %>%
    left_join(
      unique(select(dm, !!enquo(item.col), `_item`))
      , by = c("this_item" = "_item")
    ) %>%
    left_join(
      unique(select(dm, !!enquo(label.col), `_label`))
      , by = c("cat_label" = "_label")
    ) %>%
    select(!!enquo(item.col), !!enquo(label.col), est_prob) %>%
    spread(2, est_prob) %>%
    left_join(majority_votes, by = as.character(args$item.col))

  # label class prevalence
  pi_out <- tibble(
    cat_label = sort(unique(dat$y)),
    est_prob = p$pi_hat,
    prop_labels = prop_labels
  ) %>%
    left_join(
      unique(select(dm, !!enquo(label.col), `_label`))
      , by = c("cat_label" = "_label")
    ) %>%
    select(!!enquo(label.col), est_prob, prop_labels) %>%
    left_join(
      prop.table(table(z_out$majority_vote)) %>%
        enframe() %>%
        rename(!!enquo(label.col) := name, prop_mv = value)
      , by = as.character(args$label.col)
    )

  # annotator ability estimats
  theta_out <- tibble(
    this_annotator = with(dat, rep(1:J, each=K*K)),
    cat_label = with(dat, rep(rep(1:K, times=K), times=J)),
    an_labeled = with(dat, rep(rep(1:K, each=K), times=J)),
    est_prob = as.vector(p$theta_hat)
  ) %>%
    left_join(
      unique(select(dm, !!enquo(annotator.col), `_annotator`))
      , by = c("this_annotator" = "_annotator")
    ) %>%
    left_join(
      unique(select(dm, !!enquo(label.col), `_label`))
      , by = c("an_labeled" = "_label")
    ) %>%
    rename(labeled := !!enquo(label.col)) %>%
    left_join(
      unique(select(dm, !!enquo(label.col), `_label`))
      , by = c("cat_label" = "_label")
    ) %>%
    select(!!enquo(annotator.col), !!enquo(label.col), labeled, est_prob)

  # return output as list
  structure(
    list(
      # the call
      call = match.call(),
      # info
      info = get_annotations_info(dm),
      # item label class probability estimates
      est_class_probs = z_out,
      # class prevalence estimates
      est_class_prevl = pi_out,
      # annotator ablity parameter estimates (accuracies and erro rates)
      est_annotator_params = theta_out,
      # iteration history
      iterations = p$iter_log
    )
    , class = c("em.fit", "list")
  )

}

# EM-fit Plots ----

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



