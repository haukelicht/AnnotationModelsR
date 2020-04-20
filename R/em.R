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

  stopifnot(length(.ability.prior) == 1L)
  stopifnot(length(.prevalence.prior) == x$K)

  # capture labels, item and worker IDs
  item_labels_ <- sort(unique(x$y))
  item_ids_ <- unique(x$ii)
  worker_ids_ <- unique(x$jj)

  # create array with J K-by-K matrices defining annotators' confusion matrices
  # initialize with off-diagonal error-rates
  theta_hat <- with(x, array((1-.ability.prior)/(x$K-1), c(K, K, J)))
  dimnames(theta_hat) <- list(item_labels_, item_labels_, worker_ids_)
  # set on-diagonal annotator accuracies for first iteration step
  for (j in 1:x$J) diag(theta_hat[,,j]) <- .ability.prior
  # NOTE: theta_hat[k,k',j] gives the probability that annotator j assigns the
  #       label k' to an item whose true category is k

  # create items-by-categories array to store estimated class probabilities of items
  E_z <- with(x, array(1/K, c(I, K)))
  dimnames(E_z) <- list(item_ids_, item_labels_)
  # NOTE: E_z[i,k] gives the probability that item i's true category is k

  return(
    list(
      theta_hat = theta_hat
      , pi_hat = .prevalence.prior
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
    p$E_z <- array(rep(p$pi_hat, each = x$I), dim(p$E_z), dimnames(p$E_z))

    # for each annotation ...
    for (n in 1:x$N) {
      # update estimated category probabilities by taking into account annotator qualities
      p$E_z[x$ii[n], ] <- p$E_z[x$ii[n], ] * p$theta_hat[, x$y[n], x$jj[n]]
    }

    # rescale so that class probabilities sum to one
    p$E_z <- p$E_z / rowSums(p$E_z)

    # --- Maximization Step --- #
    # add beta smoothing on pi_hat
    p$pi_hat <- colSums(p$E_z) + beta.prior
    # ensure that probabilities sum to one
    p$pi_hat <- as.array(p$pi_hat / sum(p$pi_hat))

    # add alpha smoothing for theta_hat
    cnt <- array(alpha.prior, c(x$K, x$K, x$J), dimnames(p$theta_hat))
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

  # normalize
  p$E_z <- p$E_z/rowSums(p$E_z)
  p$pi_hat <- p$pi_hat/sum(p$pi_hat)

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
#' @param data a data frame with rows corresponding to annotator \eqn{j}'s annotation \ifelse{html}{\out{<em>y<sub>ij</sub></em>}}{\eqn{y_{ij}}} of item \eqn{i}.
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
#' @param .prevalence.prior \code{NULL} (default) or a named vector of numbers in interval (0, 1),
#'      specifying prior probabilities on item label classes \ifelse{html}{\out{<b><em>&pi;</em></b>&circ;}}{\eqn{\hat{\boldsymbol{\pi}}}}.
#'      Must have length equal to number of label classes (i.e., number of unique item labels in \code{data}) and sum to 1.
#'      Defaults to \code{NULL}, in which case prior class label probabilities are set all equal to \eqn{1/K} (where \eqn{K} is the number of item label classes)
#' @param .ability.prior number in interval (0, 1), specifying the annotator ability prior used to initialize ability parameters for first iteration
#'      Defaults to 0.7
#' @param .beta.prior number in interval (0, 1), specifying the beta prior used to perform smoothing on \ifelse{html}{\out{<b><em>&pi;</em></b>&circ;}}{\eqn{\hat{\boldsymbol{\pi}}}} (label class prevalence estimates)
#'      Defaults to 0.01
#' @param .alpha.prior number in interval (0, 1), specifying the alpha prior used to perform smoothing on items' label class probability estimates
#'      Defaults to 0.01
#' @param .min.relative.diff number in interval (0, 1), specifying the minimum relative difference between log-likelihoods of two subsequent iteration steps.
#'      Governs when convergence is reached (the lower, the more iterations).
#'      Defaults to 1e-8
#' @param verbose logical. print-out data description and iteration history?
#'
#' @return An \code{em.fit} object, which is a list with elements
#'     \describe{
#'       \item{call}{the \code{call} object associated with the model}
#'       \item{info}{A list object with elements
#'         \describe{
#'           \item{items}{
#'              A list with elements
#'              \code{n} (number of items),
#'              \code{mean_annotations} (average number of annotations per item),
#'              \code{sd_annotations} (std. dev. of number of annotations per item), and
#'              \code{median_annotations} (median number of annotations per item)
#'           }
#'           \item{annotators}{
#'             A list with elements
#'             \code{n} (number of annotators),
#'             \code{mean_annotations} (average number of annotations per annotator),
#'             \code{sd_annotations} (std. dev. of number of annotations per annotator), and
#'             \code{median_annotations} (median number of annotations per annotator)
#'           }
#'           \item{labels}{
#'             A list with elements
#'             \code{n} (number of label categories/classes), and
#'             \code{props} (proportions of items with label classes)
#'           }
#'         }
#'       }
#'       \item{est_class_probs}{
#'         A \code{\link[tibble]{tibble}} with one row for each item in \code{data} (as indexed by coloumn \code{label.col}).
#'
#'         The first column is named like \code{item.col} and records item identifiers.
#'         The second to \eqn{1+K}th columns record items' label class probability estimates \ifelse{html}{\out{<em>&ycirc;<sub>ij</sub></em>}}{\eqn{\hat{y}_{ij}}} (columns are named like label classes).
#'         The last column \code{majority_vote} recods majority-winner labels in item-level annotations (w/ random tie-breaing where necessary).
#'       }
#'       \item{est_class_prevl}{
#'         A \code{\link[tibble]{tibble}} with as many rows as label classes and columns
#'         \itemize{
#'           \item{\code{<label.col>}: the original label identifier column name}
#'           \item{\code{est_prob}: the label class prevalence estimate}
#'           \item{\code{prop_labels}: the proportion of items with this row's label among model-induced labelings}
#'           \item{\code{prop_mv}: the proportion of items with this row's label among majority-winner labelings}
#'         }
#'       }
#'       \item{est_annotator_params}{
#'         A \code{\link[tibble]{tibble}} with \ifelse{html}{\out{J&times;K&times;K}}{\eqn{J \times K \times K}} rows recording all \eqn{J} coders' \ifelse{html}{\out{K&times;K}}{\eqn{K \times K}} ability and error-rate estimate matrices.
#'         Columns:
#'         \itemize{
#'           \item{\code{<annotator.col>}: the original annotator identifier column name},
#'           \item{\code{<label.col>}: the "true" label class (column named like the original label identifier column)}
#'           \item{\code{labeled}: the label assigned by annotator \eqn{j}}
#'           \item{
#'             \code{est_prob}:
#'             the probability annotator \eqn{j} assigns \ifelse{html}{\out{<em>y<sub>ij</sub></em>}}{\eqn{y_{ij}}}
#'              to item \eqn{i} with true label \ifelse{html}{\out{<em>&ytilde;<sub>i</sub></em>}}{\eqn{\tilde{y}_i}}
#'           }
#'         }
#'       }
#'       \item{iter_log}{
#'         A list with as many elements as iterations until convergence.
#'         Each element is a list with elements
#'         \itemize{
#'           \item{\code{idx}: the iteration counter, starts at 1}
#'           \item{\code{log_likelihood}: the log-likelihood value of the current iteration}
#'           \item{\code{relative_diff}: the difference in log-likelihood values of the current relative to the previous iteration}
#'         }
#'       }
#'     }
#' @import rlang
#' @import dplyr
#' @importFrom tidyr pivot_wider
#' @importFrom stats runif
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # inspect the Dawid-Skene example dataset
#' head(dawidskene)
#' # fit a per-annotator model using the EM algorithm
#' fit <- em(dawidskene, item.col = patient, annotator.col = observer, label.col = diagnosis)
#' # view summary
#' summary(fit)
#' }
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

  # check input aguments

  # `data`
  if(!inherits(data, "data.frame"))
    stop("`data` must be a ", sQuote("data.frame"), " object or inherit from the ", sQuote("data.frame"), " class.")

  # `item.col`, `annotator.col`, and `label.col`
  args <- as.list(match.call())
  req_cols <- c("item.col", "annotator.col", "label.col")
  req_colns <- c(args$item.col, args$annotator.col, args$label.col)

  tmp <- rep(TRUE, 3)
  if (
    is.null(req_colns)
    || any(tmp <- !req_cols %in% names(args))
    || any(tmp <- vapply(req_colns, is.null, NA))
  ) {
    err_msg <- character()
    for (c in which(tmp))
      err_msg[length(err_msg)+1L] <- paste(req_cols[c], "cannot be NULL.")

    stop(err_msg)
  }

  if (!all(req_colns %in% colnames(data)))
    stop("`data` must be a ", sQuote("data.frame"), " object with columns: ", paste0(sQuote(req_colns), collapse = ", "))

  if (any((nm <- c("_item", "_annotator", "_label")) %in% colnames(data)))
    stop("`data` cannot have reserved columns ", paste0(sQuote(nm), collapse = ", "))

  # `.prevalence.prior`
  if (!is.null(.prevalence.prior)) {
    qlabs <- sQuote(sort(unique(data[[args$label.col]])))
    err_msg <- sprintf("`.prevalence.prior` must be NULL (default) or a named, non-negative numeric vector with names in [%s] and values summing to 1. See `?em`.", paste(qlabs, collapse = ", "))
    if (
      !is.double(.prevalence.prior)
      || any(.prevalence.prior <= 0)
      || sum(.prevalence.prior) != 1
      || is.null(names(.prevalence.prior))
    ) {
      stop(err_msg)
    }

    if (!all(w_ <- unique(data[[args$label.col]]) %in% names(.prevalence.prior))) {
      missing_labs <- paste(qlabs[!w_], collapse = ", ")
      err_msg <- sprintf(
        "Element missing from `.prevalence.prior`: No value%s provided for label class%s %s. See `?em`."
        , ifelse(sum(!w_) > 1, "s", "")
        , ifelse(sum(!w_) > 1, "es", "")
        , ifelse(sum(!w_) > 1, sprintf("[%s]", missing_labs), missing_labs)
      )
      stop(err_msg)
    }

    if (any(w_ <- !names(.prevalence.prior) %in% unique(data[[args$label.col]]))) {
      stop(paste("`.prevalence.prior` has too many elements:", err_msg))
    }
  }

  if (length(labs <- unique(data[[args$label.col]])) == 1) {
    warning(
      "There is only one label class in your data: ", sQuote(labs)
      , ". Without variation in annotations, EM estimates will be non-sensical."
      , call. = FALSE
    )
  }

  # Setup data:

  # assign unique item, annotator and label indexes
  dm <- data %>%
    ungroup() %>%
    mutate(
      `_item` = group_indices(., !!enquo(item.col)) ,
      `_annotator` = group_indices(., !!enquo(annotator.col)),
      `_label` = group_indices(., !!enquo(label.col))
    )

  # get a mapping of label.col values to indexes
  label_map <- as.list(unique(dm[, c(as.character(args$label.col), "_label")]))
  label_map <- setNames(label_map[[2]], label_map[[1]])

  # create item, annotator and label indexes
  dat <- create_em_indixes(dm, echo)

  # if no prevalence prior is provided, set equal prior item label class probabilities
  if (is.null(.prevalence.prior)){
    .prevalence.prior <- setNames(with(dat, array(1/K, K)), unique(dat$y))
  } else {
    names(.prevalence.prior) <- names(label_map[names(.prevalence.prior)])
    .prevalence.prior <- as.array(.prevalence.prior[unique(dat$y)])
  }

  # Estimation:

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
    pivot_wider(names_from = as.character(args$label.col), values_from = "est_prob") %>%
    left_join(majority_votes, by = as.character(args$item.col))


  # label class prevalence
  pi_out <- tibble(
    cat_label = sort(unique(dat$y)),
    est_prob = p$pi_hat,
    prop_labels = as.vector(prop.table(table(dat$y))[names(p$pi_hat)])
  ) %>%
    left_join(
      unique(select(dm, !!enquo(label.col), `_label`))
      , by = c("cat_label" = "_label")
    ) %>%
    select(!!enquo(label.col), est_prob, prop_labels) %>%
    left_join(
      count(z_out, majority_vote) %>%
        mutate(n = n/sum(n)) %>%
        rename(!!enquo(label.col) := majority_vote, prop_mv = n)
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



