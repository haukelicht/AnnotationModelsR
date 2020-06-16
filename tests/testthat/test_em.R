context("em()")
library(tibble)

# check em() arguments ----
test_that("em() checks arguments.", {
  # test that only data.frames are accepted
  expect_error(
    em(
      data = matrix()
      , item.col = NULL
      , annotator.col = NULL
      , label.col = NULL
    )
    , paste("`data` must be a", sQuote("data.frame"), "object")
  )

  # test that *.col are columns of `data`
  expect_error(
    em(
      data = data.frame()
      , item.col = NULL
      , annotator.col = NULL
      , label.col = NULL
    )
    , "cannot be NULL"
  )

  # test that reserved column names are protected
  expect_error(
    em(
      data = tibble("_item" = NA, "_annotator" = NA, "_label" = NA)
      , item.col = "_item"
      , annotator.col = "_annotator"
      , label.col = "_label"
    )
    , "`data` cannot have reserved columns"
  )

  # test that .prevalence.prior is checked
  expect_error(
    em(
      data = dawidskene
      , item.col = patient
      , annotator.col = observer
      , label.col = diagnosis
      # when not real-valued (double)
      , .prevalence.prior = 1:4
    )
    , "must be NULL \\(default\\) or a named, non-negative numeric vector"
  )
  expect_error(
    em(
      data = dawidskene
      , item.col = patient
      , annotator.col = observer
      , label.col = diagnosis
      # when unnamed
      , .prevalence.prior = c(.2, .2, .4, .2)
    )
    , "must be NULL \\(default\\) or a named, non-negative numeric vector"
  )
  expect_error(
    em(
      data = dawidskene
      , item.col = patient
      , annotator.col = observer
      , label.col = diagnosis
      # when not all positive
      , .prevalence.prior = c(-.2, .2, .4, .2)
    )
    , "must be NULL \\(default\\) or a named, non-negative numeric vector"
  )
  expect_error(
    em(
      data = dawidskene
      , item.col = patient
      , annotator.col = observer
      , label.col = diagnosis
      # when not summing to 1
      , .prevalence.prior = c(.2, .2, .4, .3)
    )
    , "must be NULL \\(default\\) or a named, non-negative numeric vector"
  )
  expect_error(
    em(
      data = dawidskene
      , item.col = patient
      , annotator.col = observer
      , label.col = diagnosis
      # when missing label class
      , .prevalence.prior = c("1" = .3, "2" = .3, "3" = .4)
    )
    , "Element missing from `\\.prevalence\\.prior`:"
  )
  expect_error(
    em(
      data = dawidskene
      , item.col = patient
      , annotator.col = observer
      , label.col = diagnosis
      # when superflous label class
      , .prevalence.prior = c("1" = .2, "2" = .2, "3" = .2, "4" = .2, "5" = .2)
    )
    , "`\\.prevalence\\.prior` has too many elements:"
  )
})

# test create_em_indixes() ----
test_that("create_em_indixes()", {
  td <- tibble(
    "_item" = rep(1:2, each = 2)
    , "_annotator" = rep(1:2, times = 2)
    , "_label" = rep(1:2, each = 2)
  )

  res <- create_em_indixes(td, .echo = function(...) NULL)

  # test that returns list
  expect_is(res, "list")

  # test output
  expect_identical(names(res), c("ii", "I", "jj", "J", "y", "K", "N"))
  expect_true(all(vapply(res, is.integer, NA)))
  expect_equal(length(res$y), 4L)
  expect_equal(res$I, 2L)
  expect_equal(res$J, 2L)
  expect_equal(res$N, 4L)
})

# test init_em_params() ----
test_that("init_em_params()", {
  td <- tibble(
    "_item" = rep(1:3, each = 2)
    , "_annotator" = rep(1:3, times = 2)
    , "_label" = rep(1:2, each = 3)
  )

  tmp <- create_em_indixes(td, .echo = function(...) NULL)

  # equal inital propbabilities
  prev_prior <- array(.5, 2, list(c("1", "2")))
  abl_prior <- .5

  res <- init_em_params(tmp, .ability.prior = abl_prior, .prevalence.prior = prev_prior)

  # test output
  expect_equal(length(res), 3L)
  expect_identical(names(res), c("theta_hat", "pi_hat", "E_z"))
  expect_true(all(vapply(res, is.array, NA)))

  expect_equal(length(dim(res$theta_hat)), 3L)
  expect_identical(dim(res$theta_hat), c(2L, 2L, 3L))
  expect_true(all(res$theta_hat == abl_prior))

  expect_equal(length(dim(res$pi_hat)), 1L)
  expect_equal(dim(res$pi_hat), 2L)
  expect_identical(res$pi_hat, prev_prior)

  expect_equal(length(dim(res$E_z)), 2L)
  expect_identical(dim(res$E_z), c(3L, 2L))
  expect_true(all(res$E_z == .5))

  # unequal inital propbabilities
  prev_prior <- array(c(.8, .2), 2, list(c("1", "2")))
  abl_prior <- .8

  res <- init_em_params(tmp, .ability.prior = abl_prior, .prevalence.prior = prev_prior)

  # all on-diagonal elements of annotator ability matrixes = .ability.prior
  expect_true(all(apply(res$theta_hat, 3, diag) == abl_prior))
  # all off-diagonal elements of annotator ability matrixes = (1-.ability.prior)/(K-1)
  expect_true(all(apply(res$theta_hat, 3, function(x) x[upper.tri(x)]) == 1-abl_prior/(2-1)))
})

# test em_iterate() ----
test_that("em_iterate()", {

  td <- tibble(
    "_item" = rep(1:3, each = 2)
    , "_annotator" = rep(1:3, times = 2)
    , "_label" = rep(1:2, each = 3)
  )

  x <- create_em_indixes(td, .echo = function(...) NULL)

  # equal inital propbabilities
  prev_prior <- array(.5, 2, list(c("1", "2")))
  abl_prior <- .5

  p <- init_em_params(x, .ability.prior = abl_prior, .prevalence.prior = prev_prior)

  # first iteration
  res <- em_iterate(
    x = x
    , p = p
    , max.iters = 1
    , beta.prior = 0.01
    , alpha.prior = 0.01
    , min.relative.diff = 0.01
    , .echo = function(...) NULL
  )

  expect_is(res, "list")
  expect_equal(length(res), 4L)
  expect_identical(names(res), c("theta_hat", "pi_hat", "E_z", "iter_log"))

  expect_equal(length(res$iter_log), 1L)
  expect_equal(res$iter_log[[1]]$relative_diff, NA_real_)

  expect_identical(res$theta_hat, p$theta_hat)
  expect_identical(res$pi_hat, p$pi_hat)
  expect_identical(res$E_z, p$E_z)

  # all (2) iterations
  res <- em_iterate(
    x = x
    , p = p
    , max.iters = 10
    , beta.prior = 0.01
    , alpha.prior = 0.01
    , min.relative.diff = 0.01
    , .echo = function(...) NULL
  )

  expect_equal(length(res$iter_log), 2L)
  expect_equal(res$iter_log[[2]]$relative_diff, 0)

  expect_identical(res$theta_hat, p$theta_hat)
  expect_identical(res$pi_hat, p$pi_hat)
  expect_identical(res$E_z, p$E_z)

  # unequal inital propbabilities
  prev_prior <- array(c(.8, .2), 2, list(c("1", "2")))
  abl_prior <- .8

  p <- init_em_params(x, .ability.prior = abl_prior, .prevalence.prior = prev_prior)

  # first iteration
  res <- em_iterate(
    x = x
    , p = p
    , max.iters = 1
    , beta.prior = 0.01
    , alpha.prior = 0.01
    , min.relative.diff = 0.01
    , .echo = function(...) NULL
  )

  # annotator ability estimates
  expect_equivalent(round(res$theta_hat[,,1], 1), array(c(.6, .1, .4, .9), list(2, 2)))
  expect_equivalent(round(res$theta_hat[,,2], 1), array(c(.8, 0, .2, 1), list(2, 2)))
  expect_equivalent(round(res$theta_hat[,,3], 1), array(c(.8, .2, .2, .8), list(2, 2)))

  # prevalence estimates
  expect_equivalent(round(res$pi_hat, 2), array(c(.66, .34)))

  # items' max label class probability estimates
  expect_identical(round(diag(res$E_z[1:3, apply(res$E_z, 1, which.max)]), 2), c(.98, .8, .8))
  expect_true(all(rowSums(res$E_z) == 1))

  # all iterations
  res <- em_iterate(
    x = x
    , p = p
    , max.iters = 10
    , beta.prior = 0.01
    , alpha.prior = 0.01
    , min.relative.diff = 0.01
    , .echo = function(...) NULL
  )

  n_iters = length(res$iter_log)
  expect_true(res$iter_log[[n_iters]]$relative_diff < .01)

  # summing by rows over annotator ability matrix, and then by matrix
  expect_true(all(rowSums(res$theta_hat) == x$J))
  expect_true(all(diag(res$theta_hat[,,1]) > 0.5))
  # label class prevalence estimates
  expect_true(res$pi_hat[1] > res$pi_hat[2])
  # item label class estimates
  expect_true(all(apply(res$E_z, 1, which.max) == c(1, 1, 2)))
  expect_true(all(rowSums(res$E_z) == 1))

  # with DAwid-Skene pacakge dataset
  dm <- dawidskene %>%
    ungroup() %>%
    group_by(patient) %>%
    mutate(`_item` = cur_group_id()) %>%
    group_by(observer) %>%
    mutate(`_annotator` = cur_group_id()) %>%
    group_by(diagnosis) %>%
    mutate(`_label` = cur_group_id()) %>%
    ungroup()


  x <- create_em_indixes(dm, .echo = function(...) NULL)
  prev_prior <- as.array(prop.table(table(dm[["_label"]])))

  res <- em_iterate(
    x = x
    , p = init_em_params(x, .5, prev_prior)
    , max.iters = 10
    , beta.prior = 0.01
    , alpha.prior = 0.01
    , min.relative.diff = 0.0001
    , .echo = function(...) NULL
  )

  n_iters <- length(res$iter_log)
  expect_true(res$iter_log[[n_iters]]$relative_diff < .0001)

  # summing by rows over annotator ability matrix, and then by matrix
  expect_true(all(rowSums(res$theta_hat) == x$J))
  expect_true(all(diag(res$theta_hat[,,1]) > 0.5))
  expect_true(all(apply(res$theta_hat, 3, diag)[1, ] > .5))

  # label class prevalence estimates
  expect_true(which.max(res$pi_hat) == 2L)
  # item label class estimates
  expect_true(all(1:4 %in% apply(res$E_z, 1, which.max)))
  expect_true(all(rowSums(res$E_z) == 1))
})


# test em() ----
test_that("em()", {
  res <- em(
    data = dawidskene
    , item.col = patient
    , annotator.col = observer
    , label.col = diagnosis
    , verbose = FALSE
  )

  # check class
  expect_is(res, "em.fit")
  expect_is(res, "list")
  expect_identical(class(res), c("em.fit", "list"))

  expect_is(res$call, "call")
  expect_equal(as.list(res$call)[[1]], as.name("em"))
  expect_equal(as.list(res$call)$data, as.name("dawidskene"))
  expect_equal(as.list(res$call)$item.col, as.name("patient"))
  expect_equal(as.list(res$call)$annotator.col, as.name("observer"))
  expect_equal(as.list(res$call)$label.col, as.name("diagnosis"))

  # --- items' class probability estimates --- #
  # is data frame
  expect_is(res$est_class_probs, "data.frame")
  # first column named like `item.col`
  expect_equal(names(res$est_class_probs)[1], as.character(as.list(res$call)$item.col))
  # estimates for all items
  expect_true(all(res$est_class_probs$patient %in% dawidskene$patient))
  expect_true(all(dawidskene$patient %in% res$est_class_probs$patient))

  # one column for every label class
  expect_identical(names(res$est_class_probs)[-c(1, ncol(res$est_class_probs))], names(table(dawidskene$diagnosis)))
  # label class estimates sum to ~1
  expect_true(all(rowSums(res$est_class_probs[, -c(1, ncol(res$est_class_probs))]) == 1))

  # majority vote labelings exist
  expect_identical(names(res$est_class_probs)[ncol(res$est_class_probs)], "majority_vote")
  # majority vote labelings in input data label classes
  expect_true(all(res$est_class_probs$majority_vote %in% unique(dawidskene$diagnosis)))

  # --- class prevalence estimates --- #
  # is data frame
  expect_is(res$est_class_prevl, "data.frame")
  # first column named like `label.col`
  expect_equal(names(res$est_class_prevl)[1], as.character(as.list(res$call)$label.col))
  # one estimate per label class
  expect_equal(nrow(res$est_class_prevl), length(unique(dawidskene$diagnosis)))
  # all estimates/proportions sum to 1
  expect_identical(vapply(res$est_class_prevl[,-1], typeof, NA_character_, USE.NAMES = FALSE), rep("double", 3))
  expect_true(all(colSums(res$est_class_prevl[,-1]) == 1))

  # --- annotator ability estimates --- #
  # is data frame  with four columns
  expect_is(res$est_annotator_params, "data.frame")
  expect_length(res$est_annotator_params, 4L)
  # first column named like `annotator.col`
  expect_equal(names(res$est_annotator_params)[1], as.character(as.list(res$call)$annotator.col))
  # second column named like `label.col`
  expect_equal(names(res$est_annotator_params)[2], as.character(as.list(res$call)$label.col))
  # third column named 'labeled'
  expect_equal(names(res$est_annotator_params)[3], "labeled")
  # fourth column named 'est_prob', type 'double'
  expect_equal(names(res$est_annotator_params)[4], "est_prob")
  expect_type(res$est_annotator_params[[4]], "double")

  # verbosity
  expect_message(
    em(
      data = dawidskene
      , item.col = patient
      , annotator.col = observer
      , label.col = diagnosis
      , max.iters = 1
      , verbose = TRUE
    )
    , c("Data input")
  )

})
