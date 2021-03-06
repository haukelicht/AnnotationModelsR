% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/em.R
\name{init_em_params}
\alias{init_em_params}
\title{Initialize EM-parameters for binar model}
\usage{
init_em_params(x, .ability.prior, .prevalence.prior)
}
\arguments{
\item{x}{a list as returned by \code{\link{create_em_indixes}}}

\item{.ability.prior}{number in interval (0, 1), specifying the annotator ability prior used to initialize ability parameters for first iteration
Defaults to 0.7}

\item{.prevalence.prior}{vector of numbers in interval (0, 1),
specifying prior probabilities on item label classes \eqn{\hat{\boldsymbol{pi}}}.
Most have length equal to number of label classes (i.e., number of unique item labels in \code{data}).
Defaults to \code{NULL}, in which case prior class label probabilities are set all equal to \eqn{1/K} (where \eqn{K} is the number of item label classes)}
}
\value{
a list with elements
    \itemize{
      \item 'theta_hat': annotator parameter matrices (a 3d-tensor with extension \eqn{K* K* J})
      \item 'pi_hat': prevalence estimates (a \eqn{K}-length numeric vector)
      \item 'E_z': \eqn{I* K} matrix of item label class probability estimates
    }
}
\description{
Initialize EM-parameters for binar model
}
