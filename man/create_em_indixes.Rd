% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/em.R
\name{create_em_indixes}
\alias{create_em_indixes}
\alias{em_iterate}
\title{Create EM indexes for binar model}
\usage{
create_em_indixes(x, .echo)

em_iterate(x, p, max.iters, beta.prior, alpha.prior, min.relative.diff, .echo)
}
\arguments{
\item{x}{a \code{data.frame} with columns '_item', '_annotator' and '_label'
holding corresponding integer indexes}

\item{.echo}{a function used to print messages}

\item{p}{a list of initialized parameters as returned by \code{\link{init_em_params}}}

\item{max.iters}{positive integer, specifying the maximum number of EM iterations to perform.
(CAUTION: low values may cause stopping before convergence)
Defaults to 100.}

\item{beta.prior}{number in interval (0, 1), specifying the beta prior used to perform smoothing on \eqn{\hat{pi}} (label class prevalence estimates)}

\item{alpha.prior}{number in interval (0, 1), specifying the alpha prior used to perform smoothing on items' label class probability estimates}

\item{min.relative.diff}{number in interval (0, 1), specifying the minimum relative difference between log-likelihoods of two subsequent iteration steps.
Governs when convergence is reached (the lower, the more iterations).
Defaults to 1e-8}
}
\value{
a list with elements
   \itemize{
     \item 'ii': item index values
     \item 'I': number of items
     \item 'jj': annotator index values
     \item 'J': number of annotators
     \item 'y': label index values
     \item 'K': number of label classes
     \item 'N': number of annotations
   }

a list with elements
    \itemize{
      \item 'theta_hat': annotator parameter matrices (a 3d-tensor with extension \eqn{K* K* J})
      \item 'pi_hat': prevalence estimates (a \eqn{K}-length numeric vector)
      \item 'E_z': \eqn{I* K} matrix of item label class probability estimates
      \item 'iter_log': a list with as many elements as iterations until convergence.
          Each element is a list with elements
          'idx' (the iteration counter, starts at 1),
          'log_likelihood' (the log-likelihood value of the current iteration),
          'relative_diff' (the difference in log-likelihood values of the current relative to the previous iteration)
    }
}
\description{
Create EM indexes for binar model

Expectation-Maximization iterator for binary model
}
\section{Functions}{
\itemize{
\item \code{em_iterate}: 
}}

