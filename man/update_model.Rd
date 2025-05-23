\name{update_model}

\alias{update_model}

\title{Update a BioCro model definition}

\description{
  Following an optimization, it is typically necessary to update the initial
  values and/or parameters of a base model definition with new values determined
  during the optimization. The \code{update_model} function helps to streamline
  this process. It is expected that this function will be called after
  \code{\link{objective_function}}.
}

\usage{
  update_model(
    base_model_definition,
    independent_args,
    new_ind_arg_values,
    dependent_arg_function = NULL
  )
}

\arguments{
  \item{base_model_definition}{
    The same value that was passed to \code{\link{objective_function}}.
  }

  \item{independent_args}{
    The same value that was passed to \code{\link{objective_function}}.
  }

  \item{new_ind_arg_values}{
    A numeric vector representing new values of the independent arguments,
    typically determined by an optimizer.
  }

  \item{dependent_arg_function}{
    The same value that was passed to \code{\link{objective_function}}.
  }
}

\value{
  A list based on \code{base_model_definition} but with new values of some of
  its \code{initial_values} and \code{parameters}, as specified by the elements
  of \code{independent_args} and \code{new_ind_arg_values}.
}

\examples{
if (require(BioCro)) {
  # Update the default Soybean-BioCro model with new values of `Leaf` (an
  # initial value) and `alphaStem` (a parameter)
  base_model <- BioCro::soybean

  new_model <- update_model(
    base_model,
    list(Leaf = 1, alphaLeaf = 2), # The values here will not be used
    c(1000, 2000)                  # These are the actual new values
  )

  # Compare the two models
  cat('\n\nComparing initial Leaf values:\n')
  print(base_model$initial_values$Leaf)
  print(new_model$initial_values$Leaf)

  cat('\n\nComparing alphaLeaf values:\n')
  print(base_model$parameters$alphaLeaf)
  print(new_model$parameters$alphaLeaf)
}
}
