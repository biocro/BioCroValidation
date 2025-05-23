\name{bounds_table}

\alias{bounds_table}

\title{Create a table of lower and upper bounds}

\description{
  During an optimization, it is often necessary to provide lower and upper
  bounds for the parameters that are being varied. Typically they are specified
  as numeric vectors, which often leads to confusing code, where the writer and
  reader must remember which value corresponds to each argument; for example,
  "the third element of the lower bound vector is for alphaLeaf."

  The purpose of \code{bounds_table} is to make the process of specifying bounds
  simpler and easier to follow. It is expected that this function will be called
  after \code{\link{objective_function}}.
}

\usage{
  bounds_table(
    independent_args,
    bounds_list,
    initial_ind_arg_values = NULL
  )
}

\arguments{
  \item{independent_args}{
    The same value that was passed to \code{\link{objective_function}}.
  }

  \item{bounds_list}{
    A list of named elements, where each element is a numeric vector of length
    2. The names should correspond to the independent arguments, and the values
    should indicate lower and upper bounds for the corresponding parameter (in
    any order). Any "extra" bounds (that is, bounds that do not correspond to
    any independent argument) will be ignored.
  }

  \item{initial_ind_arg_values}{
    A numeric vector of initial values for each of the independent arguments,
    supplied in the same order as in \code{independent_args}.
  }
}

\details{
  The main purpose of this function is to create vectors of lower and upper
  bounds, which are returned as the columns of a data frame. For each
  independent argument in \code{independent_args}, the bounds are supplied via
  the \code{bounds_list} input. The syntax is designed so the code calling this
  function is easy for a human to parse. (See example below.)

  It is also (optionally) possible to provide an initial guess for each
  independent argument via the \code{initial_ind_arg_values} argument. When
  provided, these will be checked to make sure they do not lie outside the
  bounds; an error will be thrown if any do lie outside the bounds. A warning
  will also be thrown if any initial guesses lie on the bounds, since this can
  be problematic for some optimizers, such as \code{\link[dfoptim]{nmkb}}.

  Some optimizers, such as \code{\link[DEoptim]{DEoptim}}, do not require an
  initial guess; in this case, there is no strong need to pass an initial guess
  to \code{bounds_table}.
}

\value{
  A data frame with three or four columns: \code{arg_name}, \code{lower},
  \code{upper}, and (optionally) \code{initial_value}.

  The \code{lower} and \code{upper} columns are the lower and upper bounds,
  determined from \code{bounds_list}. The \code{arg_name} column is the argument
  name, and the rows of the table are ordered as in \code{independent_args}.
  The \code{initial_value} column contains initial values, if they were
  provided via the \code{initial_ind_arg_values} argument.
}

\examples{
# Make a list of independent arguments; the values are not used for anything
independent_args <- list(
  alphaLeaf = 0,
  alphaRoot = 0,
  alphaStem = 0,
  betaLeaf  = 0,
  betaRoot  = 0,
  betaStem  = 0
)

# Specify bounds and initial guess for each. Note that:
#
# - The bounds will be reordered to follow the same order as the
#   `independent_args`, but the initial guess is assumed to already follow the
#   same order as the `independent_args`.
#
# - The bounds for the two extra parameters are ignored when forming the table.
#
# - The lower and upper bounds can be supplied as (upper, lower)
#   or (lower, upper) pairs.
#
b_ll <- -50 # Lower limit for beta parameters
a_ul <- 50  # Upper limit for alpha parameters

bounds <- bounds_table(
  independent_args,
  list(
    betaStem  = c(0, b_ll),
    betaRoot  = c(0, b_ll),
    betaLeaf  = c(0, b_ll),
    alphaStem = c(0, a_ul),
    alphaRoot = c(0, a_ul),
    alphaLeaf = c(0, a_ul),
    extraPar1 = c(0, 5),
    extraPar2 = c(0, 6)
  ),
  c(1, 1, 1, -1, -1, -1)
)

print(bounds)

# Now the properly-ordered lower and upper limits can be accessed as follows:
bounds$lower

bounds$lower
}
