\name{write_model}

\alias{write_model}

\title{Generate a BioCro model definition}

\description{
  To facilitate the creation and validation of new BioCro models,
  \code{write_model} stores the contents of a BioCro crop model definition
  (see \code{\link[BioCro]{crop_model_definitions}}) in an R script to make its
  modules and parameter values human-readable.
}

\usage{
  write_model(
    name,
    direct_modules,
    differential_modules,
    initial_values,
    parameters,
    ode_solver
  )
}

\arguments{
  \item{name}{A string specifying a name for the model.}

  \item{direct_modules}{
    A list or vector of fully-qualified direct module names that can be passed
    to \code{\link[BioCro]{run_biocro}} as its \code{direct_module_names}
    input argument.
  }

  \item{differential_modules}{
    A list or vector of fully-qualified differential module names that can be
    passed to \code{\link[BioCro]{run_biocro}} as its
    \code{differential_module_names} input argument.
  }

  \item{initial_values}{
    A list of named elements that can be passed to
    \code{\link[BioCro]{run_biocro}} as its \code{initial_values} input
    argument.
  }

  \item{parameters}{
    A list of named elements that can be passed to
    \code{\link[BioCro]{run_biocro}} as its \code{parameters} input
    argument.
  }

  \item{ode_solver}{
    A list of named elements that can be passed to
    \code{\link[BioCro]{run_biocro}} as its \code{ode_solver} input
    argument.
  }
}

\details{
  This function returns a string and has no file I/O. Use
  \code{\link{writeLines}} to print the output to console, or to save the
  output. See examples below. Note that it is customary to name the R script
  file as \code{name.R}, where \code{name} is the value provided to the function
  itself.
}

\value{
  A string containing a new BioCro model definition R script.
}

\examples{
if (require(BioCro)) {
  # Create a definition string for Soybean-BioCro
  out <- with(BioCro::soybean, {write_model(
    'test_soybean_model',
    direct_modules,
    differential_modules,
    initial_values,
    parameters,
    ode_solver
  )})

  # Use writeLines to print to console
  writeLines(out)

  if (interactive()) {
    # Use writeLines to save as a `.R` file
    writeLines(out, "./test_soybean_model.h")
  }
}
}
