\name{objective_function}

\alias{objective_function}

\title{Generate an objective function for BioCro model validation}

\description{
  Given a base model definition, drivers to run the model, observed values of
  model outputs, and the names of model arguments to vary,
  \code{objective_function} creates an objective function that can be used to
  parameterize the model.

  The objective function itself is based on a weighted least-squares error
  metric, with optional user-defined penalty terms, and an optional
  regularization penalty term.

  It is possible to define a multi-year or multi-location objective function by
  pairing particular sets of drivers with corresponding sets of observed model
  outputs.

  It is also possible to include "dependent" model arguments, whose values are
  determined from the "independent" model arguments that are varied during the
  parameterization procedure.
}

\usage{
  objective_function(
    base_model_definition,
    data_driver_pairs,
    independent_args,
    quantity_weights,
    data_definitions = list(),
    normalization_method = 'mean_max',
    dependent_arg_function = NULL,
    post_process_function = NULL,
    extra_penalty_function = NULL,
    regularization_method = 'none'
  )
}

\arguments{
  \item{base_model_definition}{
    A list meeting the requirements for BioCro
    \code{\link[BioCro]{crop_model_definitions}}.
  }

  \item{data_driver_pairs}{
    A list of named elements, where each element is a "data-driver pair." A
    data-driver pair is a list with two elements, \code{data} and
    \code{drivers}. The \code{data} element must be a data frame with one column
    named \code{time}, whose values follow BioCro's definition of
    \code{\link[BioCro]{time}}; the other columns should represent observed
    values of model outputs. The \code{drivers} element must be a data frame
    that can be passed to \code{\link[BioCro]{run_biocro}} as its \code{drivers}
    input argument.
  }

  \item{independent_args}{
    A list of named numeric values. The names will determine the independent
    arguments to be varied during their optimization, and the values specify
    "test" values of each argument that will be used internally to check that
    the objective function is properly defined and can be evaluated.
  }

  \item{quantity_weights}{
    A list of named numeric values, where the name of each element is one of the
    model outputs to be compared against the observed data, and the value is the
    weight for that output.
  }

  \item{data_definitions}{

  }

  \item{normalization_method}{

  }

  \item{dependent_arg_function}{

  }

  \item{post_process_function}{

  }

  \item{extra_penalty_function}{

  }

  \item{regularization_method}{

  }
}

\details{

}

\value{
  A function \code{obj_fun} with signature \code{obj_fun(x, lambda = 0)}. Here,
  \code{x} is a numeric vector of values of the independent arguments (in the
  same order as in \code{independent_arg_names}), and \code{lambda} is the value
  of the regularization parameter. The \code{obj_fun} function returns values of
  the error metric, as described above.
}

\examples{
# Example: Create an objective function that enables optimization of the
# `alphaLeaf`, `betaLeaf`, and `alphaStem` parameters of the Soybean-BioCro
# model. Additional details are provided below.


if (require(BioCro)) {
  # We will use Soybean-BioCro as the base model definition, but we will change
  # the ODE solver to use the Euler method so the model runs faster.
  base_model_definition            <- BioCro::soybean
  base_model_definition$ode_solver <- BioCro::default_ode_solvers[['homemade_euler']]

  # We will use the `soyface_biomass` data set (included with the
  # `BioCroValidation` package) for the observed values; this set includes
  # observations of leaf, stem, and pod biomass from two years, which are stored
  # in two data tables. However, these data tables each have a `DOY` column
  # rather than a `time` column, so we need to alter them.
  ambient_2002_data <- within(soyface_biomass[['ambient_2002']], {
    time = (DOY - 1) * 24.0 # Define new `time` column
    DOY  = NULL             # Remove unneeded `DOY` column
  })

  ambient_2005_data <- within(soyface_biomass[['ambient_2005']], {
    time = (DOY - 1) * 24.0 # Define new `time` column
    DOY  = NULL             # Remove unneeded `DOY` column
  })

  # The data-driver pairs can now be created by associating each data set with
  # its corresponding weather data.
  data_driver_pairs <- list(
    ambient_2002 = list(
      data    = ambient_2002_data,
      drivers = BioCro::soybean_weather[['2002']]
    ),
    ambient_2005 = list(
      data    = ambient_2005_data,
      drivers = BioCro::soybean_weather[['2005']]
    )
  )

  # In the data, the leaf biomass is in the `Leaf_Mg_per_ha` column, but in the
  # simulation output, it is in the `Leaf` column. Similar naming differences
  # occur for the stem and pod mass. To address this, we can provide a data
  # definition list.
  data_definitions <- list(
    Leaf_Mg_per_ha = 'Leaf',
    Stem_Mg_per_ha = 'Stem',
    Pod_Mg_per_ha = 'Pod'
  )

  # The data contains values of pod mass, but the model does not calculate pod
  # mass; instead, it returns separate values of `Grain` (seed) and `Shell`
  # mass, two components which form the pod together. To address this, we can
  # provide a post-processing function to calculate the pod mass.
  post_process_function <- function(x) {
    within(x, {Pod = Grain + Shell})
  }

  # Here we wish to independently vary the `alphaLeaf` and `betaLeaf`
  # parameters. We also wish to vary `alphaStem`, but require that its value is
  # always equal to `alphaLeaf`. To do this, we can specify independent
  # arguments, and a function for determining dependent argument values. We will
  # choose "test" values of the independent arguments as their values in the
  # original Soybean-BioCro model.
  independent_args <- BioCro::soybean[['parameters']][c('alphaLeaf', 'betaLeaf')]

  dependent_arg_function <- function(x) {
    list(alphaStem = x[['alphaLeaf']])
  }

  # When determining the error metric value, we wish to weight the pod twice as
  # large as the leaf and stem to ensure a close fit to the observed pod masses.
  quantity_weights <- list(
    Leaf = 0.5,
    Stem = 0.5,
    Pod = 1
  )

  # Now we can finally create the objective function.
  obj_fun <- objective_function(
    base_model_definition,
    data_driver_pairs,
    independent_args,
    quantity_weights,
    data_definitions = data_definitions,
    post_process_function = post_process_function
  )

  # This function could now be passed to an optimizer; here we will simply
  # evaluate it for two sets of parameter values.

  # First try the initial values.
  obj_fun(as.numeric(independent_args))

  # Now try doubling each parameter value; in this case, the value of the
  # objective function increases, indicating a lower degree of agreement between
  # the model and the observed data.
  obj_fun(2 * as.numeric(independent_args))
}
}
