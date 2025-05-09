\name{objective_function}

\alias{objective_function}

\title{Generate an objective function for BioCro model validation}

\description{
  Given a base model definition, drivers to run the model, observed values of
  model outputs, and the names of model arguments to vary,
  \code{objective_function} creates an objective function that can be passed to
  a minimization algorithm in order to find optimal parameter values that
  produce the best agreement between the model and the observed data.

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
    regularization_method = 'none',
    verbose_startup = FALSE
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
    weight for that output. Each weight can be a single number, or a pair of
    numbers. When the weight is a pair, the first number is the weight that will
    be used for underestimates (when the modeled value is smaller than the
    observed value), and the second is the weight for overestimates.
  }

  \item{data_definitions}{
    A list of named string values, where the name of each element is one of the
    data columns in the \code{data_driver_pairs} and the value is that column's
    corresponding name in the model output. For all other data columns in the
    \code{data_driver_pairs}, it is assumed that the data column name matches a
    column in the model output.
  }

  \item{normalization_method}{
    A string indicating the normalization method to be used when calculating the
    error metric; see below for more details.
  }

  \item{dependent_arg_function}{
    A function whose input argument is a named list of independent argument
    values, and which returns a named list of dependent argument values. If the
    \code{dependent_arg_function} is \code{NULL}, no dependent argument values
    will be calculated.
  }

  \item{post_process_function}{
    A function whose input argument is a data frame representing the output from
    \code{\link[BioCro]{run_biocro}}, and which returns a data frame, typically
    based on the input but with one or more new columns. If the
    \code{post_process_function} is \code{NULL}, no post-processing will be
    applied to the raw simulation output.
  }

  \item{extra_penalty_function}{
    A function whose input argument is a data frame representing the output from
    \code{\link[BioCro]{run_biocro}}, and which returns a numeric penalty to be
    added to the least-squares term when calculating the error metric. If the
    \code{extra_penalty_function} is \code{NULL}, no extra penalties will be
    added.
  }

  \item{regularization_method}{
    A string indicating the regularization method to be used when calculating
    the regularization penalty term; see below for more details.
  }

  \item{verbose_startup}{
    A logical (\code{TRUE} or \code{FALSE}) value indicating whether to print
    additional information to the R terminal when creating the objective
    function.
  }
}

\details{
  \strong{Overview}

  When parameterizing a BioCro model, the general idea is to vary a subset of
  the model's parameters to achieve the best agreement with a set of observed
  data. The degree of agreement is expressed as an "error metric," which
  includes terms derived from the agreement between the modeled and observed
  values, as well as (optional) penalty terms. A function that calculates an
  error metric when given a set of parameter values is called an "objective
  function." Defining an objective function suitable for BioCro parameterization
  can be complicated, but the \code{objective_function} function helps to
  simplify the process of creating such a function. It is designed to
  accommodate the following scenarios, which often occur in the context of
  BioCro model parameterization:

  \itemize{
    \item \strong{Multi-year or multi-location data:} Often the model needs to
          be run several times with different drivers corresponding to multiple
          years or locations, and then the results from each individual
          simulation must be compared to associated sets of observed data. Here,
          this is handled through the \code{data_driver_pairs}, which allows the
          user to specify which drivers and data sets should be compared to each
          other.

    \item \strong{Complicated normalization:} Care must be taken to ensure that
          certain years or output variables are not over-valued in the error
          metric; for example, one year may have more observations of leaf mass
          than another year, or the stem mass may be much larger than the leaf
          mass. Here, this is handled through pre-set normalization approaches,
          which can be specified through the \code{normalization_method} input.
          See below for more information.

    \item \strong{Extra penalties:} Sometimes an optimizer chooses parameters
          that produce close agreement with the observed data, but are
          nevertheless not biologically resonable. For example, it may produce
          a sharp peak at a high leaf mass in between two measured points, when
          in reality, leaf mass should vary monotonically between them. In this
          case, it may be necessarily to add an extra penalty to the objective
          function that prevents the optimizer from choosing such values. Here,
          this is handled through the \code{extra_penalty_function} input.

    \item \strong{Flexible weights:} Typically a user would like to specify a
          weight for each variable being considered in the error metric, either
          to represent an uncertainty or to emphasize agreement with one output
          at the expense of another. For example, the seed mass may need a high
          weight to prioritize accurate yield predictions. Further, these
          weights may need to differ for underestimates as compared to
          overestimates; for example, measured root mass is often lower than the
          true root mass, so a user may wish to penalize underestimates of root
          mass more severely than overestimates. Here, this is handled through
          the \code{quantity_weights} input.

    \item \strong{Dependent parameters:} Sometimes one model parameter must be
          determined from one or more other parameters; for example, a user may
          wish to require that the leaf and stem maintenance respiration
          coefficients are identical. Here, this is handled through the
          \code{dependent_arg_function}, which allows the user to specify how
          the values of such "dependent" parameters should be determined from
          the values of the "independent" parameters.

    \item \strong{Name mismatches:} Often a particular variable has different
          names in the data set and simulation output. Here, this is handled
          through the \code{data_definitions}, which allows the user to specify
          which columns in the model output should be compared to particular
          columns in the observed data.

    \item \strong{Incomplete outputs:} Sometimes a model may not produce outputs
          that are directly comparable to the observed values; for example, a
          model may calculate seed and shell mass, while a data set includes pod
          mass, which is the sum of seed and shell. Here, this is handled by an
          optional \code{post_process_function}, which allows users to specify
          additional operations to perform on the model output; in this example,
          it would be used to calculate the pod mass so it can be compared to
          the observations.
  }

  \strong{Error metric calculations}

  As mentioned above, the overall error metric \eqn{E} is calculated as

  \deqn{E = E_{agreement} + P_{user} + P_{regularization},}
  where \eqn{E_{agreement}} is determined by the agreement between the model and
  observations, \eqn{P_{user}} is an optional user-defined penalty, and
  \eqn{P_{regularization}} is an optional regularization penalty. These terms
  are explained in more detail below:

  \itemize{
    \item \strong{Agreement term:} The agreement term \eqn{E_{agreement}} is
          calculated using a least-squares approach. In other words,

          \deqn{E_{agreement} =%
          \sum_i \frac{\left(Y_{obs}^i - Y_{mod}^i \right)^2}{N_i} \cdot w_i,}
          where the sum runs over all \eqn{n} observations, \eqn{Y_{obs}^i} and
          \eqn{Y_{mod}^i} are observed and modeled values of variable \eqn{Y_i},
          \eqn{N_i} is a normalization factor that depends on the name of
          \eqn{Y_i} and the data set that includes the \eqn{i^{th}} observation,
          and \eqn{w_i} is a weight factor that depends on the name of \eqn{Y}.

          Each value of \eqn{Y_{obs}^i} is specified at a particular time
          \eqn{t_i}. The corresponding modeled value, \eqn{Y_{mod}^i}, is found
          by retrieving the value of the \eqn{Y_i} variable at the closest time
          to \eqn{t_i} that is included in the model output. It is assumed that
          the model always outputs the same sequence of time values each time it
          is run with a particular set of drivers, regardless of the input
          argument values.

          The normalization factors \eqn{N_i} are determined by the choice of
          \code{normalization_method}; the available methods are discussed
          below.

          The weight factors \eqn{w_i} are directly specified by the user via
          the \code{quantity_weights} input. For example, if
          \code{quantity_weights} has an element named \code{Leaf} that is equal
          to 0.5, then \eqn{w_i} will be equal to 0.5 whenever \eqn{Y_i}
          represents a leaf mass value, regardless of which set of drivers or
          time point corresponds to \eqn{Y_i}. The weights can also be supplied
          as \eqn{(w_{under}, w_{over})} pairs instead of single values; in this
          case, the value of \eqn{w_i} depends on whether the model makes an
          underprediction or an overprediction: \eqn{w_i = w_{under}} when
          \eqn{Y_{mod}^i < Y_{obs}^i} and \eqn{w_i = w_{over}} otherwise.

          There are a few special cases where \eqn{E_{agreement}} is set to a
          very high value (\code{BioCroValidation:::FAILURE_VALUE}). This is
          done when a simulation fails to run, when the \eqn{E_{agreement}} term
          would otherwise evaluate to \code{NA}, or when the \eqn{E_{agreement}}
          term would otherwise evaluate to an infinite value.

    \item \strong{User-defined penalty term:} The user-defined penalty term
          \eqn{P_{user}} is calculated by applying a function \eqn{f_{user}} to
          the full simulation output from each set of drivers. In other words,

          \deqn{P_{user} = \sum_k f_{user} \left( M_k \right),}
          where the sum runs over all \eqn{k} sets of drivers and \eqn{M_k} is
          the model output when it is run with the \eqn{k^{th}} set of drivers.

          The function \eqn{f_{user}} must accept a single data frame as an
          input and return a single numeric value as its output, but has no
          other requirements. It is specified via the
          \code{extra_penalty_function} argument. When
          \code{extra_penalty_function} is \code{NULL}, \eqn{P_{user}} is zero.

    \item \strong{Regularization penalty term:} The regularization penalty term
          \eqn{P_{regularization}} is calculated from the values of the
          arguments being varied during the optimization by applying a function
          \eqn{R}. In other words,

          \deqn{P_{regularization} = R \left( X \right),}
          where \eqn{X} represents the model argument values.

          The function \eqn{R} is determined by the choice of
          \code{regularization_method}; the available methods are discussed
          below.
  }

  \strong{Normalization methods}

  The following normalization methods are available:

  \itemize{
    \item \code{'none'}: For this method, \eqn{N_i} is always set to 1. In other
          words, no normalization is performed.

    \item \code{'mean'}: For this method, when \eqn{Y_i} is named \code{vtype}
          and the observation is from a set called \code{vdata}, then

          \deqn{N_i = n_{vtype}^{vdata},}
          where \eqn{n_{vtype}^{vdata}} is the number of observations of type
          \code{vtype} that are included in \code{vdata}. In this case, the
          error term \eqn{E_{agreement}} becomes a mean error for each set of
          drivers, hence the name for this method. This approach avoids
          over-representing drivers with larger numbers of associated
          observations when determining \eqn{E_{agreement}}.

    \item \code{'max'}: For this method, when \eqn{Y_i} is named \code{vtype}
          and the observation is from a set called \code{vdata}, then

          \deqn{N_i = \left( max_{vtype}^{vdata} \right)^2,}
          where \eqn{max_{vtype}^{vdata}} is the maximum observed value of
          \code{vtype} across \code{vdata}. In this case, the observed and
          modeled values that appear in the equation for \eqn{E_{agreement}} are
          essentially normalized by their maximum value, hence the name for this
          method. This approach avoids over-representing variable types with
          larger magnitude when determining \eqn{E_{agreement}}.

    \item \code{'mean_max'}: For this method, the "mean" and "max" methods are
          combined so that

          \deqn{N_i =%
          n_{vtype}^{vdrivers} \cdot \left( max_{vtype}^{vdata} \right)^2.}
          This approach avoids over-representing drivers with larger numbers of
          associated observations, and variable types with larger magnitudes.
          This method is used for parameterizing Soybean-BioCro, and is
          recommended for most situations.
  }

  \strong{Regularization methods}

  The following regularization methods are available:

  \itemize{
    \item \code{'none'}: For this method, \eqn{P_{regularization}} is always set
          to 0. In other words, no regularization is performed.

    \item \code{'L1'} or \code{'lasso'}: For this method,
          \eqn{P_{regularization}} is given by the sum of the absolute values of
          each independent argument, multiplied by a "regularization parameter"
          \eqn{\lambda} that sets the overall weight of the penalty:

          \deqn{P_{regularization} = \lambda \sum_j | X_j |,}
          where the sum runs over all \eqn{j} independent arguments, and
          \eqn{X_j} is the value of the \eqn{j^{th}} argument. See the "Value"
          section below for details of how to specify \eqn{\lambda}.

    \item \code{'L2'} or \code{'ridge'}: For this method,
          \eqn{P_{regularization}} is given by the sum of the squared values of
          each independent argument, multiplied by a "regularization parameter"
          \eqn{\lambda} that sets the overall weight of the penalty:

          \deqn{P_{regularization} = \lambda \sum_j X_j^2,}
          where the sum runs over all \eqn{j} independent arguments, and
          \eqn{X_j} is the value of the \eqn{j^{th}} argument. See the "Value"
          section below for details of how to specify \eqn{\lambda}.
  }

  \strong{Input checks}

  Several checks are made to ensure that the objective function is properly
  defined. These checks include, but are not limited to, the following:

  \itemize{
    \item Ensuring that each set of drivers in \code{data_driver_pairs} defines
          a valid dynamical system along with the \code{base_model_definition}.
          This is accomplished using
          \code{\link[BioCro]{validate_dynamical_system_inputs}}.

    \item Ensuring that the model output corresponding to each set of drivers
          spans the times at which the observations were made.

    \item Ensuring that each variable type in the data elements of
          \code{data_driver_pairs} matches a corresponding column in the model
          output, when accounting for the \code{data_definitions} and
          \code{post_process_function}.

    \item Ensuring that each independent and dependent argument name is either
          a parameter or initial value of the model. This is accomplished using
          \code{\link[BioCro]{partial_run_biocro}}. Note that argument names
          passed to \code{partial_run_biocro} can technically include drivers,
          but it is unlikely that the value of a driver would be varied during
          an optimization, so the argument names are not allowed to include
          columns in the drivers.

    \item Ensuring that the optional \code{dependent_arg_function},
          \code{post_process_function}, and \code{extra_penalty_function}
          functions can be run without causing errors.
  }

  If any issues are detected, an informative error message will be sent. Note
  that several of these checks require running the model with each set of
  drivers. For these checks, the argument values specified by
  \code{independent_args} will be used, so they should be valid or otherwise
  reasonable values.
}

\value{
  A function \code{obj_fun} with signature
  \code{obj_fun(x, lambda = 0, return_terms = FALSE)}.

  Here, \code{x} is a numeric vector of values of the independent arguments (in
  the same order as in \code{independent_arg_names}), and \code{lambda} is the
  value of the regularization parameter.

  The \code{return_terms} argument determines the return value of
  \code{obj_fun}. When \code{return_terms} is \code{FALSE}, \code{obj_fun}
  returns values of the error metric \eqn{E}. When \code{return_terms} is
  \code{TRUE}, \code{obj_fun} returns a list including each individual term of
  the total error metric.

  During optimization, \code{return_terms} should always be \code{FALSE}.
  Setting it to \code{TRUE} can be useful for troubleshooting, or for
  diagnostics such as checking the quality of fit for each of the data-driver
  pairs.
}

\examples{
# Example: Create an objective function that enables optimization of the
# `alphaLeaf`, `betaLeaf`, and `alphaStem` parameters of the Soybean-BioCro
# model. Additional details are provided below. Important note: This example is
# designed to highlight key features of `objective_function`, and is not
# necessarily realistic.

if (require(BioCro)) {
  # We will use Soybean-BioCro as the base model definition, but we will change
  # the ODE solver to use the Euler method so the model runs faster.
  base_model_definition            <- BioCro::soybean
  base_model_definition$ode_solver <- BioCro::default_ode_solvers[['homemade_euler']]

  # We will use the `soyface_biomass` data set (included with the
  # `BioCroValidation` package) for the observed values; this set includes
  # observations of leaf, stem, and pod biomass from two years, which are stored
  # in two data tables; it also includes the standard deviations of the measured
  # biomasses, which are included in two separate tables. However, these data
  # tables each have a `DOY` column rather than a `time` column, so we need to
  # alter them. Here we will do this by defining a short helper function.
  convert_time <- function(x) {
    within(x, {
      time = (DOY - 1) * 24.0 # Define new `time` column
      DOY  = NULL             # Remove unneeded `DOY` column
    })
  }

  # The data-driver pairs can now be created by associating each data set with
  # its corresponding weather data. Here we will weight the 2005 data twice as
  # heavily as the 2002 data.
  data_driver_pairs <- list(
    ambient_2002 = list(
      data       = convert_time(soyface_biomass[['ambient_2002']]),
      data_stdev = convert_time(soyface_biomass[['ambient_2002_std']]),
      drivers    = BioCro::soybean_weather[['2002']],
      weight     = 1
    ),
    ambient_2005 = list(
      data       = convert_time(soyface_biomass[['ambient_2005']]),
      data_stdev = convert_time(soyface_biomass[['ambient_2005_std']]),
      drivers    = BioCro::soybean_weather[['2005']],
      weight     = 2
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
  post_process_function <- function(sim_res) {
    within(sim_res, {Pod = Grain + Shell})
  }

  # Here we wish to independently vary the `alphaLeaf` and `betaLeaf`
  # parameters. We also wish to vary `alphaStem`, but require that its value is
  # always equal to `alphaLeaf`. To do this, we can specify independent
  # arguments, and a function for determining dependent argument values. We will
  # choose "test" values of the independent arguments as their values in the
  # original Soybean-BioCro model.
  independent_args <- BioCro::soybean[['parameters']][c('alphaLeaf', 'betaLeaf')]

  dependent_arg_function <- function(ind_args) {
    list(alphaStem = ind_args[['alphaLeaf']])
  }

  # When determining the error metric value, we wish to weight the pod highest
  # to ensure a close fit to the observed pod masses. We also wish to decrease
  # the penalty for overestimates of the stem mass, since we believe our
  # observations to be underestimates.
  quantity_weights <- list(
    Leaf = 0.5,
    Stem = c(0.5, 0.25),
    Pod = 1
  )

  # We want to prevent the optimizer from choosing parameters that produce
  # unreasonably high leaf mass
  extra_penalty_function <- function(sim_res) {
    max_leaf <- max(sim_res[['Leaf']], na.rm = TRUE)

    if (is.na(max_leaf) || max_leaf > 4) {
      1e5 # Add a steep penalty
    } else {
      0
    }
  }

  # Now we can finally create the objective function; here we choose to print
  # extra startup information, but this is not strictly required.
  obj_fun <- objective_function(
    base_model_definition,
    data_driver_pairs,
    independent_args,
    quantity_weights,
    data_definitions = data_definitions,
    post_process_function = post_process_function,
    extra_penalty_function = extra_penalty_function,
    verbose_startup = TRUE
  )

  # This function could now be passed to an optimizer; here we will simply
  # evaluate it for two sets of parameter values.

  # First try the initial values.
  cat('\nError metric calculated from initial argument values:\n')
  print(
    obj_fun(as.numeric(independent_args))
  )

  # Now try doubling each parameter value; in this case, the value of the
  # objective function increases, indicating a lower degree of agreement between
  # the model and the observed data.
  cat('\nError metric calculated by doubling the initial argument values:\n')
  print(
    obj_fun(2 * as.numeric(independent_args))
  )

  # We can also see the values of each term that makes up the error metric
  cat('\nError metric terms calculated by doubling the initial argument values:\n')
  str(
    obj_fun(2 * as.numeric(independent_args), return_terms = TRUE)
  )
}
}
