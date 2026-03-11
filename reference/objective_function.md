# Generate an objective function for BioCro model validation

Given a base model definition, drivers to run the model, observed values
of model outputs, and the names of model arguments to vary,
`objective_function` creates an objective function that can be passed to
a minimization algorithm in order to find optimal parameter values that
produce the best agreement between the model and the observed data.

The objective function itself is based on a weighted least-squares error
metric, with optional user-defined penalty terms, and an optional
regularization penalty term.

It is possible to define a multi-year or multi-location objective
function by pairing particular sets of drivers with corresponding sets
of observed model outputs.

It is also possible to include "dependent" model arguments, whose values
are determined from the "independent" model arguments that are varied
during the parameterization procedure.

For a detailed example of using `objective_function`, see the
"Parameterizing Soybean-BioCro" vignette.

## Usage

``` r
objective_function(
    base_model_definition,
    data_driver_pairs,
    independent_args,
    quantity_weights,
    data_definitions = list(),
    normalization_method = 'mean_max',
    normalization_param = NULL,
    stdev_weight_method = 'equal',
    stdev_weight_param = NULL,
    regularization_method = 'none',
    dependent_arg_function = NULL,
    post_process_function = NULL,
    extra_penalty_function = NULL,
    verbose_startup = TRUE
  )
```

## Arguments

- base_model_definition:

  A list meeting the requirements for BioCro
  [`crop_model_definitions`](https://rdrr.io/pkg/BioCro/man/crop_model_definitions.html).

- data_driver_pairs:

  A list of named elements, where each element is a "data-driver pair."
  A data-driver pair is a list with three required elements: `data`,
  `drivers`, and `weight`. Optionally, it may also have elements named
  `data_stdev`, `initial_values`, and `parameters`.

  The `data` element must be a data frame with one column named `time`,
  whose values follow BioCro's definition of
  [`time`](https://rdrr.io/pkg/BioCro/man/time_variable.html); the other
  columns should represent observed values of model outputs. Any `NA`
  values in `data` will be ignored when calculating the error metric,
  but all non-`NA` values of all columns (except `time`) will be
  compared to the model output.

  The `drivers` element must be a data frame that can be passed to
  [`run_biocro`](https://rdrr.io/pkg/BioCro/man/run_biocro.html) as its
  `drivers` input argument.

  The `weight` element must be a single numeric value indicating a
  weight to be used when calculating the error metric.

  The optional `data_stdev` element must be a data frame with the same
  column names as `data`, and the same time values; its other entries
  should represent the standard deviation associated with each entry in
  `data`. If `data_stdev` is not supplied, all standard deviations will
  be set to 1.

  The optional `initial_values` and `parameters` elements must be named
  lists of driver-specific initial values and parameters, respectively,
  that will overwrite the default values specified in the
  `base_model_definition`. When driver-specific `initial_values` are
  provided, they must be provided for each of the data-driver pairs; the
  same rule applies for driver-specific `parameters`.

- independent_args:

  A list of named numeric values. The names will determine the
  independent arguments to be varied during their optimization, and the
  values specify "test" values of each argument that will be used
  internally to check that the objective function is properly defined
  and can be evaluated.

- quantity_weights:

  A list of named numeric values, where the name of each element is one
  of the model outputs to be compared against the observed data, and the
  value is the weight for that output. Each weight can be a single
  number, or a pair of numbers. When the weight is a pair, the first
  number is the weight that will be used for underestimates (when the
  modeled value is smaller than the observed value), and the second is
  the weight for overestimates.

- data_definitions:

  A list of named string values, where the name of each element is one
  of the data columns in the `data_driver_pairs` and the value is that
  column's corresponding name in the model output. For all other data
  columns in the `data_driver_pairs`, it is assumed that the data column
  name matches a column in the model output.

- normalization_method:

  A string indicating the normalization method to be used when
  calculating the error metric; see below for more details.

- normalization_param:

  An (optional) parameter value used by some normalization methods. When
  `normalization_param` is `NULL`, a default value will be used, which
  depends on the particular normalization method. Otherwise, the
  user-specified value will be used. See details below.

- stdev_weight_method:

  A string indicating the method to be used when calculating the
  variance-based weights used in the error metric; see below for more
  details.

- stdev_weight_param:

  An (optional) parameter value used by some normalization methods. When
  `stdev_weight_param` is `NULL`, a default value will be used, which
  depends on the particular normalization method. Otherwise, the
  user-specified value will be used. See details below.

- regularization_method:

  A string indicating the regularization method to be used when
  calculating the regularization penalty term, or a function that
  calculates the penalty; see below for more details.

- dependent_arg_function:

  A function whose input argument is a named list of independent
  argument values, and which returns a named list of dependent argument
  values. If the `dependent_arg_function` is `NULL`, no dependent
  argument values will be calculated.

- post_process_function:

  A function whose input argument is a data frame representing the
  output from
  [`run_biocro`](https://rdrr.io/pkg/BioCro/man/run_biocro.html), and
  which returns a data frame, typically based on the input but with one
  or more new columns. If the `post_process_function` is `NULL`, no
  post-processing will be applied to the raw simulation output.

- extra_penalty_function:

  A function whose input argument is a data frame representing the
  output from
  [`run_biocro`](https://rdrr.io/pkg/BioCro/man/run_biocro.html), and
  which returns a numeric penalty to be added to the least-squares term
  when calculating the error metric. If the `extra_penalty_function` is
  `NULL`, no extra penalties will be added.

- verbose_startup:

  A logical (`TRUE` or `FALSE`) value indicating whether to print
  additional information to the R terminal when creating the objective
  function.

## Details

**Overview**

When parameterizing a BioCro model, the general idea is to vary a subset
of the model's `parameters` and/or `initial_values` to achieve the best
agreement with a set of observed data. The degree of agreement is
expressed as an "error metric," which includes terms derived from the
agreement between the modeled and observed values, as well as (optional)
penalty terms. A function that calculates an error metric when given a
set of argument values (`parameters` and/or `initial_values`) is called
an "objective function." Defining an objective function suitable for
BioCro parameterization can be complicated, but the `objective_function`
function helps to simplify the process of creating such a function. It
is designed to accommodate the following scenarios, which often occur in
the context of BioCro model parameterization:

- **Multi-year or multi-location data:** Often the model needs to be run
  several times with different drivers corresponding to multiple years
  or locations, and then the results from each individual simulation
  must be compared to associated sets of observed data. Here, this is
  handled through the `data_driver_pairs`, which allows the user to
  specify which drivers and data sets should be compared to each other.
  Optionally, it is also possible to specify different initial values or
  parameters for each set of drivers; for example, the atmospheric CO2
  concentration may need to change for different years, or soil
  properties may need to change for different locations.

- **Complicated normalization:** Care must be taken to ensure that
  certain years or output variables are not over-valued in the error
  metric; for example, one year may have more observations of leaf mass
  than another year, or the stem mass may be much larger than the leaf
  mass. Here, this is handled through pre-set normalization approaches,
  which can be specified through the `normalization_method` input. See
  below for more information.

- **Extra penalties:** Sometimes an optimizer chooses parameters that
  produce close agreement with the observed data, but are nevertheless
  not biologically resonable. For example, it may produce a sharp peak
  at a high leaf mass in between two measured points, when in reality,
  leaf mass should be nearly constant between them. In this case, it may
  be necessarily to add an extra penalty to the objective function that
  prevents the optimizer from choosing such values. Here, this is
  handled through the `extra_penalty_function` input.

- **Flexible weights:** Often a user would like to specify a weight for
  each variable being considered in the error metric, either to
  represent an uncertainty or to emphasize agreement with one output at
  the expense of another. For example, the seed mass may need a high
  weight to prioritize accurate yield predictions. Further, these
  weights may need to differ for underestimates as compared to
  overestimates; for example, measured root mass is often lower than the
  true root mass, so a user may wish to penalize underestimates of root
  mass more severely than overestimates. Here, these situations are
  handled through the `quantity_weights` and `stdev_weight_method`
  inputs.

- **Dependent parameters:** Sometimes one model parameter must be
  determined from one or more other parameters; for example, a user may
  wish to require that the leaf and stem maintenance respiration
  coefficients are identical. Here, this is handled through the
  `dependent_arg_function`, which allows the user to specify how the
  values of such "dependent" parameters should be determined from the
  values of the "independent" parameters.

- **Name mismatches:** Often a particular variable has different names
  in the data set and simulation output. Here, this is handled through
  the `data_definitions`, which allows the user to specify which columns
  in the model output should be compared to particular columns in the
  observed data.

- **Incomplete outputs:** Sometimes a model may not produce outputs that
  are directly comparable to the observed values; for example, a model
  may calculate seed and shell mass, while a data set includes pod mass,
  which is the sum of seed and shell. Here, this is handled by an
  optional `post_process_function`, which allows users to specify
  additional operations to perform on the model output; in this example,
  it would be used to calculate the pod mass so it can be compared to
  the observations.

**Error metric calculations**

As mentioned above, the overall error metric \\E\\ is calculated as

\$\$E = E\_{agreement} + P\_{user} + P\_{regularization},\$\$ where
\\E\_{agreement}\\ is determined by the agreement between the model and
observations, \\P\_{user}\\ is an optional user-defined penalty, and
\\P\_{regularization}\\ is an optional regularization penalty. These
terms are explained in more detail below:

- **Agreement term:** The agreement term \\E\_{agreement}\\ is
  calculated using a least-squares approach. In other words,

  \$\$E\_{agreement} = \sum_i \left(Y\_{obs}^i - Y\_{mod}^i \right)^2%
  \cdot \frac{w_i^{quantity} w_i^{data} w_i^{stdev}}{N_i},\$\$ where the
  sum runs over all \\n\\ observations; \\Y\_{obs}^i\\ and
  \\Y\_{mod}^i\\ are observed and modeled values of variable \\Y_i\\;
  \\w_i^{quantity}\\, \\w_i^{data}\\, and \\w_i^{stdev}\\ are weight
  factors that depend on the name of \\Y_i\\, the data set that includes
  the \\i^{th}\\ observation, and the standard deviation associated with
  \\Y\_{obs}^i\\, respectively; and \\N_i\\ is a normalization factor.

  Each value of \\Y\_{obs}^i\\ is specified at a particular time
  \\t_i\\. The corresponding modeled value, \\Y\_{mod}^i\\, is found by
  retrieving the value of the \\Y_i\\ variable at the closest time to
  \\t_i\\ that is included in the model output. It is assumed that the
  model always outputs the same sequence of time values each time it is
  run with a particular set of drivers, regardless of the provided
  argument values.

  The quantity-based weight factors \\w_i^{quantity}\\ are directly
  specified by the user via the `quantity_weights` input. For example,
  if `quantity_weights` has an element named `Leaf` that is equal to
  0.5, then \\w_i\\ will be equal to 0.5 whenever \\Y_i\\ represents a
  leaf mass value, regardless of which set of drivers or time point
  corresponds to \\Y_i\\. The weights can also be supplied as
  \\(w\_{under}, w\_{over})\\ pairs instead of single values; in this
  case, the value of \\w_i\\ depends on whether the model makes an
  underprediction or an overprediction: \\w_i = w\_{under}\\ when
  \\Y\_{mod}^i \< Y\_{obs}^i\\ and \\w_i = w\_{over}\\ otherwise.

  The data-set-based weight factors \\w_i^{data}\\ are directly
  specified by the user via the `weight` element of each data-driver
  pair. For example, if the second element of `data_driver_pairs` has a
  `weight` of 2.0, then \\w_i^{data}\\ will be equal to 2.0 for all
  observations from the corresponding data set.

  The standard-deviation-based weight factors \\w_i^{stdev}\\ are
  determined by the choice of `stdev_weight_method`; the available
  methods are discussed below.

  The normalization factors \\N_i\\ are determined by the choice of
  `normalization_method`; the available methods are discussed below.

  There are a few special cases where \\E\_{agreement}\\ is set to a
  very high value (`BioCroValidation:::FAILURE_VALUE`). This is done
  when a simulation fails to run, when the \\E\_{agreement}\\ term would
  otherwise evaluate to `NA`, or when the \\E\_{agreement}\\ term would
  otherwise evaluate to an infinite value.

- **User-defined penalty term:** The user-defined penalty term
  \\P\_{user}\\ is calculated by applying a function \\f\_{user}\\ to
  the full simulation output from each set of drivers. In other words,

  \$\$P\_{user} = \sum_k f\_{user} \left( M_k \right),\$\$ where the sum
  runs over all \\k\\ sets of drivers and \\M_k\\ is the model output
  when it is run with the \\k^{th}\\ set of drivers.

  The function \\f\_{user}\\ must accept a single data frame as an input
  and return a single numeric value as its output, but has no other
  requirements. It is specified via the `extra_penalty_function`. When
  `extra_penalty_function` is `NULL`, \\P\_{user}\\ is zero.

- **Regularization penalty term:** The regularization penalty term
  \\P\_{regularization}\\ is calculated from the values of the arguments
  being varied during the optimization by applying a function \\R\\. In
  other words,

  \$\$P\_{regularization} = R \left( X \right),\$\$ where \\X\\
  represents the model argument values.

  The function \\R\\ is determined by the choice of
  `regularization_method`; the available methods are discussed below.

**Standard-deviation-based weight methods**

The following pre-set methods are available for determining weight
factors from values of the standard deviation (\\\sigma\\), which can be
(optionally) supplied via the `data_stdev` elements of the
`data_driver_pairs`:

- `'equal'`: For this method, \\w_i^{stdev}\\ is always set to 1. In
  other words, all variances are treated as being equal, regardless of
  any user-supplied values. This is usually the best choice when values
  of \\\sigma\\ are unavailable or cannot be estimated.

- `'logarithm'`: For this method, \\w_i^{stdev}\\ is calculated as

  \$\$w_i^{stdev} =% ln \left( \frac{1}{\sigma_i + \epsilon}
  \right),\$\$ where \\ln\\ denotes a logarithm with base \\e\\ and
  \\\epsilon\\ is a small number included to prevent numerical errors
  that would otherwise occur when \\\sigma_i = 0\\. This method was used
  in the original Soybean-BioCro paper.

  The value of \\\epsilon\\ is specified by `stdev_weight_param`, which
  defaults to `1e-5` if `stdev_weight_param` is `NULL` when using this
  method. With the default value of \\\epsilon\\, \\w_i^{stdev} \approx
  11.512\\ when \\\sigma = 0\\.

  Note: this method should be used with caution, because \\w_i^{stdev}\\
  is zero for \\\sigma_i = 1 - \epsilon\\, and it becomes negative for
  larger values of \\\sigma_i\\.

- `'inverse'`: For this method, \\w_i^{stdev}\\ is calculated as

  \$\$w_i^{stdev} = \frac{1}{\sigma_i^2 + \epsilon},\$\$ where
  \\\epsilon\\ is a small number included to prevent numerical errors
  that would otherwise occur when \\\sigma_i = 0\\.

  The value of \\\epsilon\\ is specified by `stdev_weight_param`, which
  defaults to `1e-1` if `stdev_weight_param` is `NULL` when using this
  method. With the default value of \\\epsilon\\, \\w_i^{stdev} = 10\\
  when \\\sigma_i = 0\\.

If any values of \\w_i^{stdev}\\ are undefined, negative, or infinite,
an error message will occur (see the "Input checks" section below).

**Normalization methods**

The following pre-set normalization methods are available:

- `'equal'`: For this method, \\N_i\\ is always set to 1. In other
  words, no normalization is performed.

- `'mean'`: For this method, when \\Y_i\\ is named `vtype` and the
  observation is from a set called `vdata`, then

  \$\$N_i = n\_{vtype}^{vdata} \cdot n\_{data},\$\$ where
  \\n\_{vtype}^{vdata}\\ is the number of observations of type `vtype`
  that are included in `vdata` and \\n\_{data}\\ is the total number of
  data-driver pairs. In this case, the error term \\E\_{agreement}\\
  becomes a mean error across the full set of drivers, hence the name
  for this method. This approach avoids over-representing drivers with
  larger numbers of associated observations when determining
  \\E\_{agreement}\\. It also preserves the overall magnitude of
  \\E\_{agreement}\\ when data-driver pairs are added.

- `'max'`: For this method, when \\Y_i\\ is named `vtype` and the
  observation is from a set called `vdata`, then

  \$\$N_i = \left( max\_{vtype}^{vdata} \right)^2 + \epsilon,\$\$ where
  \\max\_{vtype}^{vdata}\\ is the maximum observed absolute value of
  `vtype` across `vdata` and \\\epsilon\\ is a small number included to
  prevent numerical errors that would otherwise occur when
  \\max\_{vtype}^{vdata} = 0\\. In this case, the observed and modeled
  values that appear in the equation for \\E\_{agreement}\\ are
  essentially normalized by their maximum magnitude, hence the name for
  this method. This approach avoids over-representing variable types
  with larger magnitude when determining \\E\_{agreement}\\.

  The value of \\\epsilon\\ is specified by `normalization_param`, which
  defaults to `1e-1` if `normalization_param` is `NULL` when using this
  method. With the default value of \\\epsilon\\, \\N_i = 10\\ when
  \\max\_{vtype}^{vdata} = 0\\.

- `'obs'`: For this method,

  \$\$N_i = \left( Y\_{obs}^i \right)^2 + \epsilon,\$\$ where
  \\\epsilon\\ is a small number included to prevent numerical errors
  that would otherwise occur when \\Y\_{obs}^i = 0\\. In this case, the
  equation for \\E\_{agreement}\\ essentially uses relative differences
  rather than absolute differences, which is achieved by normalizing by
  the observed values, hence the name. This approach avoids
  over-representing time points where a particular quantity takes its
  largest values when determining \\E\_{agreement}\\.

  The value of \\\epsilon\\ is specified by `normalization_param`, which
  defaults to `1e-1` if `normalization_param` is `NULL` when using this
  method. With the default value of \\\epsilon\\, \\N_i = 10\\ when
  \\Y\_{obs}^i = 0\\.

- `'mean_max'`: For this method, the "mean" and "max" methods are
  combined so that

  \$\$N_i = n\_{vtype}^{vdrivers} \cdot n\_{data}% \cdot \left\[ \left(
  max\_{vtype}^{vdata} \right)^2 + \epsilon \right\].\$\$ This approach
  avoids over-representing drivers with larger numbers of associated
  observations, and variable types with larger magnitudes. This method
  is used for parameterizing Soybean-BioCro.

  The value \\\epsilon\\ is specified by `normalization_param`, using
  the same default value as in the `'max'` method.

- `'mean_obs'`: For this method, the "mean" and "obs" methods are
  combined so that

  \$\$N_i = n\_{vtype}^{vdrivers} \cdot n\_{data}% \cdot \left\[ \left(
  Y\_{obs}^i \right)^2 + \epsilon \right\].\$\$ This approach avoids
  over-representing drivers with larger numbers of associated
  observations, and time points with larger observed values.

  The value \\\epsilon\\ is specified by `normalization_param`, using
  the same default value as in the `'obs'` method.

In most situations, it is recommended to use either `'mean_max'` or
`'mean_obs'` depending on user preference or performance.

**Regularization methods**

The following pre-set regularization methods are available:

- `'none'`: For this method, \\P\_{regularization}\\ is always set to 0.
  In other words, no regularization is performed.

- `'L1'` or `'lasso'`: For this method, \\P\_{regularization}\\ is given
  by the sum of the absolute values of each independent argument,
  multiplied by a "regularization parameter" \\\lambda\\ that sets the
  overall weight of the penalty:

  \$\$P\_{regularization} = \lambda \sum_j \| X_j \|,\$\$ where the sum
  runs over all \\j\\ independent arguments, and \\X_j\\ is the value of
  the \\j^{th}\\ argument. See the "Value" section below for details of
  how to specify \\\lambda\\.

- `'L2'` or `'ridge'`: For this method, \\P\_{regularization}\\ is given
  by the sum of the squared values of each independent argument,
  multiplied by a "regularization parameter" \\\lambda\\ that sets the
  overall weight of the penalty:

  \$\$P\_{regularization} = \lambda \sum_j X_j^2,\$\$ where the sum runs
  over all \\j\\ independent arguments, and \\X_j\\ is the value of the
  \\j^{th}\\ argument. See the "Value" section below for details of how
  to specify \\\lambda\\.

It is also possible to supply a function that accepts two input
arguments (`x` and `lambda`, as described in the "Value" section below)
and returns a numeric penalty value.

**Input checks**

Several checks are made to ensure that the objective function is
properly defined. These checks include, but are not limited to, the
following:

- Ensuring that each set of drivers in `data_driver_pairs` defines a
  valid dynamical system along with the `base_model_definition`. This is
  accomplished using
  [`validate_dynamical_system_inputs`](https://rdrr.io/pkg/BioCro/man/dynamical_system.html).

- Ensuring that the model output corresponding to each set of drivers
  spans the times at which the observations were made.

- Ensuring that each variable type in the data elements of
  `data_driver_pairs` matches a corresponding column in the model
  output, when accounting for the `data_definitions` and
  `post_process_function`.

- Ensuring that each independent and dependent argument name is a member
  of the model's `parameters` or `initial_values`. Internally, argument
  names are passed to
  [`partial_run_biocro`](https://rdrr.io/pkg/BioCro/man/partial_application.html)
  via its `arg_names` input. Note that argument names passed to
  `partial_run_biocro` can technically include members of the `drivers`,
  but it is unlikely that the value of a driver would be varied during
  an optimization, so the argument names are not allowed to include
  columns in the `drivers`.

- Ensuring that the optional `dependent_arg_function`,
  `post_process_function`, and `extra_penalty_function` functions can be
  run without causing errors.

- Ensuring that certain values are finite (such as \\Y\_{obs}\\,
  \\\sigma_i\\, \\w_i^{stdev}\\, and \\N_i\\), and that certain values
  are not negative (such as \\\sigma_i\\, \\w_i^{stdev}\\, and \\N_i\\).

If any issues are detected, an informative error message will be sent.
Note that several of these checks require running the model with each
set of drivers. For these checks, the argument values specified by
`independent_args` will be used, so they should be valid or otherwise
reasonable values.

If an error message occurs when `verbose_startup` was set to `FALSE`, it
is recommended to call this function again with `verbose_startup` set to
`TRUE`, since the additional output can be helpful for troubleshooting.

## Value

A function `obj_fun` with signature
`obj_fun(x, lambda = 0, return_terms = FALSE, debug_mode = 'minimal')`.

Here, `x` is a numeric vector of values of the independent arguments (in
the same order as in `independent_arg_names`), and `lambda` is the value
of the regularization parameter.

When `return_terms` is `FALSE`, `obj_fun` returns values of the error
metric \\E\\. When `return_terms` is `TRUE`, `obj_fun` returns a list
including each individual term of the total error metric.

During optimization, `return_terms` should always be `FALSE`. Setting it
to `TRUE` can be useful for troubleshooting, or for diagnostics such as
checking the quality of fit for each of the data-driver pairs.

When `debug_mode` is `'minimal'` or `'everything'`, `obj_fun` operates
in "debug mode." In the "minimal" version, problematic values of `x`
will be printed to the terminal, along with an explanation of the
problem that was caused. In the "everything" version, the values of `x`
and the error metric will be printed to the terminal every time
`obj_fun` is called. When `debug_mode` is `'none'` (or any other value
not mentioned above), no debug printing will occur.

Debug mode can be useful when troubleshooting a problem with an
optimization, since it provides a record of any problematic parameter
combinations. Once a problematic set of argument values is identified,
it can be investigated further by calling `obj_fun` again with `x` set
to the problematic values and `return_terms` set to `TRUE`.

When setting `debug_mode` to `TRUE`, also consider using
[`sink`](https://rdrr.io/r/base/sink.html) to write the outputs to a
file instead of the R terminal. In that case, there will still be a
record even if R crashes.

## Examples

``` r
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
  # alter them. The tables also include other columns we do not wish to use in
  # this example. So, we will define a short helper function that can be used to
  # pre-process each table.
  process_table <- function(x) {
    within(x, {
      # Define new `time` column
      time = (DOY - 1) * 24.0

      # Remove unneeded columns
      DOY                 = NULL
      Seed_Mg_per_ha      = NULL
      Litter_Mg_per_ha    = NULL
      CumLitter_Mg_per_ha = NULL
    })
  }

  # The data-driver pairs can now be created by associating each data set with
  # its corresponding weather data. Here we will weight the 2005 data twice as
  # heavily as the 2002 data. Note that we also specify different atmospheric
  # CO2 concentrations for each year.
  data_driver_pairs <- list(
    ambient_2002 = list(
      data       = process_table(soyface_biomass[['ambient_2002']]),
      data_stdev = process_table(soyface_biomass[['ambient_2002_std']]),
      drivers    = BioCro::soybean_weather[['2002']],
      parameters = list(Catm = with(BioCro::catm_data, {Catm[year == '2002']})),
      weight     = 1
    ),
    ambient_2005 = list(
      data       = process_table(soyface_biomass[['ambient_2005']]),
      data_stdev = process_table(soyface_biomass[['ambient_2005_std']]),
      drivers    = BioCro::soybean_weather[['2005']],
      parameters = list(Catm = with(BioCro::catm_data, {Catm[year == '2005']})),
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
    Rep_Mg_per_ha = 'Pod'
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

  initial_guess <- as.numeric(independent_args)

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

  # We want to use the regularization term to penalize deviations away from the
  # initial guess, so we will define a custom L2 regularization function
  regularization_function <- function(x, lambda) {
    lambda * sum((x - initial_guess)^2)
  }

  # Now we can finally create the objective function
  obj_fun <- objective_function(
    base_model_definition,
    data_driver_pairs,
    independent_args,
    quantity_weights,
    data_definitions = data_definitions,
    stdev_weight_method = 'logarithm',
    regularization_method = regularization_function,
    dependent_arg_function = dependent_arg_function,
    post_process_function = post_process_function,
    extra_penalty_function = extra_penalty_function
  )

  # This function could now be passed to an optimizer; here we will simply
  # evaluate it for two sets of parameter values.

  # Try doubling each parameter value and setting lambda to a nonzero value; in
  # this case, the value of the objective function increases, indicating a lower
  # degree of agreement between the model and the observed data. Here we will
  # call `obj_fun` in debug mode, which will automatically print the value of
  # the error metric.
  cat('\nError metric calculated by doubling the original argument values:\n')
  error_metric <- obj_fun(2 * initial_guess, 0.001, debug_mode = 'everything')

  # We can also see the values of each term that makes up the error metric;
  # again, we will call `obj_fun` in debug mode for automatic printing.
  cat('\nError metric terms calculated by doubling the original argument values:\n')
  error_terms <-
    obj_fun(2 * initial_guess, 0.001, return_terms = TRUE, debug_mode = 'everything')
}
#> Loading required package: BioCro
#> 
#> Driver-specific initial values:
#> 
#>   None
#> 
#> Driver-specific parameters:
#> 
#> List of 2
#>  $ ambient_2002:List of 1
#>   ..$ Catm: num 373
#>  $ ambient_2005:List of 1
#>   ..$ Catm: num 379
#> 
#> The independent arguments and their initial values:
#> 
#> List of 2
#>  $ alphaLeaf: num 23.4
#>  $ betaLeaf : num -18.1
#> 
#> The dependent arguments and their initial values:
#> 
#> List of 1
#>  $ alphaStem: num 23.4
#> 
#> The full data definitions:
#> 
#> List of 3
#>  $ Leaf_Mg_per_ha: chr "Leaf"
#>  $ Stem_Mg_per_ha: chr "Stem"
#>  $ Rep_Mg_per_ha : chr "Pod"
#> 
#> Normalization method: MEAN_MAX with eps = 0.1 
#> 
#> Standard-deviation-based weight method: LOGARITHM with eps = 1e-05 
#> 
#> The user-supplied data in long form:
#> 
#> $ambient_2002
#>    time quantity_name quantity_value quantity_stdev time_index expected_npts
#> 1  4272          Leaf   0.1802843394   0.0408155501        649          3288
#> 2  4512          Leaf   0.5544619422   0.1638632739        889          3288
#> 3  4848          Leaf   1.3265529308   0.1337744335       1225          3288
#> 4  5184          Leaf   1.6979440069   0.2283266576       1561          3288
#> 5  5520          Leaf   1.8077427820   0.2024754215       1897          3288
#> 6  5880          Leaf   1.5788136482   0.0754751654       2257          3288
#> 7  6192          Leaf   0.9475377733   0.3445500325       2569          3288
#> 8  6888          Leaf   0.0000000000   0.0000000000       3265          3288
#> 9  4272          Stem   0.0852449694   0.0170797372        649          3288
#> 10 4512          Stem   0.4188538932   0.1384490248        889          3288
#> 11 4848          Stem   1.7110673664   0.1837107594       1225          3288
#> 12 5184          Stem   2.8928258965   0.4487440652       1561          3288
#> 13 5520          Stem   3.6859142604   0.4534474707       1897          3288
#> 14 5880          Stem   3.7452607171   0.2753213561       2257          3288
#> 15 6192          Stem   3.6184015745   0.1510453777       2569          3288
#> 16 6888          Stem   2.3057012247   0.1483892609       3265          3288
#> 17 4272           Pod   0.0000000000   0.0000000000        649          3288
#> 18 4512           Pod   0.0000000000   0.0000000000        889          3288
#> 19 4848           Pod   0.0003171479   0.0005493162       1225          3288
#> 20 5184           Pod   0.0793963255   0.0309899985       1561          3288
#> 21 5520           Pod   1.5545713035   0.2184025435       1897          3288
#> 22 5880           Pod   3.9760135605   0.5735488281       2257          3288
#> 23 6192           Pod   6.5446506556   0.7434407011       2569          3288
#> 24 6888           Pod   7.0089676285   0.1418286169       3265          3288
#>         norm      w_var
#> 1   53.88694  3.1984472
#> 2   53.88694  1.8086619
#> 3   53.88694  2.0115255
#> 4   53.88694  1.4769342
#> 5   53.88694  1.5970874
#> 6   53.88694  2.5838191
#> 7   53.88694  1.0654869
#> 8   53.88694 11.5129255
#> 9  226.03165  4.0692772
#> 10 226.03165  1.9771808
#> 11 226.03165  1.6943383
#> 12 226.03165  0.8012803
#> 13 226.03165  0.7908538
#> 14 226.03165  1.2897800
#> 15 226.03165  1.8901088
#> 16 226.03165  1.9078489
#> 17 787.61004 11.5129255
#> 18 787.61004 11.5129255
#> 19 787.61004  7.4887956
#> 20 787.61004  3.4737681
#> 21 787.61004  1.5213696
#> 22 787.61004  0.5558948
#> 23 787.61004  0.2964528
#> 24 787.61004  1.9530654
#> 
#> $ambient_2005
#>    time quantity_name quantity_value quantity_stdev time_index expected_npts
#> 1  4104          Leaf     0.22227188     0.03289659        577          2952
#> 2  4440          Leaf     0.84603750     0.14679830        913          2952
#> 3  4776          Leaf     1.18446563     0.33807429       1249          2952
#> 4  5112          Leaf     2.21805938     0.15217591       1585          2952
#> 5  5448          Leaf     2.14744687     0.11907759       1921          2952
#> 6  5784          Leaf     1.51948125     0.51280870       2257          2952
#> 7  6120          Leaf     0.06575625     0.06168624       2593          2952
#> 8  6456          Leaf     0.00000000     0.00000000       2929          2952
#> 9  4104          Stem     0.18880312     0.01431814        577          2952
#> 10 4440          Stem     0.85220625     0.19883006        913          2952
#> 11 4776          Stem     1.61896875     0.60528625       1249          2952
#> 12 5112          Stem     4.04361563     0.55987405       1585          2952
#> 13 5448          Stem     4.47772500     0.30674464       1921          2952
#> 14 5784          Stem     3.89208750     0.37910849       2257          2952
#> 15 6120          Stem     2.89905000     0.22082398       2593          2952
#> 16 6456          Stem     2.17560000     0.24325473       2929          2952
#> 17 4104           Pod     0.00000000     0.00000000        577          2952
#> 18 4440           Pod     0.00000000     0.00000000        913          2952
#> 19 4776           Pod     0.00000000     0.00000000       1249          2952
#> 20 5112           Pod     0.29925000     0.16427520       1585          2952
#> 21 5448           Pod     2.30455312     0.43414807       1921          2952
#> 22 5784           Pod     5.53277813     0.58847698       2257          2952
#> 23 6120           Pod     5.37107813     0.52004438       2593          2952
#> 24 6456           Pod     6.37225313     0.63309086       2929          2952
#>        norm      w_var
#> 1   80.3166  3.4140824
#> 2   80.3166  1.9186276
#> 3   80.3166  1.0844600
#> 4   80.3166  1.8826524
#> 5   80.3166  2.1278960
#> 6   80.3166  0.6678329
#> 7   80.3166  2.7855322
#> 8   80.3166 11.5129255
#> 9  322.4003  4.2455301
#> 10 322.4003  1.6152545
#> 11 322.4003  0.5020373
#> 12 322.4003  0.5800256
#> 13 322.4003  1.1817071
#> 14 322.4003  0.9699065
#> 15 322.4003  1.5103441
#> 16 322.4003  1.4136050
#> 17 651.2898 11.5129255
#> 18 651.2898 11.5129255
#> 19 651.2898 11.5129255
#> 20 651.2898  1.8061514
#> 21 651.2898  0.8343466
#> 22 651.2898  0.5302005
#> 23 651.2898  0.6538219
#> 24 651.2898  0.4571255
#> 
#> The user-supplied quantity weights:
#> 
#> List of 3
#>  $ Leaf: num [1:2] 0.5 0.5
#>  $ Stem: num [1:2] 0.5 0.25
#>  $ Pod : num [1:2] 1 1
#> 
#> The user-supplied data-driver pair weights:
#> 
#> List of 2
#>  $ ambient_2002: num 1
#>  $ ambient_2005: num 2
#> 
#> Regularization method: user-supplied function:
#> 
#> function (x, lambda) 
#> {
#>     lambda * sum((x - initial_guess)^2)
#> }
#> <environment: 0x555f3cc2dc38>
#> 
#> Dependent argument function: user-supplied function:
#> 
#> function (ind_args) 
#> {
#>     list(alphaStem = ind_args[["alphaLeaf"]])
#> }
#> <environment: 0x555f3cc2dc38>
#> 
#> Post-processing function: user-supplied function:
#> 
#> function (sim_res) 
#> {
#>     within(sim_res, {
#>         Pod = Grain + Shell
#>     })
#> }
#> <environment: 0x555f3cc2dc38>
#> 
#> Extra penalty function: user-supplied function:
#> 
#> function (sim_res) 
#> {
#>     max_leaf <- max(sim_res[["Leaf"]], na.rm = TRUE)
#>     if (is.na(max_leaf) || max_leaf > 4) {
#>         1e+05
#>     }
#>     else {
#>         0
#>     }
#> }
#> <environment: 0x555f3cc2dc38>
#> 
#> The initial error metric terms:
#> 
#> List of 2
#>  $ terms_from_data_driver_pairs:List of 2
#>   ..$ ambient_2002:List of 2
#>   .. ..$ least_squares_terms:List of 3
#>   .. .. ..$ Leaf: num 0.0908
#>   .. .. ..$ Pod : num 0.0182
#>   .. .. ..$ Stem: num 0.00474
#>   .. ..$ extra_penalty      : num 0
#>   ..$ ambient_2005:List of 2
#>   .. ..$ least_squares_terms:List of 3
#>   .. .. ..$ Leaf: num 0.124
#>   .. .. ..$ Pod : num 0.0261
#>   .. .. ..$ Stem: num 0.0376
#>   .. ..$ extra_penalty      : num 0
#>  $ regularization_penalty      : num 0
#> 
#> The initial error metric value:
#> 
#> [1] 0.3012756
#> 
#> Error metric calculated by doubling the original argument values:
#> 
#> Time: 2026-03-11 23:03:09.565925      Independent argument values : 46.73542898017520030862215207889676, -36.22026162863519971324421931058168
#> 
#> Time: 2026-03-11 23:03:09.972865      Error metric : 2.99754729469563185872971189382952
#> 
#> Error metric terms calculated by doubling the original argument values:
#> 
#> Time: 2026-03-11 23:03:09.973026      Independent argument values : 46.73542898017520030862215207889676, -36.22026162863519971324421931058168
#> 
#> Time: 2026-03-11 23:03:10.367931      Error metric terms : 
#> 
#> List of 2
#>  $ terms_from_data_driver_pairs:List of 2
#>   ..$ ambient_2002:List of 2
#>   .. ..$ least_squares_terms:List of 3
#>   .. .. ..$ Leaf: num 0.122
#>   .. .. ..$ Pod : num 0.154
#>   .. .. ..$ Stem: num 0.39
#>   .. ..$ extra_penalty      : num 0
#>   ..$ ambient_2005:List of 2
#>   .. ..$ least_squares_terms:List of 3
#>   .. .. ..$ Leaf: num 0.139
#>   .. .. ..$ Pod : num 0.179
#>   .. .. ..$ Stem: num 1.14
#>   .. ..$ extra_penalty      : num 0
#>  $ regularization_penalty      : num 0.874
#> 
```
