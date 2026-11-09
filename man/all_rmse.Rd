\name{all_rmse}

\alias{all_rmse}

\title{Calculate RMSE values}

\description{
  Given a base model definition, drivers to run the model, and observed values
  of model outputs, \code{all_rmse} runs the model and determines the root mean
  squared error (RMSE) value for each driver and observed quantity.

  This function is typically used after minimizing the value of an objective
  function created by \code{\link{objective_function}}. The RMSE values can be
  calculated for the original drivers used for training, and potentially other
  sets of validation drivers. This enables quantification of model performance.
}

\usage{
  all_rmse(
    base_model_definition,
    data_driver_pairs,
    data_definitions = list(),
    post_process_function = NULL,
    verbose_startup = TRUE
  )
}

\arguments{
  \item{base_model_definition}{
    See \code{\link{objective_function}}.
  }

  \item{data_driver_pairs}{
    See \code{\link{objective_function}}. Note that any weights will be reset to
    1 automatically.
  }

  \item{data_definitions}{
    See \code{\link{objective_function}}.
  }

  \item{post_process_function}{
    See \code{\link{objective_function}}.
  }

  \item{verbose_startup}{
    A logical (\code{TRUE} or \code{FALSE}) value indicating whether to print
    additional information to the R terminal when evaluating the RMSE terms.
  }
}

\details{
  This function is essentially a wrapper for \code{\link{objective_function}}
  where all weights are set to 1 and the \code{mean} normalization method is
  used. The objective function is evaluated, and each of the resulting squared
  error terms can be used to calculate the RMSE.
}

\value{
  A data frame where each row corresponds to one of the data-driver pairs. The
  \code{data_driver_pair} column indicates the name of the data-driver pair,
  and all other columns contain RMSE values for each quantity type.
}

\examples{
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
      LAI                 = NULL
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

  # Now we can calculate RMSE terms
  rmse_table <- all_rmse(
    base_model_definition,
    data_driver_pairs,
    data_definitions,
    post_process_function
  )

  print(rmse_table)
}
}
