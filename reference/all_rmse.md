# Calculate RMSE values

Given a base model definition, drivers to run the model, and observed
values of model outputs, `all_rmse` runs the model and determines the
root mean squared error (RMSE) value for each driver and observed
quantity.

This function is typically used after minimizing the value of an
objective function created by
[`objective_function`](https://biocro.github.io/BioCroValidation/reference/objective_function.md).
The RMSE values can be calculated for the original drivers used for
training, and potentially other sets of validation drivers. This enables
quantification of model performance.

## Usage

``` r
all_rmse(
    base_model_definition,
    data_driver_pairs,
    data_definitions = list(),
    post_process_function = NULL,
    verbose_startup = TRUE
  )
```

## Arguments

- base_model_definition:

  See
  [`objective_function`](https://biocro.github.io/BioCroValidation/reference/objective_function.md).

- data_driver_pairs:

  See
  [`objective_function`](https://biocro.github.io/BioCroValidation/reference/objective_function.md).
  Note that any weights will be reset to 1 automatically.

- data_definitions:

  See
  [`objective_function`](https://biocro.github.io/BioCroValidation/reference/objective_function.md).

- post_process_function:

  See
  [`objective_function`](https://biocro.github.io/BioCroValidation/reference/objective_function.md).

- verbose_startup:

  A logical (`TRUE` or `FALSE`) value indicating whether to print
  additional information to the R terminal when evaluating the RMSE
  terms.

## Details

This function is essentially a wrapper for
[`objective_function`](https://biocro.github.io/BioCroValidation/reference/objective_function.md)
where all weights are set to 1 and the `mean` normalization method is
used. The objective function is evaluated, and each of the resulting
squared error terms can be used to calculate the RMSE.

## Value

A data frame where each row corresponds to one of the data-driver pairs.
The `data_driver_pair` column indicates the name of the data-driver
pair, and all other columns contain RMSE values for each quantity type.

## Examples

``` r
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
#> List of 1
#>  $ alpha1: num 0
#> 
#> The dependent arguments and their initial values:
#> 
#>   None
#> 
#> The full data definitions:
#> 
#> List of 3
#>  $ Leaf_Mg_per_ha: chr "Leaf"
#>  $ Stem_Mg_per_ha: chr "Stem"
#>  $ Rep_Mg_per_ha : chr "Pod"
#> 
#> Normalization method: MEAN 
#> 
#> Standard-deviation-based weight method: EQUAL 
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
#> 23 4272          Stem   0.0852449694   0.0170797372        649          3288
#> 24 4512          Stem   0.4188538932   0.1384490248        889          3288
#> 25 4848          Stem   1.7110673664   0.1837107594       1225          3288
#> 26 5184          Stem   2.8928258965   0.4487440652       1561          3288
#> 27 5520          Stem   3.6859142604   0.4534474707       1897          3288
#> 28 5880          Stem   3.7452607171   0.2753213561       2257          3288
#> 29 6192          Stem   3.6184015745   0.1510453777       2569          3288
#> 30 6888          Stem   2.3057012247   0.1483892609       3265          3288
#> 45 4272           Pod   0.0000000000   0.0000000000        649          3288
#> 46 4512           Pod   0.0000000000   0.0000000000        889          3288
#> 47 4848           Pod   0.0003171479   0.0005493162       1225          3288
#> 48 5184           Pod   0.0793963255   0.0309899985       1561          3288
#> 49 5520           Pod   1.5545713035   0.2184025435       1897          3288
#> 50 5880           Pod   3.9760135605   0.5735488281       2257          3288
#> 51 6192           Pod   6.5446506556   0.7434407011       2569          3288
#> 52 6888           Pod   7.0089676285   0.1418286169       3265          3288
#>    norm w_var
#> 1    16     1
#> 2    16     1
#> 3    16     1
#> 4    16     1
#> 5    16     1
#> 6    16     1
#> 7    16     1
#> 8    16     1
#> 23   16     1
#> 24   16     1
#> 25   16     1
#> 26   16     1
#> 27   16     1
#> 28   16     1
#> 29   16     1
#> 30   16     1
#> 45   16     1
#> 46   16     1
#> 47   16     1
#> 48   16     1
#> 49   16     1
#> 50   16     1
#> 51   16     1
#> 52   16     1
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
#> 25 4104          Stem     0.18880312     0.01431814        577          2952
#> 26 4440          Stem     0.85220625     0.19883006        913          2952
#> 27 4776          Stem     1.61896875     0.60528625       1249          2952
#> 28 5112          Stem     4.04361563     0.55987405       1585          2952
#> 29 5448          Stem     4.47772500     0.30674464       1921          2952
#> 30 5784          Stem     3.89208750     0.37910849       2257          2952
#> 31 6120          Stem     2.89905000     0.22082398       2593          2952
#> 32 6456          Stem     2.17560000     0.24325473       2929          2952
#> 49 4104           Pod     0.00000000     0.00000000        577          2952
#> 50 4440           Pod     0.00000000     0.00000000        913          2952
#> 51 4776           Pod     0.00000000     0.00000000       1249          2952
#> 52 5112           Pod     0.29925000     0.16427520       1585          2952
#> 53 5448           Pod     2.30455312     0.43414807       1921          2952
#> 54 5784           Pod     5.53277813     0.58847698       2257          2952
#> 55 6120           Pod     5.37107813     0.52004438       2593          2952
#> 56 6456           Pod     6.37225313     0.63309086       2929          2952
#>    norm w_var
#> 1    16     1
#> 2    16     1
#> 3    16     1
#> 4    16     1
#> 5    16     1
#> 6    16     1
#> 7    16     1
#> 8    16     1
#> 25   16     1
#> 26   16     1
#> 27   16     1
#> 28   16     1
#> 29   16     1
#> 30   16     1
#> 31   16     1
#> 32   16     1
#> 49   16     1
#> 50   16     1
#> 51   16     1
#> 52   16     1
#> 53   16     1
#> 54   16     1
#> 55   16     1
#> 56   16     1
#> 
#> The user-supplied quantity weights:
#> 
#> List of 3
#>  $ Leaf: num [1:2] 1 1
#>  $ Stem: num [1:2] 1 1
#>  $ Pod : num [1:2] 1 1
#> 
#> The user-supplied data-driver pair weights:
#> 
#> List of 2
#>  $ ambient_2002: num 1
#>  $ ambient_2005: num 1
#> 
#> Regularization method: NONE
#> 
#> Dependent argument function: none
#> 
#> Post-processing function: user-supplied function:
#> 
#> function (sim_res) 
#> {
#>     within(sim_res, {
#>         Pod = Grain + Shell
#>     })
#> }
#> <environment: 0x55630d470378>
#> 
#> Extra penalty function: none
#> 
#> The initial error metric terms:
#> 
#> List of 2
#>  $ terms_from_data_driver_pairs:List of 2
#>   ..$ ambient_2002:List of 2
#>   .. ..$ least_squares_terms:List of 3
#>   .. .. ..$ Leaf: num 0.0576
#>   .. .. ..$ Pod : num 0.35
#>   .. .. ..$ Stem: num 0.1
#>   .. ..$ extra_penalty      : num 0
#>   ..$ ambient_2005:List of 2
#>   .. ..$ least_squares_terms:List of 3
#>   .. .. ..$ Leaf: num 0.0215
#>   .. .. ..$ Pod : num 0.861
#>   .. .. ..$ Stem: num 0.0842
#>   .. ..$ extra_penalty      : num 0
#>  $ regularization_penalty      : num 0
#> 
#> The initial error metric value:
#> 
#> [1] 1.474685
#>        Leaf       Pod      Stem data_driver_pair
#> 1 0.3394553 0.8370377 0.4482304     ambient_2002
#> 2 0.2072064 1.3119349 0.4104760     ambient_2005
```
