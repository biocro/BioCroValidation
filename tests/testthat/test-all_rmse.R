# Specify key inputs to use for these tests
model <- BioCro::soybean
model$ode_solver <- BioCro::default_ode_solvers[['homemade_euler']]

process_table <- function(x) {
    within(x, {
        time                = (DOY - 1) * 24.0
        DOY                 = NULL
        Seed_Mg_per_ha      = NULL
        Litter_Mg_per_ha    = NULL
        CumLitter_Mg_per_ha = NULL
    })
}

ddps <- list(
        ambient_2002 = list(
        data       = process_table(soyface_biomass[['ambient_2002']]),
        data_stdev = process_table(soyface_biomass[['ambient_2002_std']]),
        drivers    = BioCro::soybean_weather[['2002']],
        weight     = 1
    ),
        ambient_2005 = list(
        data       = process_table(soyface_biomass[['ambient_2005']]),
        drivers    = BioCro::soybean_weather[['2005']],
        weight     = 2
    )
)

data_definitions <- list(
    Leaf_Mg_per_ha = 'Leaf',
    Stem_Mg_per_ha = 'Stem',
    Rep_Mg_per_ha = 'Pod',
    LAI = 'lai'
)

post_process_function <- function(x) {
    within(x, {Pod = Grain + Shell})
}

verbose_startup <- FALSE

test_that('RMSE values can be calculated', {
    rmse <- expect_silent(
        all_rmse(
            model,
            ddps,
            data_definitions,
            post_process_function,
            verbose_startup
        )
    )
})
