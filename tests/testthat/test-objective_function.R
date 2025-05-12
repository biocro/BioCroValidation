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

independent_args <- with(BioCro::soybean[['parameters']], {
    list(alphaLeaf = alphaLeaf, betaLeaf = betaLeaf)
})

data_definitions <- list(
    Leaf_Mg_per_ha = 'Leaf',
    Stem_Mg_per_ha = 'Stem',
    Rep_Mg_per_ha = 'Pod'
)

dependent_arg_function <- function(x) {
    list(alphaStem = x[['alphaLeaf']])
}

post_process_function <- function(x) {
    within(x, {Pod = Grain + Shell})
}

quantity_weights <- list(
    Leaf = 0.5,
    Stem = c(0.5, 0.25),
    Pod = 1
)

verbose_startup <- FALSE

# Run tests
test_that('Objective functions can be created and behave as expected', {
    # Two data-driver pairs, no dependent arguments
    obj_fun <- expect_silent(
        objective_function(
            model,
            ddps,
            independent_args,
            quantity_weights,
            data_definitions = data_definitions,
            post_process_function = post_process_function,
            verbose_startup = verbose_startup
        )
    )

    expect_silent(
        obj_fun(as.numeric(independent_args))
    )

    # One data-driver pair, no dependent arguments
    obj_fun <- expect_silent(
        objective_function(
            model,
            ddps[1],
            independent_args,
            quantity_weights,
            data_definitions = data_definitions,
            post_process_function = post_process_function,
            verbose_startup = verbose_startup
        )
    )

    # Two data-driver pairs, with dependent arguments and L2 regularization
    obj_fun <- expect_silent(
        objective_function(
            model,
            ddps,
            independent_args,
            quantity_weights,
            data_definitions = data_definitions,
            dependent_arg_function = dependent_arg_function,
            post_process_function = post_process_function,
            regularization_method = 'L2',
            verbose_startup = verbose_startup
        )
    )

    expect_silent(
        obj_fun(as.numeric(independent_args), lambda = 0.5)
    )

    expect_true(
        is.list(
            obj_fun(as.numeric(independent_args), lambda = 0.5, return_terms = TRUE)
        )
    )
})

test_that('Bad definitions are detected', {
    expect_error(
        objective_function(
            model,
            within(ddps, {ambient_2005$drivers$temp = NULL}),
            independent_args,
            quantity_weights,
            data_definitions = data_definitions,
            post_process_function = post_process_function,
            verbose_startup = verbose_startup
        ),
        'The following drivers did not form a valid dynamical system: ambient_2005'
    )
})

test_that('Independent and dependent arguments must have names', {
    expect_error(
        objective_function(
            model,
            ddps,
            as.numeric(independent_args),
            quantity_weights,
            data_definitions = data_definitions,
            post_process_function = post_process_function,
            verbose_startup = verbose_startup
        ),
        '`independent_args` must have names'
    )

    expect_error(
        objective_function(
            model,
            ddps,
            independent_args,
            quantity_weights,
            data_definitions = data_definitions,
            post_process_function = post_process_function,
            dependent_arg_function = function(x) {1.0},
            verbose_startup = verbose_startup
        ),
        'The return value of `dependent_arg_function` must have names'
    )

    expect_error(
        objective_function(
            model,
            ddps,
            c(independent_args, list(solar = 1000)),
            quantity_weights,
            data_definitions = data_definitions,
            post_process_function = post_process_function,
            dependent_arg_function = function(x) {list(precip = 0.1)},
            verbose_startup = verbose_startup
        ),
        'Some independent or dependent argument names refer to columns in the drivers: solar, precip'
    )
})

test_that('Bad argument names are detected', {
    expect_error(
        objective_function(
            model,
            ddps,
            c(independent_args, list(bad_arg_name = 1)),
            quantity_weights,
            data_definitions = data_definitions,
            post_process_function = post_process_function,
            verbose_startup = verbose_startup
        ),
        'Model runners could not be created for the following drivers:
ambient_2002: Error: `bad_arg_name` from `arg_names` is not in the `initial_values`, `parameters`, or `drivers`
ambient_2005: Error: `bad_arg_name` from `arg_names` is not in the `initial_values`, `parameters`, or `drivers`',
        fixed = TRUE
    )
})

test_that('Model failures are detected', {
    expect_error(
        objective_function(
            within(model, {parameters$lnfun = 1}),
            ddps,
            independent_args,
            quantity_weights,
            data_definitions = data_definitions,
            post_process_function = post_process_function,
            verbose_startup = verbose_startup
        ),
        'The model could not be run with the following drivers:
ambient_2002: Error in as.data.frame(.Call(R_run_biocro, initial_values, parameters, : Caught exception in R_run_biocro: Thrown by the multilayer_canopy_properties module: lnfun != 0 is not yet supported.
ambient_2005: Error in as.data.frame(.Call(R_run_biocro, initial_values, parameters, : Caught exception in R_run_biocro: Thrown by the multilayer_canopy_properties module: lnfun != 0 is not yet supported.',
        fixed = TRUE
    )
})

test_that('Data-driver pairs must have correct elements', {
    expect_error(
        objective_function(
            model,
            list(),
            independent_args,
            quantity_weights,
            data_definitions = data_definitions,
            post_process_function = post_process_function,
            verbose_startup = verbose_startup
        ),
        '`data_driver_pairs` must have at least one element'
    )

    expect_error(
        objective_function(
            model,
            within(ddps, {ambient_2002$data = NULL; ambient_2005$drivers = NULL}),
            independent_args,
            quantity_weights,
            data_definitions = data_definitions,
            post_process_function = post_process_function,
            verbose_startup = verbose_startup
        ),
        'The following data-driver pairs are missing at least one required element (drivers, data, weight): ambient_2002, ambient_2005',
        fixed = TRUE
    )

    expect_error(
        objective_function(
            model,
            within(ddps, {ambient_2002$extra_element = 5}),
            independent_args,
            quantity_weights,
            data_definitions = data_definitions,
            post_process_function = post_process_function,
            verbose_startup = verbose_startup
        ),
        'The following data-driver pairs have unexpected elements: ambient_2002. The allowed elements are: drivers, data, weight, data_stdev.'
    )

    expect_error(
        objective_function(
            model,
            within(ddps, {ambient_2002$data_stdev = 5}),
            independent_args,
            quantity_weights,
            data_definitions = data_definitions,
            post_process_function = post_process_function,
            verbose_startup = verbose_startup
        ),
        'The following data-driver pairs have a `data_stdev` element that does not match the columns and/or times of their `data` element: ambient_2002'
    )
})

test_that('Data must have a `time` column', {
    expect_error(
        objective_function(
            model,
            within(ddps, {ambient_2002$data$time = NULL}),
            independent_args,
            quantity_weights,
            data_definitions = data_definitions,
            post_process_function = post_process_function,
            verbose_startup = verbose_startup
        ),
        'The following data-driver pairs are missing a `time` column in their `data` element: ambient_2002'
    )
})

test_that('Missing simulation outputs are detected', {
    expect_error(
        objective_function(
            model,
            ddps,
            independent_args,
            quantity_weights,
            data_definitions = data_definitions,
            verbose_startup = verbose_startup
        ),
        'Some data columns were missing from runner outputs:
ambient_2002: Pod
ambient_2005: Pod',
        fixed = TRUE
    )
})

test_that('Out-of-range times are detected', {
    time_offset <- 1e5

    expect_error(
        objective_function(
            model,
            within(ddps, {
                ambient_2002$data$time       <- ambient_2002$data$time       + time_offset
                ambient_2002$data_stdev$time <- ambient_2002$data_stdev$time + time_offset
            }),
            independent_args,
            quantity_weights,
            data_definitions = data_definitions,
            post_process_function = post_process_function,
            verbose_startup = verbose_startup
        ),
        'Some observed times were missing from runner outputs:
ambient_2002: 104272, 104512, 104848, 105184, 105520, 105880, 106192, 106888 (min_time = 3624, max_time = 6911)',
        fixed = TRUE
    )
})

test_that('Multiple time matches are handled', {
    # The drivers have a time step of 1, so if we specify half-integer times,
    # there will actually be two "closest" points to each observed time.
    time_offset <- 0.5

    expect_silent(
        objective_function(
            model,
            within(ddps, {
                ambient_2002$data$time       <- ambient_2002$data$time       + time_offset
                ambient_2002$data_stdev$time <- ambient_2002$data_stdev$time + time_offset
            }),
            independent_args,
            quantity_weights,
            data_definitions = data_definitions,
            post_process_function = post_process_function,
            verbose_startup = verbose_startup
        )
    )
})

test_that('Weights must be supplied for all measured quantities', {
    expect_error(
        objective_function(
            model,
            ddps,
            independent_args,
            list(),
            data_definitions = data_definitions,
            post_process_function = post_process_function,
            verbose_startup = verbose_startup
        ),
        'Weights were not supplied for the following measured quantities: Leaf, Stem, Pod'
    )
})

test_that('Bad normalization methods are detected', {
    expect_error(
        objective_function(
            model,
            ddps,
            independent_args,
            quantity_weights,
            normalization_method = 'bad_normalization_method',
            data_definitions = data_definitions,
            post_process_function = post_process_function,
            verbose_startup = verbose_startup
        ),
        'Unsupported normalization_method: bad_normalization_method'
    )
})

test_that('Bad variance methods are detected', {
    expect_error(
        objective_function(
            model,
            ddps,
            independent_args,
            quantity_weights,
            stdev_weight_method = 'bad_stdev_method',
            data_definitions = data_definitions,
            post_process_function = post_process_function,
            verbose_startup = verbose_startup
        ),
        'Unsupported stdev_weight_method: bad_stdev_method'
    )
})

test_that('Bad return values are detected', {
    expect_error(
        objective_function(
            model,
            ddps,
            independent_args,
            quantity_weights,
            data_definitions = data_definitions,
            post_process_function = post_process_function,
            extra_penalty_function = function(x) {NA},
            verbose_startup = verbose_startup
        ),
        'The objective function did not return a finite value when using the initial argument values; instead, it returned: NA'
    )
})

test_that('Bad regularization methods are detected', {
    expect_error(
        objective_function(
            model,
            ddps,
            independent_args,
            quantity_weights,
            data_definitions = data_definitions,
            post_process_function = post_process_function,
            regularization_method = 'bad_regularization_method',
            verbose_startup = verbose_startup
        ),
        'Unsupported regularization method: bad_regularization_method'
    )
})

test_that('Bad data values and weights are detected', {
    expect_error(
        objective_function(
            model,
            within(ddps, {
                ambient_2005$data_stdev = process_table(soyface_biomass[['ambient_2005_std']])
                ambient_2005$data_stdev[['Leaf_Mg_per_ha']] <- -0.1
            }),
            independent_args,
            quantity_weights,
            stdev_weight_method = 'inverse',
            data_definitions = data_definitions,
            post_process_function = post_process_function,
            verbose_startup = verbose_startup
        ),
        'Issues were found with the following data sets:
  ambient_2002:
  The following columns contained non-finite values: w_var
  ambient_2005:
  The following columns contained non-finite values: w_var
  The following columns contained negative values: quantity_stdev',
        fixed = TRUE
    )
})
