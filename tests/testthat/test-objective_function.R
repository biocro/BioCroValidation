# Specify key inputs to use for these tests
model <- BioCro::soybean
model$ode_solver <- BioCro::default_ode_solvers[['homemade_euler']]

ddps <- list(
    ambient_2002 = list(
        data = within(soyface_biomass[['ambient_2002']], {time = (DOY - 1) * 24.0; DOY = NULL}),
        drivers = BioCro::soybean_weather[['2002']]
    ),
    ambient_2005 = list(
        data = within(soyface_biomass[['ambient_2005']], {time = (DOY - 1) * 24.0; DOY = NULL}),
        drivers = BioCro::soybean_weather[['2005']]
    )
)

independent_arg_names <- c('alphaLeaf', 'betaLeaf')

initial_independent_arg_values <- with(BioCro::soybean[['parameters']], {
    list(alphaLeaf = alphaLeaf, betaLeaf = betaLeaf)
})

data_definitions <- list(
    Leaf_Mg_per_ha = 'Leaf',
    Stem_Mg_per_ha = 'Stem',
    Pod_Mg_per_ha = 'Pod'
)

dependent_arg_function <- function(x) {
    list(alphaStem = x[['alphaLeaf']])
}

post_process_function <- function(x) {
    within(x, {Pod = Grain + Shell})
}

quantity_weights <- list(
    Leaf = 0.5,
    Stem = 0.5,
    Pod = 1
)

normalization_method <- 'mean_max'

# Run tests
test_that('An objective function can be created', {
    obj_fun <- expect_silent(
        objective_function(
            model,
            ddps,
            independent_arg_names,
            initial_independent_arg_values,
            data_definitions,
            quantity_weights,
            normalization_method,
            post_process_function = post_process_function
        )
    )

    obj_fun <- expect_silent(
        objective_function(
            model,
            ddps[1],
            independent_arg_names,
            initial_independent_arg_values,
            data_definitions,
            quantity_weights,
            normalization_method,
            post_process_function = post_process_function
        )
    )

    obj_fun <- expect_silent(
        objective_function(
            model,
            ddps,
            independent_arg_names,
            initial_independent_arg_values,
            data_definitions,
            quantity_weights,
            normalization_method,
            dependent_arg_function = dependent_arg_function,
            post_process_function = post_process_function
        )
    )
})

test_that('Bad definitions are detected', {
    expect_error(
        objective_function(
            model,
            within(ddps, {ambient_2005$drivers$temp = NULL}),
            independent_arg_names,
            initial_independent_arg_values,
            data_definitions,
            quantity_weights,
            normalization_method,
            post_process_function = post_process_function
        ),
        'The following drivers did not form a valid dynamical system: ambient_2005'
    )
})

test_that('Independent argument names must be consistent', {
    expect_error(
        objective_function(
            model,
            ddps,
            independent_arg_names,
            c(),
            data_definitions,
            quantity_weights,
            normalization_method,
            post_process_function = post_process_function
        ),
        '`independent_arg_names` and `initial_independent_arg_values` must have the same length'
    )

    expect_error(
        objective_function(
            model,
            ddps,
            independent_arg_names,
            list(arg1 = 1, arg2 = 2),
            data_definitions,
            quantity_weights,
            normalization_method,
            post_process_function = post_process_function
        ),
        'The following arguments are included in `initial_independent_arg_values` but not `independent_arg_names`: arg1, arg2'
    )

    expect_error(
        objective_function(
            model,
            ddps,
            independent_arg_names,
            as.numeric(initial_independent_arg_values),
            data_definitions,
            quantity_weights,
            normalization_method,
            post_process_function = post_process_function
        ),
        '`initial_independent_arg_values` must have names'
    )
})

test_that('Bad argument names are detected', {
    expect_error(
        objective_function(
            model,
            ddps,
            c(independent_arg_names, 'bad_arg_name'),
            c(initial_independent_arg_values, list(bad_arg_name = 1)),
            data_definitions,
            quantity_weights,
            normalization_method,
            post_process_function = post_process_function
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
            independent_arg_names,
            initial_independent_arg_values,
            data_definitions,
            quantity_weights,
            normalization_method,
            post_process_function = post_process_function
        ),
        'The model could not be run with the following drivers:
ambient_2002: Error in as.data.frame(.Call(R_run_biocro, initial_values, parameters, : Caught exception in R_run_biocro: Thrown by the multilayer_canopy_properties module: lnfun != 0 is not yet supported.
ambient_2005: Error in as.data.frame(.Call(R_run_biocro, initial_values, parameters, : Caught exception in R_run_biocro: Thrown by the multilayer_canopy_properties module: lnfun != 0 is not yet supported.',
        fixed = TRUE
    )
})

test_that('Data-driver pairs must be complete', {
    expect_error(
        objective_function(
            model,
            within(ddps, {ambient_2002$data = NULL; ambient_2005$drivers = NULL}),
            independent_arg_names,
            initial_independent_arg_values,
            data_definitions,
            quantity_weights,
            normalization_method,
            post_process_function = post_process_function
        ),
        'The following data-driver pairs are missing a `drivers` element, a `data` element, or both: ambient_2002, ambient_2005'
    )
})

test_that('Data must have a `time` column', {
    expect_error(
        objective_function(
            model,
            within(ddps, {ambient_2002$data$time = NULL}),
            independent_arg_names,
            initial_independent_arg_values,
            data_definitions,
            quantity_weights,
            normalization_method,
            post_process_function = post_process_function
        ),
        'The following data-driver pairs are missing a `time` column in their `data` element: ambient_2002'
    )
})

test_that('Missing simulation outputs are detected', {
    expect_error(
        objective_function(
            model,
            ddps,
            independent_arg_names,
            initial_independent_arg_values,
            data_definitions,
            quantity_weights,
            normalization_method
        ),
        'Some data columns were missing from the following runner outputs:
ambient_2002: Pod
ambient_2005: Pod',
        fixed = TRUE
    )
})

test_that('Out-of-range times are detected', {
    expect_error(
        objective_function(
            model,
            within(ddps, {ambient_2002$data$time <- ambient_2002$data$time + 1e5}),
            independent_arg_names,
            initial_independent_arg_values,
            data_definitions,
            quantity_weights,
            normalization_method,
            post_process_function = post_process_function
        ),
        'Some observed times were missing from the following runner outputs:
ambient_2002: 104272, 104512, 104848, 105184, 105520, 105880, 106192, 106888',
        fixed = TRUE
    )
})

test_that('Weights must be supplied for all measured quantities', {
    expect_error(
        objective_function(
            model,
            ddps,
            independent_arg_names,
            initial_independent_arg_values,
            data_definitions,
            list(),
            normalization_method,
            post_process_function = post_process_function
        ),
        'Weights were not supplied for the following measured quantities: Leaf, Stem, Pod'
    )
})

test_that('Bad normalization methods are detected', {
    expect_error(
        objective_function(
            model,
            ddps,
            independent_arg_names,
            initial_independent_arg_values,
            data_definitions,
            quantity_weights,
            'bad_normalization_method',
            post_process_function = post_process_function
        ),
        'Unsupported normalization_method: bad_normalization_method'
    )
})

test_that('Bad return values are detected', {
    expect_error(
        objective_function(
            model,
            ddps,
            independent_arg_names,
            initial_independent_arg_values,
            data_definitions,
            quantity_weights,
            normalization_method,
            post_process_function = post_process_function,
            extra_penalty_function = function(x) {NA}
        ),
        'The objective function did not return a finite value when using the initial argument values; instead, it returned: NA'
    )
})
