# Specify key inputs to use for these tests
model <- BioCro::soybean

ddps <- list(
    ambient_2002 = list(
        data = soyface_biomass[['ambient_2002']],
        drivers = BioCro::soybean_weather[['2002']]
    ),
    ambient_2005 = list(
        data = soyface_biomass[['ambient_2005']],
        drivers = BioCro::soybean_weather[['2005']]
    )
)

independent_arg_names <- c('alphaLeaf', 'betaLeaf')

initial_independent_arg_values <- with(BioCro::soybean[['parameters']], {
    list(alphaLeaf = alphaLeaf, betaLeaf = betaLeaf)
})

dependent_arg_function <- function(x) {
    list(alphaStem = x[['alphaLeaf']])
}

# Run tests
test_that('An objective function can be created', {
    obj_fun <- expect_silent(
        objective_function(
            model,
            ddps,
            independent_arg_names,
            initial_independent_arg_values
        )
    )

    obj_fun <- expect_silent(
        objective_function(
            model,
            ddps,
            independent_arg_names,
            initial_independent_arg_values,
            dependent_arg_function
        )
    )
})

test_that('Bad definitions are detected', {
    expect_error(
        objective_function(
            model,
            within(ddps, {ambient_2005$drivers$temp = NULL}),
            independent_arg_names,
            initial_independent_arg_values
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
            c()
        ),
        '`independent_arg_names` and `initial_independent_arg_values` must have the same length'
    )

    expect_error(
        objective_function(
            model,
            ddps,
            independent_arg_names,
            list(arg1 = 1, arg2 = 2)
        ),
        'The following arguments are included in `initial_independent_arg_values` but not `independent_arg_names`: arg1, arg2'
    )

    expect_error(
        objective_function(
            model,
            ddps,
            independent_arg_names,
            as.numeric(initial_independent_arg_values)
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
            c(initial_independent_arg_values, list(bad_arg_name = 1))
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
            initial_independent_arg_values
        ),
        'The model could not be run with the following drivers:
ambient_2002: Error in as.data.frame(.Call(R_run_biocro, initial_values, parameters, : Caught exception in R_run_biocro: Thrown by the multilayer_canopy_properties module: lnfun != 0 is not yet supported.
ambient_2005: Error in as.data.frame(.Call(R_run_biocro, initial_values, parameters, : Caught exception in R_run_biocro: Thrown by the multilayer_canopy_properties module: lnfun != 0 is not yet supported.',
        fixed = TRUE
    )
})
