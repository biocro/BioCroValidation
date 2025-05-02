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
    c(alphaLeaf, betaLeaf)
})

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

test_that('Bad argument names are detected', {
    expect_error(
        objective_function(
            model,
            ddps,
            c(independent_arg_names, 'bad_arg_name'),
            initial_independent_arg_values
        ),
        'Model runners could not be created for the following drivers:
ambient_2002: Error: `bad_arg_name` from `arg_names` is not in the `initial_values`, `parameters`, or `drivers`
ambient_2005: Error: `bad_arg_name` from `arg_names` is not in the `initial_values`, `parameters`, or `drivers`',
        fixed = TRUE
    )
})

test_that('Bad initial argument values are detected', {
    expect_error(
        objective_function(
            model,
            ddps,
            independent_arg_names,
            c()
        ),
        'The model could not be run with the following drivers:
ambient_2002: Error in runner(initial_independent_arg_values): The unlisted `x` argument (`unlist(x)`) does not have the correct number of elements: required = 2, actual = 0
ambient_2005: Error in runner(initial_independent_arg_values): The unlisted `x` argument (`unlist(x)`) does not have the correct number of elements: required = 2, actual = 0',
        fixed = TRUE
    )
})
