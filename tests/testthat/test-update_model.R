NEWVAL <- 10000

base_model_definition  <- BioCro::soybean
independent_args       <- list(alphaLeaf = NEWVAL, Leaf = NEWVAL)
new_ind_arg_values     <- c(NEWVAL, NEWVAL)
dependent_arg_function <- function(x) {list(alphaStem = x$alphaLeaf)}

test_that('A base model definition can be updated', {
    # Without dependent arguments
    expect_silent(
        update_model(
            base_model_definition,
            independent_args,
            new_ind_arg_values
        )
    )

    # With dependent arguments
    new_model <- expect_silent(
        update_model(
            base_model_definition,
            independent_args,
            new_ind_arg_values,
            dependent_arg_function = dependent_arg_function
        )
    )

    # Initial values are updated
    expect_equal(
        new_model[['initial_values']][['Leaf']],
        NEWVAL
    )

    # Other initial values remain the same
    expect_equal(
        new_model[['initial_values']][['Stem']],
        base_model_definition[['initial_values']][['Stem']]
    )

    # Dependent parameters are updated
    expect_equal(
        new_model[['parameters']][['alphaStem']],
        NEWVAL
    )

    # Other parameters remain the same
    expect_equal(
        new_model[['parameters']][['alphaRoot']],
        base_model_definition[['parameters']][['alphaRoot']]
    )
})

test_that('Base model definition must be valid', {
    expect_error(
        update_model(
            base_model_definition[c('direct_modules', 'differential_modules')],
            independent_args,
            new_ind_arg_values,
            dependent_arg_function = dependent_arg_function
        ),
        'The `base_model_definition` must have the following elements: initial_values, parameters'
    )
})

test_that('Supplied arguments must be part of the base model definition', {
    expect_error(
        update_model(
            base_model_definition,
            c(independent_args, list(bad_arg = 10)),
            new_ind_arg_values
        ),
        '`independent_args` and `new_ind_arg_values` must have the same length'
    )

    expect_error(
        update_model(
            base_model_definition,
            c(independent_args, list(bad_arg = 10)),
            c(new_ind_arg_values, 25)
        ),
        'Values were supplied for the following quantities, but they are not `initial_values` or `parameters` of the `base_model_definition`: bad_arg'
    )
})
