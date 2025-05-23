independent_args <- list(param1 = 0.1, param2 = 0.2, param3 = 0.3)

bounds_list <- list(
    param1 = c(0, 1),
    param3 = c(4, 3),
    param2 = c(1, 2),
    param4 = c(7, 8)
)

test_that('Bounds tables can be created', {
    # No initial values
    expect_silent(
        bounds_table(
            independent_args,
            bounds_list
        )
    )

    # Initial values
    bounds <- expect_silent(
        bounds_table(
            independent_args,
            bounds_list,
            c(0.5, 1.5, 3.5)
        )
    )

    expect_equal(
        bounds[['arg_name']],
        c('param1', 'param2', 'param3')
    )

    expect_equal(
        bounds[['lower']],
        c(0, 1, 3)
    )

    expect_equal(
        bounds[['upper']],
        c(1, 2, 4)
    )
})

test_that('Values outside the bounds are detected', {
    expect_error(
        bounds_table(
            independent_args,
            bounds_list,
            c(0.5, 2.5, 3.5)
        ),
        'The initial values for the following arguments lie outside the bounds: param2'
    )
})

test_that('Values on the bounds are detected', {
    expect_warning(
        bounds_table(
            independent_args,
            bounds_list,
            c(0.5, 2.0, 3.5)
        ),
        'The initial values for the following arguments lie on the bounds, which can be problematic for some optimizers: param2'
    )
})

test_that('Missing bounds are detected', {
    expect_error(
        bounds_table(
            independent_args,
            1.0
        ),
        '`bounds_list` must be a list of named elements'
    )

    expect_error(
        bounds_table(
            independent_args,
            list(1, 2, 3)
        ),
        '`bounds_list` must be a list of named elements'
    )

    expect_error(
        bounds_table(
            independent_args,
            bounds_list[1:2]
        ),
        'The following elements were included in `independent_args` but not `bounds_list`: param2'
    )

    expect_error(
        bounds_table(
            independent_args,
            within(bounds_list, {param1 = 1.0})
        ),
        'The following elements of `bounds_list` do not have a length of 2: param1'
    )
})
