# A helper function for checking the bounds list for mistakes; if an issue is
# found, this function will throw an error; otherwise it will be silent with no
# return value.
check_bounds_list <- function(bounds_list, independent_args) {
    # Must be a list of named elements
    if (!is.list(bounds_list) | is.null(names(bounds_list))) {
        stop('`bounds_list` must be a list of named elements')
    }

    # Must contain all elements in independent_args
    missing_element <- sapply(names(independent_args), function(x) {
        !x %in% names(bounds_list)
    })

    if (any(missing_element)) {
        msg <- paste0(
            'The following elements were included in ',
            '`independent_args` but not `bounds_list`: ',
            paste(names(independent_args)[missing_element], collapse = ', ')
        )
        stop(msg)
    }

    # Each element must have length 2
    length_two <- sapply(bounds_list, function(x) {
        xlen <- length(x)

        if (is.finite(xlen)) {
            length(x) == 2
        } else {
            FALSE
        }
    })

    if (any(!length_two)) {
        msg <- paste0(
            'The following elements of `bounds_list` do not have a length of 2: ',
            paste(names(bounds_list)[!length_two], collapse = ', ')
        )
        stop(msg)
    }

    # Each element must be numeric
    not_numeric <- sapply(bounds_list, function(x) {!is.numeric(x)})

    if (any(not_numeric)) {
        msg <- paste0(
            'The following elements of `bounds_list` are not numeric: ',
            paste(names(bounds_list)[not_numeric], collapse = ', ')
        )
        stop(msg)
    }

    return(invisible(NULL))
}

# A helper function for checking the initial guess for mistakes; if an issue is
# found, this function will throw an error or a warning; otherwise it will be
# silent with no return value.
check_initial_ind_arg_values <- function(
    independent_args,
    lbounds,
    ubounds,
    initial_ind_arg_values
)
{
    # Check the length
    if (length(initial_ind_arg_values) != length(independent_args)) {
        stop('`initial_ind_arg_values` must have the same length as `independent_args`')
    }

    # Check to make sure the initial values are not outside the bounds
    outside_bounds <- sapply(seq_along(initial_ind_arg_values), function(i) {
        initial_ind_arg_values[i] < lbounds[i] | initial_ind_arg_values[i] > ubounds[i]
    })

    if (any(outside_bounds)) {
        msg <- paste0(
            'The initial values for the following arguments lie outside the bounds: ',
            paste(names(independent_args)[outside_bounds], collapse = ', ')
        )
        stop(msg)
    }

    # Check to see if any initial values are on the bounds
    eps <- sqrt(.Machine$double.eps)

    on_bounds <- sapply(seq_along(initial_ind_arg_values), function(i) {
        abs(initial_ind_arg_values[i] - lbounds[i]) < eps |
            abs(initial_ind_arg_values[i] - ubounds[i]) < eps
    })

    if (any(on_bounds)) {
        msg <- paste0(
            'The initial values for the following arguments lie on the ',
            'bounds, which can be problematic for some optimizers: ',
            paste(names(independent_args)[on_bounds], collapse = ', ')
        )
        warning(msg)
    }

    return(invisible(NULL))
}

bounds_table <- function(
    independent_args,
    bounds_list,
    initial_ind_arg_values = NULL
)
{
    # Check the bounds_list
    check_bounds_list(bounds_list, independent_args)

    # Get an ordering for the elements of `bounds_list` so they match the order
    # of elements in `independent_args`; note that this will also exclude any
    # elements of `bounds_list` that are not included in `independent_args`.
    ordering <- match(
        names(independent_args),
        names(bounds_list)
    )

    bounds_list <- bounds_list[ordering]

    # Get the lower and upper bounds
    lbounds <- sapply(bounds_list, min)
    ubounds <- sapply(bounds_list, max)

    # Form the bounds table
    bounds_table <- data.frame(
        arg_name = names(independent_args),
        lower = lbounds,
        upper = ubounds,
        stringsAsFactors = FALSE
    )

    # Include initial values in the table if they were provided
    if (!is.null(initial_ind_arg_values)) {
        # Check the values
        check_initial_ind_arg_values(
            independent_args,
            lbounds,
            ubounds,
            initial_ind_arg_values
        )

        # Include the values
        bounds_table$initial_value <- initial_ind_arg_values
    }

    # Remove row names and return the table
    rownames(bounds_table) <- NULL

    bounds_table
}
