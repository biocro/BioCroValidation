## All the functions defined in this file are intended to check inputs for
## certain errors. Each function will throw an error if a problem is detected,
## and will otherwise be silent with no return value.

# Helping function for checking the data-driver pairs; will throw an error if
# a problem is detected, and will otherwise be silent with no return value.
check_data_driver_pairs <- function(base_model_definition, data_driver_pairs) {
    # There must be at least one data-driver pair
    if (length(data_driver_pairs) < 1) {
        stop('`data_driver_pairs` must have at least one element')
    }

    # Data-driver pairs must have names
    if (is.null(names(data_driver_pairs))) {
        stop('`data_driver_pairs` must have names')
    }

    # Each data-driver pair must have the required elements
    required_elements <- c('drivers', 'data', 'weight')

    has_elements <- sapply(data_driver_pairs, function(x) {
        all(required_elements %in% names(x))
    })

    if (any(!has_elements)) {
        missing_elements <- names(data_driver_pairs)[!has_elements]

        msg <- paste0(
            'The following data-driver pairs are missing at least one ',
            'required element (',
            paste(required_elements, collapse = ', '),
            '): ',
            paste(missing_elements, collapse = ', ')
        )

        stop(msg)
    }

    # Only required or optional elements should be provided
    optional_elements <- 'data_stdev'

    acceptable_elements <- c(required_elements, optional_elements)

    has_extra_elements <- sapply(data_driver_pairs, function(x) {
        any(!names(x) %in% acceptable_elements)
    })

    if (any(has_extra_elements)) {
        bad_elements <- names(data_driver_pairs)[has_extra_elements]

        msg <- paste0(
            'The following data-driver pairs have unexpected elements: ',
            paste(bad_elements, collapse = ', '),
            '. The allowed elements are: ',
            paste(acceptable_elements, collapse = ', '),
            '.'
        )

        stop(msg)
    }

    # Each data table must have a time column
    has_time <- sapply(data_driver_pairs, function(x) {
        'time' %in% colnames(x[['data']])
    })

    if (any(!has_time)) {
        missing_time <- names(data_driver_pairs)[!has_time]

        msg <- paste(
            'The following data-driver pairs are missing a `time` column',
            'in their `data` element:',
            paste(missing_time, collapse = ', ')
        )

        stop(msg)
    }

    # If provided, stdev tables must have the same columns and time values as
    # their corresponding data tables. Time values must also be in the same
    # order.
    stdev_okay <- sapply(data_driver_pairs, function(x) {
        if ('data_stdev' %in% names(x)) {
            data_table  <- x[['data']]
            stdev_table <- x[['data_stdev']]

            if (is.null(colnames(stdev_table))) {
                FALSE
            } else {
                colnames_match <- identical(
                    sort(colnames(data_table)),
                    sort(colnames(stdev_table))
                )

                times_match <- isTRUE(all.equal(
                    data_table[['time']],
                    stdev_table[['time']]
                ))

                colnames_match && times_match
            }
        } else {
            TRUE
        }
    })

    if (any(!stdev_okay)) {
        bad_stdev <- names(data_driver_pairs)[!stdev_okay]

        msg <- paste(
            'The following data-driver pairs have a `data_stdev` element',
            'that does not match the columns and/or times of their',
            '`data` element:',
            paste(bad_stdev, collapse = ', ')
        )

        stop(msg)
    }

    # Each set of drivers must form a valid dynamical system along with the
    # base model definition
    valid_definitions <- sapply(data_driver_pairs, function(ddp) {
        BioCro::validate_dynamical_system_inputs(
            base_model_definition[['initial_values']],
            base_model_definition[['parameters']],
            ddp[['drivers']],
            base_model_definition[['direct_modules']],
            base_model_definition[['differential_modules']],
            verbose = FALSE
        )
    })

    if (any(!valid_definitions)) {
        invalid_ddp <- names(data_driver_pairs)[!valid_definitions]

        msg <- paste(
            'The following drivers did not form a valid dynamical system:',
            paste(invalid_ddp, collapse = ', ')
        )

        stop(msg)
    }

    return(invisible(NULL))
}

# Helping function for checking the independent arguments
check_args_to_vary <- function(
    independent_args,
    dependent_arg_function,
    data_driver_pairs
)
{
    # Make sure the independent arguments have names
    ind_arg_names <- names(independent_args)

    if (is.null(ind_arg_names)) {
        stop('`independent_args` must have names')
    }

    # Make sure the dependent argument function returns a named list
    if (!is.null(dependent_arg_function)) {
        dep_arg_names <-
            names(dependent_arg_function(independent_args))

        if (is.null(dep_arg_names)) {
            stop('The return value of `dependent_arg_function` must have names')
        }
    }

    # Make sure no drivers were specified
    arg_names <- get_full_arg_names(independent_args, dependent_arg_function)

    args_in_drivers <- lapply(data_driver_pairs, function(ddp) {
        driver_names <- names(ddp[['drivers']])
        arg_names[arg_names %in% driver_names]
    })

    args_in_drivers <- unique(unlist(args_in_drivers))

    if (length(args_in_drivers) > 0) {
        msg <- paste(
            'Some independent or dependent argument names refer to columns',
            'in the drivers:',
            paste(args_in_drivers, collapse = ', ')
        )

        stop(msg)
    }

    return(invisible(NULL))
}

# Helping function for checking the model runners; will throw an error if a
# problem is detected, and will otherwise be silent with no return value.
check_runners <- function(model_runners) {
    bad_runners <- sapply(model_runners, is.character)

    if (any(bad_runners)) {
        bad_runner_names <- names(model_runners)[bad_runners]
        bad_runner_msg <- as.character(model_runners[bad_runners])

        msg <- paste0(
            'Model runners could not be created for the following drivers:\n',
            paste0(bad_runner_names, ': ', bad_runner_msg, collapse = '')
        )

        stop(msg)
    }

    return(invisible(NULL))
}

# Helping function for checking the initial model runner results; will throw an
# error if a problem is detected, and will otherwise be silent with no return
# value.
check_runner_results <- function(
    initial_runner_res,
    full_data_definitions,
    data_driver_pairs
)
{
    # Check for runners that could not be evaluated
    bad_res <- sapply(initial_runner_res, is.character)

    if (any(bad_res)) {
        bad_runner_names <- names(initial_runner_res)[bad_res]
        bad_runner_msg   <- initial_runner_res[bad_res]

        msg <- paste0(
            'The model could not be run with the following drivers:\n',
            paste0(bad_runner_names, ': ', bad_runner_msg, collapse = '')
        )

        stop(msg)
    }

    # Make sure each runner produces the necessary columns in its output
    expected_columns <- as.character(full_data_definitions)

    missing_columns <- lapply(initial_runner_res, function(res) {
        expected_columns[!expected_columns %in% colnames(res)]
    })

    bad_outputs <- sapply(missing_columns, function(x) {
        length(x) > 0
    })

    if (any(bad_outputs)) {
        msg <- 'Some data columns were missing from runner outputs:'

        for (i in seq_along(bad_outputs)) {
            if (bad_outputs[i]) {
                msg <- append(
                    msg,
                    paste0(
                        names(initial_runner_res)[i], ': ',
                        paste(missing_columns[[i]], collapse = ', ')
                    )
                )
            }
        }

        stop(paste(msg, collapse = '\n'))
    }

    # Make sure the output from each runner includes the observed times
    times_out_of_range <- lapply(seq_along(initial_runner_res), function(i) {
        res <- initial_runner_res[[i]]

        min_time <- min(res[['time']])
        max_time <- max(res[['time']])

        data_times <- data_driver_pairs[[i]][['data']][['time']]

        oor <- sapply(data_times, function(datat) {
            datat < min_time || datat > max_time
        })

        data_times[oor]
    })

    bad_times <- sapply(times_out_of_range, function(x) {
        length(x) > 0
    })

    if (any(bad_times)) {
        msg <- 'Some observed times were missing from runner outputs:'

        for (i in seq_along(bad_times)) {
            if (bad_times[i]) {
                msg <- append(
                    msg,
                    paste0(
                        names(initial_runner_res)[i], ': ',
                        paste(times_out_of_range[[i]], collapse = ', ')
                    )
                )
            }
        }

        stop(paste(msg, collapse = '\n'))
    }

    return(invisible(NULL))
}

# Helping function for checking the long-form data; will throw an error if a
# problem is detected, and will otherwise be silent with no return value.
check_long_form_data <- function(long_form_data) {
    # Check each element for issues
    messages <- sapply(long_form_data, function(lfd) {
        msg <- character()

        # Check that certain columns have finite values
        check_for_not_finite <- c(
            'time',
            'quantity_value',
            'quantity_stdev',
            'time_index',
            'expected_npts',
            'norm',
            'w_var'
        )

        not_finite <- sapply(check_for_not_finite, function(cn) {
            any(!is.finite(lfd[, cn]))
        })

        if (any(not_finite)) {
            not_finite_col <- check_for_not_finite[not_finite]

            new_msg <- paste(
                'The following columns contained non-finite values:',
                paste(not_finite_col, collapse = ', ')
            )

            msg <- append(msg, new_msg)
        }

        # Check that certain columns do not have negative values
        check_for_negative <- c(
            'quantity_stdev',
            'time_index',
            'expected_npts',
            'norm',
            'w_var'
        )

        negative <- sapply(check_for_negative, function(cn) {
            any(lfd[, cn] < 0)
        })

        if (any(negative)) {
            negative_col <- check_for_negative[negative]

            new_msg <- paste(
                'The following columns contained negative values:',
                paste(negative_col, collapse = ', ')
            )

            msg <- append(msg, new_msg)
        }

        # Return any messages
        paste(msg, collapse = '\n  ')
    })

    # Send a message if problems were detected
    error_found <- messages != ''

    if (any(error_found)) {
        data_names     <- names(long_form_data)[error_found]
        error_messages <- messages[error_found]

        msg <- paste0(
            'Issues were found with the following data sets:\n  ',
            paste0(
                data_names, ':\n  ',
                error_messages,
                collapse = '\n  '
            )
        )

        stop(msg)
    }

    return(invisible(NULL))
}

# Helping function for checking the objective function; will throw an error if a
# problem is detected, and will otherwise be silent with no return value.
check_obj_fun <- function(obj_fun, initial_ind_arg_values, verbose) {
    initial_error_terms <-
        obj_fun(as.numeric(initial_ind_arg_values), return_terms = TRUE)

    if (verbose) {
        cat('\nThe initial error metric terms:\n\n')
        utils::str(initial_error_terms)
    }

    initial_error <- sum(unlist(initial_error_terms))

    if (!is.finite(initial_error)) {
        stop(
            'The objective function did not return a finite value when using ',
            'the initial argument values; instead, it returned: ',
            initial_error
        )
    }

    return(invisible(NULL))
}
