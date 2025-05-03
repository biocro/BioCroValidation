# Helping function for checking the data-driver pairs; will throw an error if
# a problem is detected, and will otherwise be silent with no return value.
check_data_driver_pairs <- function(base_model_definition, data_driver_pairs) {
    if (is.null(names(data_driver_pairs))) {
        stop('`data_driver_pairs` must have names')
    }

    has_elements <- sapply(data_driver_pairs, function(x) {
        'drivers' %in% names(x) && 'data' %in% names(x)
    })

    if (any(!has_elements)) {
        missing_elements <- names(data_driver_pairs)[!has_elements]

        msg <- paste(
            'The following data-driver pairs are missing a `drivers` element,',
            'a `data` element, or both:',
            paste(missing_elements, collapse = ', ')
        )

        stop(msg)
    }

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

# Helping function for checking the independent argument names and initial
# values
check_independent_arguments <- function(
    independent_arg_names,
    initial_ind_arg_values
)
{
    if (length(independent_arg_names) != length(initial_ind_arg_values)) {
        stop(
            '`independent_arg_names` and `initial_ind_arg_values` ',
            'must have the same length'
        )
    }

    if (is.null(names(initial_ind_arg_values))) {
        stop('`initial_ind_arg_values` must have names')
    }

    if (any(!names(initial_ind_arg_values) %in% independent_arg_names)) {
        bad_arg <- !names(initial_ind_arg_values) %in% independent_arg_names

        msg <- paste(
            'The following arguments are included in `initial_ind_arg_values`',
            'but not `independent_arg_names`:',
            paste(names(initial_ind_arg_values)[bad_arg], collapse = ', ')
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

# Helping function for checking the objective function; will throw an error if a
# problem is detected, and will otherwise be silent with no return value.
check_obj_fun <- function(obj_fun, initial_ind_arg_values) {
    initial_error <- obj_fun(as.numeric(initial_ind_arg_values))

    if (!is.finite(initial_error)) {
        stop(
            'The objective function did not return a finite value when using ',
            'the initial argument values; instead, it returned: ',
            initial_error
        )
    }

    return(invisible(NULL))
}
