# Helping function for checking the data-driver pairs; will throw an error if
# a problem is detected, and will otherwise be silent with no return value.
check_data_driver_pairs <- function(model_definition, data_driver_pairs) {
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
            model_definition[['initial_values']],
            model_definition[['parameters']],
            ddp[['drivers']],
            model_definition[['direct_modules']],
            model_definition[['differential_modules']],
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
    initial_independent_arg_values
)
{
    if (length(independent_arg_names) != length(initial_independent_arg_values)) {
        stop('`independent_arg_names` and `initial_independent_arg_values` must have the same length')
    }

    if (is.null(names(initial_independent_arg_values))) {
        stop('`initial_independent_arg_values` must have names')
    }

    if (any(!names(initial_independent_arg_values) %in% independent_arg_names)) {
        bad_arg <- !names(initial_independent_arg_values) %in% independent_arg_names

        msg <- paste(
            'The following arguments are included in `initial_independent_arg_values`',
            'but not `independent_arg_names`:',
            paste(names(initial_independent_arg_values)[bad_arg], collapse = ', ')
        )

        stop(msg)
    }

    return(invisible(NULL))
}

# Helping function for getting a model runner; if the runner cannot be created,
# an error message will be returned instead
get_model_runner <- function(
    model_definition,
    independent_arg_names,
    initial_independent_arg_values,
    dependent_arg_function,
    post_process_function,
    ddp
)
{
    # Get the full list of arg_names
    arg_names <- if (is.null(dependent_arg_function)) {
        independent_arg_names
    } else {
        dependent_arg_values <-
            dependent_arg_function(initial_independent_arg_values)

        c(independent_arg_names, names(dependent_arg_values))
    }

    # Build the runner
    tryCatch({
            partial_func <- BioCro::partial_run_biocro(
                model_definition[['initial_values']],
                model_definition[['parameters']],
                ddp[['drivers']],
                model_definition[['direct_modules']],
                model_definition[['differential_modules']],
                model_definition[['ode_solver']],
                arg_names
            )

            function(x) {
                if (!is.numeric(x)) {
                    stop('`x` must be numeric')
                }

                x_for_partial <- if (is.null(dependent_arg_function)) {
                    x
                } else {
                    x_for_dependent_arg_func <-
                        stats::setNames(x, independent_arg_names)

                    c(x, as.numeric(dependent_arg_function(x_for_dependent_arg_func)))
                }

                initial_res <- partial_func(x_for_partial)

                if (is.null(post_process_function)) {
                    initial_res
                } else {
                    post_process_function(initial_res)
                }
            }
        },
        error = function(e) {as.character(e)}
    )
}

# Helping function for checking the model runners; will throw an error if a
# problem is detected, and will otherwise be silent with no return value.
check_runners <- function(
    data_driver_pairs,
    model_runners,
    initial_independent_arg_values,
    full_data_definitions
)
{
    # First check for runners that could not be created
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

    # Now check for runners that cannot be evaluated
    runner_eval_msg <- sapply(model_runners, function(runner) {
        runner_result <- tryCatch(
            runner(as.numeric(initial_independent_arg_values)),
            error = function(e) {as.character(e)}
        )

        if (is.character(runner_result)) {
            runner_result
        } else {
            ''
        }
    })

    if (any(runner_eval_msg != '')) {
        bad_runner_names <- names(model_runners)[runner_eval_msg != '']
        bad_runner_msg <- runner_eval_msg[runner_eval_msg != '']

        msg <- paste0(
            'The model could not be run with the following drivers:\n',
            paste0(bad_runner_names, ': ', bad_runner_msg, collapse = '')
        )

        stop(msg)
    }

    # Run each runner
    runner_results <- lapply(model_runners, function(runner) {
        runner(as.numeric(initial_independent_arg_values))
    })

    # Now make sure each runner produces the necessary columns in its output
    expected_columns <- as.character(full_data_definitions)

    missing_columns <- lapply(runner_results, function(res) {
        expected_columns[!expected_columns %in% colnames(res)]
    })

    bad_outputs <- sapply(missing_columns, function(x) {
        length(x) > 0
    })

    if (any(bad_outputs)) {
        msg <- 'Some data columns were missing from the following runner outputs:'

        for (i in seq_along(bad_outputs)) {
            if (bad_outputs[i]) {
                msg <- append(
                    msg,
                    paste0(
                        names(model_runners)[i], ': ',
                        paste(missing_columns[[i]], collapse = ', ')
                    )
                )
            }
        }

        stop(paste(msg, collapse = '\n'))
    }

    # Make sure the output from each runner includes the observed times
    times_out_of_range <- lapply(seq_along(runner_results), function(i) {
        res <- runner_results[[i]]

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
        msg <- 'Some observed times were missing from the following runner outputs:'

        for (i in seq_along(bad_times)) {
            if (bad_times[i]) {
                msg <- append(
                    msg,
                    paste0(
                        names(model_runners)[i], ': ',
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
check_obj_fun <- function(obj_fun, initial_independent_arg_values) {
    initial_error <- obj_fun(as.numeric(initial_independent_arg_values))

    if (!is.finite(initial_error)) {
        stop(
            'The objective function did not return a finite value when using ',
            'the initial argument values; instead, it returned: ',
            initial_error
        )
    }

    return(invisible(NULL))
}

# Helping function for getting a "data definition list," which specifies the
# names of the `data` columns as they appear in the simulation output
get_data_definition_list <- function(data_driver_pairs, user_data_definitions)
{
    # First get all the column names found in the observed data
    all_data_colnames <-
        lapply(data_driver_pairs, function(x) {colnames(x[['data']])})

    all_data_colnames <- unlist(all_data_colnames)

    all_data_colnames <- all_data_colnames[!duplicated(all_data_colnames)]

    # Remove the `time` column
    all_data_colnames <- all_data_colnames[all_data_colnames != 'time']

    # Build the data definition list
    data_definitions <- lapply(all_data_colnames, function(cn) {
        if (cn %in% names(user_data_definitions)) {
            user_data_definitions[[cn]]
        } else {
            cn
        }
    })

    names(data_definitions) <- all_data_colnames

    data_definitions
}

# Helping function for converting each data table to a "long form."
get_long_form_data <- function(data_driver_pairs, full_data_definitions) {
    lapply(data_driver_pairs, function(ddp) {
        short_form_data <- ddp[['data']]

        data_column_names <- colnames(short_form_data)
        data_column_names <- data_column_names[data_column_names != 'time']

        long_form_data_list <- lapply(data_column_names, function(cn) {
            data.frame(
                time = short_form_data[, 'time'],
                quantity_name = full_data_definitions[[cn]],
                quantity_value = short_form_data[, cn],
                stringsAsFactors = FALSE
            )
        })

        long_form_data <- do.call(rbind, long_form_data_list)

        long_form_data[!is.na(long_form_data[['quantity_value']]), ]
    })
}

# Helping function for getting time indices
add_time_indices <- function(
    model_runners,
    initial_independent_arg_values,
    long_form_data
)
{
    for (i in seq_along(long_form_data)) {
        runner <- model_runners[[i]]
        res    <- runner(as.numeric(initial_independent_arg_values))

        dataf <- long_form_data[[i]]

        indices <- sapply(dataf[, 'time'], function(x) {
            tdiff <- abs(res[, 'time'] - x)
            which(tdiff == min(tdiff))
        })

        long_form_data[[i]][, 'time_index']    <- indices
        long_form_data[[i]][, 'expected_npts'] <- nrow(res)
    }

    long_form_data
}

# Helping function that processes and checks the quantity weights
process_quantity_weights <- function(quantity_weights, long_form_data) {
    # First make sure that weights have been provided for all measured
    # quantities
    all_data_colnames <- lapply(long_form_data, function(x) {
        unique(x[, 'quantity_name'])
    })

    all_data_colnames <- unlist(all_data_colnames)

    all_data_colnames <- unique(all_data_colnames)

    weight_was_supplied <- sapply(all_data_colnames, function(cn) {
        cn %in% names(quantity_weights)
    })

    if (any(!weight_was_supplied)) {
        missing_weights <- all_data_colnames[!weight_was_supplied]

        msg <- paste(
            'Weights were not supplied for the following measured quantities:',
            paste(missing_weights, collapse = ', ')
        )

        stop(msg)
    }

    # Now make sure all the weights have length 2
    lapply(quantity_weights, function(wt) {
        rep_len(wt, 2)
    })
}

# Helping function that calculates one normalization factor
one_norm <- function(long_form_data_table, qname, normalization_method) {
    if (tolower(normalization_method) == 'mean_max') {
        npts <- sum(long_form_data_table[, 'quantity_name'] == qname)
        qmax <- max(long_form_data_table[long_form_data_table[, 'quantity_name'] == qname, 'quantity_value'])
        npts * qmax^2
    } else {
        stop('Unsupported normalization_method: ', normalization_method)
    }
}

# Helping function that calculates one error
one_error <- function(observed, predicted, weight, normalization) {
    weight_multiplier <- if (observed <= predicted) {
        weight[1]
    } else {
        weight[2]
    }

    (observed - predicted)^2 * weight_multiplier / normalization
}

# Helping function that calculates an error value from a simulation result
error_from_res <- function(
    simulation_result,
    long_form_data_table,
    quantity_weights,
    normalization_method,
    extra_penalty_function
)
{
    # If the simulation did not finish, return a very high value
    expected_npts <- long_form_data_table[1, 'expected_npts']

    if (nrow(simulation_result) < expected_npts) {
        return(1e6)
    }

    # Calculate any user-specified penalties
    penalty <- if (is.null(extra_penalty_function)) {
        0.0
    } else {
        extra_penalty_function(simulation_result)
    }

    # Calculate the error terms
    n_obs <- nrow(long_form_data_table)

    errors <- sapply(seq_len(n_obs), function(i) {
        qname <- as.character(long_form_data_table[i, 'quantity_name'])
        obs   <- long_form_data_table[i, 'quantity_value']
        indx  <- long_form_data_table[i, 'time_index']
        pred  <- simulation_result[indx, qname]
        wt    <- quantity_weights[[qname]]
        norm  <- one_norm(long_form_data_table, qname, normalization_method)

        one_error(obs, pred, wt, norm)
    })

    # Return the sum of the penalty and error terms
    penalty + sum(errors)
}

# Helping function that forms the overall objective function
get_obj_fun <- function(
    model_runners,
    long_form_data,
    processed_weights,
    normalization_method,
    extra_penalty_function
)
{
    function(x) {
        errors <- sapply(seq_along(model_runners), function(i) {
            runner <- model_runners[[i]]
            res    <- runner(x)

            error_from_res(
                res,
                long_form_data[[i]],
                processed_weights,
                normalization_method,
                extra_penalty_function
            )
        })

        sum(errors)
    }
}

objective_function <- function(
    model_definition,
    data_driver_pairs,
    independent_arg_names,
    initial_independent_arg_values,
    quantity_weights,
    data_definitions = list(),
    normalization_method = 'mean_max',
    dependent_arg_function = NULL,
    post_process_function = NULL,
    extra_penalty_function = NULL
)
{
    # Check the data-driver pairs
    check_data_driver_pairs(model_definition, data_driver_pairs)

    # Check the independent arguments
    check_independent_arguments(
        independent_arg_names,
        initial_independent_arg_values
    )

    # Get the model runners
    model_runners <- lapply(data_driver_pairs, function(ddp) {
        get_model_runner(
            model_definition,
            independent_arg_names,
            initial_independent_arg_values,
            dependent_arg_function,
            post_process_function,
            ddp
        )
    })

    # Get the full data definition list
    full_data_definitions <-
        get_data_definition_list(data_driver_pairs, data_definitions)

    # Check the model runners
    check_runners(
        data_driver_pairs,
        model_runners,
        initial_independent_arg_values,
        full_data_definitions
    )

    # Get the long-form data
    long_form_data <-
        get_long_form_data(data_driver_pairs, full_data_definitions)

    # Find indices corresponding to the measured time points
    long_form_data <- add_time_indices(
        model_runners,
        initial_independent_arg_values,
        long_form_data
    )

    # Process the quantity weights
    processed_weights <-
        process_quantity_weights(quantity_weights, long_form_data)

    # Create the objective function
    obj_fun <- get_obj_fun(
        model_runners,
        long_form_data,
        processed_weights,
        normalization_method,
        extra_penalty_function
    )

    # Check the objective function
    check_obj_fun(obj_fun, initial_independent_arg_values)

    # Return it
    obj_fun
}
