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

# Helping function for getting a model runner; if the runner cannot be created,
# an error message will be returned instead
get_model_runner <- function(
    base_model_definition,
    independent_arg_names,
    initial_ind_arg_values,
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
            dependent_arg_function(initial_ind_arg_values)

        c(independent_arg_names, names(dependent_arg_values))
    }

    # Build the runner
    tryCatch({
            partial_func <- BioCro::partial_run_biocro(
                base_model_definition[['initial_values']],
                base_model_definition[['parameters']],
                ddp[['drivers']],
                base_model_definition[['direct_modules']],
                base_model_definition[['differential_modules']],
                base_model_definition[['ode_solver']],
                arg_names
            )

            function(x) {
                if (!is.numeric(x)) {
                    stop('`x` must be numeric')
                }

                x_for_partial <- if (is.null(dependent_arg_function)) {
                    x
                } else {
                    x_for_dep_arg_func <-
                        stats::setNames(x, independent_arg_names)

                    c(x, as.numeric(dependent_arg_function(x_for_dep_arg_func)))
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

# Helping function for running each runner with the initial argument values
get_initial_runner_res <- function(model_runners, initial_ind_arg_values) {
    lapply(model_runners, function(runner) {
        tryCatch(
            runner(as.numeric(initial_ind_arg_values)),
            error = function(e) {as.character(e)}
        )
    })
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
    initial_runner_res,
    initial_ind_arg_values,
    long_form_data
)
{
    for (i in seq_along(long_form_data)) {
        res   <- initial_runner_res[[i]]
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

# Helping function for getting normalization factors
add_norm <- function(long_form_data, normalization_method) {
    for (i in seq_along(long_form_data)) {
        data_table <- long_form_data[[i]]

        data_table[['norm']] <- sapply(seq_len(nrow(data_table)), function(j) {
            qname <- data_table[j, 'quantity_name']

            qname_subset <-
                    data_table[data_table[['quantity_name']] == qname, ]

            if (tolower(normalization_method) == 'mean_max') {
                npts <- nrow(qname_subset)
                qmax <- max(qname_subset[['quantity_value']])
                npts * qmax^2
            } else {
                stop('Unsupported normalization_method: ', normalization_method)
            }
        })

        long_form_data[[i]] <- data_table
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

# Helping function that calculates one error
one_error <- function(observed, predicted, weight, normalization) {
    weight_multiplier <- if (predicted < observed) {
        weight[1] # Underprediction
    } else {
        weight[2] # Overprediction
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
        norm  <- long_form_data_table[i, 'norm']

        one_error(obs, pred, wt, norm)
    })

    # Return the sum of the penalty and error terms
    penalty + sum(errors)
}

# Helping function for calculating a regularization penalty term
regularization_penalty <- function(
    ind_arg_vals,
    regularization_method,
    regularization_lambda
)
{
    if (toupper(regularization_method) == 'NONE') {
        0.0
    } else if (toupper(regularization_method) == 'LASSO' || toupper(regularization_method) == 'L1') {
        regularization_lambda * sum(abs(ind_arg_vals))
    } else if (toupper(regularization_method) == 'RIDGE' || toupper(regularization_method) == 'L2') {
        regularization_lambda * sum(ind_arg_vals^2)
    } else {
        stop('Unsupported regularization method: ', regularization_method)
    }
}

# Helping function that forms the overall objective function
get_obj_fun <- function(
    model_runners,
    long_form_data,
    processed_weights,
    normalization_method,
    extra_penalty_function,
    regularization_method
)
{
    function(x, lambda = 0) {
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

        reg_penalty <- regularization_penalty(x, regularization_method, lambda)

        sum(errors)
    }
}

objective_function <- function(
    base_model_definition,
    data_driver_pairs,
    independent_arg_names,
    initial_ind_arg_values,
    quantity_weights,
    data_definitions = list(),
    normalization_method = 'mean_max',
    dependent_arg_function = NULL,
    post_process_function = NULL,
    extra_penalty_function = NULL,
    regularization_method = 'none'
)
{
    # Check the data-driver pairs
    check_data_driver_pairs(base_model_definition, data_driver_pairs)

    # Check the independent arguments
    check_independent_arguments(
        independent_arg_names,
        initial_ind_arg_values
    )

    # Get the model runners
    model_runners <- lapply(data_driver_pairs, function(ddp) {
        get_model_runner(
            base_model_definition,
            independent_arg_names,
            initial_ind_arg_values,
            dependent_arg_function,
            post_process_function,
            ddp
        )
    })

    # Get the full data definition list
    full_data_definitions <-
        get_data_definition_list(data_driver_pairs, data_definitions)

    # Check the model runners
    check_runners(model_runners)

    # Get initial model runner results
    initial_runner_res <-
        get_initial_runner_res(model_runners, initial_ind_arg_values)

    # Check the initial model runner results
    check_runner_results(
        initial_runner_res,
        full_data_definitions,
        data_driver_pairs
    )

    # Get the long-form data
    long_form_data <-
        get_long_form_data(data_driver_pairs, full_data_definitions)

    # Find indices corresponding to the measured time points
    long_form_data <- add_time_indices(
        initial_runner_res,
        initial_ind_arg_values,
        long_form_data
    )

    # Add normalization factors
    long_form_data <- add_norm(long_form_data, normalization_method)

    # Process the quantity weights
    processed_weights <-
        process_quantity_weights(quantity_weights, long_form_data)

    # Create the objective function
    obj_fun <- get_obj_fun(
        model_runners,
        long_form_data,
        processed_weights,
        normalization_method,
        extra_penalty_function,
        regularization_method
    )

    # Check the objective function
    check_obj_fun(obj_fun, initial_ind_arg_values)

    # Return it
    obj_fun
}
