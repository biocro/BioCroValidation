## All the functions defined in this file are intended to perform key operations
## required by `objective_function`.

# Value to return when a simulation fails to run
FAILURE_VALUE <- 1e10

# Helping function for getting a full list of argument names
get_full_arg_names <- function(independent_args, dependent_arg_function) {
    # Get the independent argument names
    independent_arg_names <- names(independent_args)

    # Get the full list of arg_names
    if (is.null(dependent_arg_function)) {
        independent_arg_names
    } else {
        dependent_arg_values <-
            dependent_arg_function(independent_args)

        c(independent_arg_names, names(dependent_arg_values))
    }
}

# Helping function for getting a model runner; if the runner cannot be created,
# an error message will be returned instead
get_model_runner <- function(
    base_model_definition,
    independent_args,
    dependent_arg_function,
    post_process_function,
    ddp
)
{
    # Get the independent argument names
    independent_arg_names <- names(independent_args)

    # Get the full list of arg_names
    arg_names <- get_full_arg_names(independent_args, dependent_arg_function)

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

# Helping function for running each runner with the initial argument values
get_initial_runner_res <- function(model_runners, independent_args) {
    lapply(model_runners, function(runner) {
        tryCatch(
            runner(as.numeric(independent_args)),
            error = function(e) {as.character(e)}
        )
    })
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

# Helping function for converting each data table to a "long form," including
# stdev values when available
get_long_form_data <- function(data_driver_pairs, full_data_definitions) {
    lapply(data_driver_pairs, function(ddp) {
        short_form_data <- ddp[['data']]

        has_std <- 'data_stdev' %in% names(ddp)

        short_form_stdev <- if (has_std) {
            ddp[['data_stdev']]
        } else {
            NA
        }

        data_column_names <- colnames(short_form_data)
        data_column_names <- data_column_names[data_column_names != 'time']

        long_form_data_list <- lapply(data_column_names, function(cn) {
            data.frame(
                time           = short_form_data[, 'time'],
                quantity_name  = full_data_definitions[[cn]],
                quantity_value = short_form_data[, cn],
                quantity_stdev = if (has_std) {short_form_stdev[, cn]} else {1},
                stringsAsFactors = FALSE
            )
        })

        long_form_data <- do.call(rbind, long_form_data_list)

        long_form_data[!is.na(long_form_data[['quantity_value']]), ]
    })
}

# Helping function for getting time indices
add_time_indices <- function(initial_runner_res, long_form_data) {
    for (i in seq_along(long_form_data)) {
        res   <- initial_runner_res[[i]]
        dataf <- long_form_data[[i]]

        indices <- sapply(dataf[, 'time'], function(x) {
            tdiff <- abs(res[, 'time'] - x)

            # Take only the first match, in case there are more
            which(tdiff == min(tdiff))[1]
        })

        long_form_data[[i]][, 'time_index']    <- indices
        long_form_data[[i]][, 'expected_npts'] <- nrow(res)
    }

    long_form_data
}

# Helping function for getting normalization factors
add_norm <- function(long_form_data, normalization_method, n_ddp) {
    for (i in seq_along(long_form_data)) {
        data_table <- long_form_data[[i]]

        data_table[['norm']] <- sapply(seq_len(nrow(data_table)), function(j) {
            qname <- data_table[j, 'quantity_name']

            qname_subset <-
                    data_table[data_table[['quantity_name']] == qname, ]

            if (tolower(normalization_method) == 'equal') {
                1.0
            } else if (tolower(normalization_method) == 'mean') {
                nrow(qname_subset) * n_ddp
            } else if (tolower(normalization_method) == 'max') {
                max(qname_subset[['quantity_value']])^2
            } else if (tolower(normalization_method) == 'mean_max') {
                npts <- nrow(qname_subset)
                qmax <- max(qname_subset[['quantity_value']])
                npts * n_ddp * qmax^2
            } else {
                stop('Unsupported normalization_method: ', normalization_method)
            }
        })

        long_form_data[[i]] <- data_table
    }

    long_form_data
}

# Helping function for getting variance-based weights
add_w_var <- function(long_form_data, stdev_weight_method) {
    for (i in seq_along(long_form_data)) {
        data_table <- long_form_data[[i]]
        data_stdev <- data_table[['quantity_stdev']]

        data_table[['w_var']] <-
            if (tolower(stdev_weight_method) == 'equal') {
                1.0
            } else if (tolower(stdev_weight_method) == 'logarithm') {
                log(1 / (data_stdev + 1e-5))
            } else if (tolower(stdev_weight_method) == 'inverse') {
                1 / data_stdev^2
            } else {
                stop('Unsupported stdev_weight_method: ', stdev_weight_method)
            }

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

# Helping function for getting the data-driver pair weights
get_ddp_weights <- function(data_driver_pairs) {
    lapply(data_driver_pairs, function(ddp) {
        ddp[['weight']]
    })
}

# Helping function that calculates one error
one_error <- function(
    observed,
    predicted,
    quantity_weight,
    ddp_weight,
    var_weight,
    normalization
)
{
    qw <- if (predicted < observed) {
        quantity_weight[1] # Underprediction
    } else {
        quantity_weight[2] # Overprediction
    }

    (observed - predicted)^2 * qw * ddp_weight * var_weight / normalization
}

# Helping function for returning a failure value
failure_value <- function(error_sum, return_terms) {
    if (return_terms) {
        list(
            least_squares_term = error_sum,
            extra_penalty = FAILURE_VALUE
        )
    } else {
        FAILURE_VALUE
    }
}

# Helping function that calculates an error value from a simulation result
error_from_res <- function(
    simulation_result,
    long_form_data_table,
    quantity_weights,
    ddp_weight,
    normalization_method,
    extra_penalty_function,
    return_terms
)
{
    # If the simulation did not finish, return a very high value
    expected_npts <- long_form_data_table[1, 'expected_npts']

    if (nrow(simulation_result) < expected_npts) {
        return(
            failure_value(NA, return_terms)
        )
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
        indx  <- long_form_data_table[i, 'time_index']

        one_error(
            long_form_data_table[i, 'quantity_value'], # obs
            simulation_result[indx, qname],            # pred
            quantity_weights[[qname]],                 # quantity_weight
            ddp_weight,                                # ddp_weight
            long_form_data_table[i, 'w_var'],          # var_weight
            long_form_data_table[i, 'norm']            # norm
        )
    })

    error_sum <- sum(errors)

    # If the error sum is not finite, return a very high value
    if (!is.finite(error_sum)) {
        return(
            failure_value(error_sum, return_terms)
        )
    }

    # Return the sum of the penalty and error terms, or the individual errors
    if (return_terms) {
        error_terms_by_quantity <- as.list(tapply(
            errors,
            long_form_data_table[['quantity_name']],
            sum
        ))

        list(
            least_squares_terms = error_terms_by_quantity,
            extra_penalty = penalty
        )
    } else {
        penalty + error_sum
    }
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
    full_quantity_weights,
    ddp_weights,
    normalization_method,
    extra_penalty_function,
    regularization_method
)
{
    function(x, lambda = 0, return_terms = FALSE) {
        errors <- lapply(seq_along(model_runners), function(i) {
            runner <- model_runners[[i]]
            res    <- runner(x)

            error_from_res(
                res,
                long_form_data[[i]],
                full_quantity_weights,
                ddp_weights[[i]],
                normalization_method,
                extra_penalty_function,
                return_terms
            )
        })

        reg_penalty <- regularization_penalty(x, regularization_method, lambda)

        if (return_terms) {
            list(
                terms_from_data_driver_pairs = stats::setNames(
                    errors,
                    names(model_runners)
                ),
                regularization_penalty = reg_penalty
            )
        } else {
            sum(as.numeric(errors)) + reg_penalty
        }
    }
}
