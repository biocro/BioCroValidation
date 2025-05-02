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

    # Now make sure each runner produces the necessary columns in its output
    expected_columns <- as.character(full_data_definitions)

    missing_columns <- lapply(model_runners, function(runner) {
        runner_result <- runner(as.numeric(initial_independent_arg_values))
        expected_columns[!expected_columns %in% colnames(runner_result)]
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

    return(invisible(NULL))
}

# Helping function for getting a "data definition list," which specifies the
# names of the `data` columns as they appear in the simulation output
get_data_definition_list <- function(data_driver_pairs, user_data_definitions) {
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

objective_function <- function(
    model_definition,
    data_driver_pairs,
    independent_arg_names,
    initial_independent_arg_values,
    data_definitions,
    dependent_arg_function = NULL,
    post_process_function = NULL
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
        model_runners,
        initial_independent_arg_values,
        full_data_definitions
    )
}
