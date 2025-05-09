## Here we use several functions that are defined in
## `objective_function_input_checks.R` and `objective_function_helpers.R`

objective_function <- function(
    base_model_definition,
    data_driver_pairs,
    independent_args,
    quantity_weights,
    data_definitions = list(),
    normalization_method = 'mean_max',
    variance_weight_method = 'equal',
    dependent_arg_function = NULL,
    post_process_function = NULL,
    extra_penalty_function = NULL,
    regularization_method = 'none',
    verbose_startup = FALSE
)
{
    # Check the data-driver pairs
    check_data_driver_pairs(base_model_definition, data_driver_pairs)

    # Check the arguments to be varied
    check_args_to_vary(
        independent_args,
        dependent_arg_function,
        data_driver_pairs
    )

    # Get the model runners
    model_runners <- lapply(data_driver_pairs, function(ddp) {
        get_model_runner(
            base_model_definition,
            independent_args,
            dependent_arg_function,
            post_process_function,
            ddp
        )
    })

    # Get the full data definition list
    full_data_definitions <-
        get_data_definition_list(data_driver_pairs, data_definitions)

    if (verbose_startup) {
        cat('\nThe full data definitions:\n\n')
        utils::str(full_data_definitions)
    }

    # Check the model runners
    check_runners(model_runners)

    # Get initial model runner results
    initial_runner_res <-
        get_initial_runner_res(model_runners, independent_args)

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
    long_form_data <- add_time_indices(initial_runner_res, long_form_data)

    # Add normalization factors
    long_form_data <- add_norm(
        long_form_data,
        normalization_method,
        length(data_driver_pairs)
    )

    # Add variance-based weights
    long_form_data <- add_w_var(
        long_form_data,
        variance_weight_method
    )

    if (verbose_startup) {
        cat('\nThe user-supplied data in long form:\n\n')
        print(long_form_data)
    }

    # Process the quantity weights
    full_quantity_weights <-
        process_quantity_weights(quantity_weights, long_form_data)

    if (verbose_startup) {
        cat('The user-supplied quantity weights:\n\n')
        utils::str(full_quantity_weights)
    }

    # Get the data-driver pair weights
    ddp_weights <- get_ddp_weights(data_driver_pairs)

    if (verbose_startup) {
        cat('\nThe user-supplied data-driver pair weights:\n\n')
        utils::str(ddp_weights)
    }

    # Create the objective function
    obj_fun <- get_obj_fun(
        model_runners,
        long_form_data,
        full_quantity_weights,
        ddp_weights,
        normalization_method,
        extra_penalty_function,
        regularization_method
    )

    # Check the objective function
    check_obj_fun(obj_fun, independent_args, verbose_startup)

    # Return it
    obj_fun
}
