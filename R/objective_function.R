## Here we use several functions that are defined in
## `objective_function_input_checks.R` and `objective_function_helpers.R`

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
