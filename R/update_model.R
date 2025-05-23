update_model <- function(
    base_model_definition,
    independent_args,
    new_ind_arg_values,
    dependent_arg_function = NULL
)
{
    # Make sure the model definition has initial_values and parameters
    req_elements <- c('initial_values', 'parameters')
    if (!all(req_elements %in% names(base_model_definition))) {
        stop(
            'The `base_model_definition` must have the following elements: ',
            paste(req_elements, collapse = ', ')
        )
    }

    # Make sure the argument lists have the same length
    if (length(independent_args) != length(new_ind_arg_values)) {
        stop('`independent_args` and `new_ind_arg_values` must have the same length')
    }

    # Update the values of the independent arguments
    new_independent_args <- stats::setNames(
        as.list(new_ind_arg_values),
        names(independent_args)
    )

    # Also get the values of the dependent arguments
    all_args <- if (!is.null(dependent_arg_function)) {
        c(new_independent_args, dependent_arg_function(new_independent_args))
    } else {
        new_independent_args
    }

    # Find all quantities in the initial values and parameters and store them in
    # a list
    iv <- as.list(
        rep_len(
            'initial_values',
            length(base_model_definition[['initial_values']])
        )
    )
    names(iv) <- names(base_model_definition[['initial_values']])

    param <- as.list(
        rep_len(
            'parameters',
            length(base_model_definition[['parameters']])
        )
    )
    names(param) <- names(base_model_definition[['parameters']])

    model_quantities <- c(iv, param)

    # Make sure each supplied argument is included in the model
    not_in_model <- !names(all_args) %in% names(model_quantities)

    if (any(not_in_model)) {
        msg <- paste0(
            'Values were supplied for the following quantities, but they ',
            'are not `initial_values` or `parameters` of ',
            'the `base_model_definition`: ',
            paste(names(all_args)[not_in_model], collapse = ', ')
        )
        stop(msg)
    }

    # Make a copy of the model with the new argument values and return it
    new_model_definition <- base_model_definition

    for (i in seq_along(all_args)) {
        arg_name  <- names(all_args)[i]
        arg_type  <- model_quantities[[arg_name]]
        arg_value <- all_args[[i]]

        new_model_definition[[arg_type]][[arg_name]] <- arg_value
    }

    new_model_definition
}
