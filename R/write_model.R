write_model <- function(
    name,
    direct_modules,
    differential_modules,
    initial_values,
    parameters,
    ode_solver
)
{
    # Get the longest name among the initial values and parameters
    all_names <- c(
        names(initial_values),
        names(parameters)
    )

    name_width <- max(sapply(all_names, nchar))

    # Alphabetize the initial values and parameters
    initial_values <- initial_values[order(tolower(names(initial_values)))]
    parameters     <- parameters[order(tolower(names(parameters)))]

    # Fill in the module definition template (defined below)
    model_text <- sprintf(
        model_definition_template,
        name,
        paste(paste0('        "', direct_modules, '"'), collapse = ',\n'),
        paste(paste0('        "', differential_modules, '"'), collapse = ',\n'),
        ode_solver[['type']],
        ode_solver[['output_step_size']],
        ode_solver[['adaptive_rel_error_tol']],
        ode_solver[['adaptive_abs_error_tol']],
        ode_solver[['adaptive_max_steps']],
        paste(paste0('        ', format(names(initial_values), width = name_width), ' = ', initial_values), collapse = ',\n'),
        paste(paste0('        ', format(names(parameters), width = name_width), ' = ', parameters), collapse = ',\n')
    )
}

model_definition_template <- '%s <- list(
    direct_modules = list(
%s
    ),
    differential_modules = list(
%s
    ),
    ode_solver = list(
        type                   = "%s",
        output_step_size       = %f,
        adaptive_rel_error_tol = %e,
        adaptive_abs_error_tol = %e,
        adaptive_max_steps     = %i
    ),
    initial_values = list(
%s
    ),
    parameters = list(
%s
    )
)'
