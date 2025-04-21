write_model <- function(
    name,
    direct_modules,
    differential_modules,
    initial_values,
    parameters,
    ode_solver
)
{
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
        paste(paste0('        "', names(initial_values), '" = ', initial_values), collapse = ',\n'),
        paste(paste0('        "', names(parameters), '" = ', parameters), collapse = ',\n')
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
        type = "%s",
        output_step_size = %f,
        adaptive_rel_error_tol = %e,
        adaptive_abs_error_tol = %e,
        adaptive_max_steps = %i
    ),
    initial_values = list(
%s
    ),
    parameters = list(
%s
    )
)'
