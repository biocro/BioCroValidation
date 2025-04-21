test_model <- list(
    direct_modules = list(
        "BioCro:harmonic_energy"
    ),
    differential_modules = list(
        "BioCro:harmonic_oscillator"
    ),
    ode_solver = list(
        type = "boost_rkck54",
        output_step_size = 1.000000,
        adaptive_rel_error_tol = 1.000000e-04,
        adaptive_abs_error_tol = 1.000000e-04,
        adaptive_max_steps = 200
    ),
    initial_values = list(
        "position" = 1,
        "velocity" = 0
    ),
    parameters = list(
        "timestep" = 1,
        "mass" = 1,
        "spring_constant" = 0.5
    )
)
