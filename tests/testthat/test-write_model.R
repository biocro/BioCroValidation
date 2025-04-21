test_that('Model definitions are written correctly', {
    # Get stored results
    stored_fname <- file.path('..', 'test_data', 'test_model.R')
    stored_file <- paste(readLines(stored_fname), collapse='\n')

    # Re-run the same commands used to generate the stored results
    new_file <- write_model(
        'test_model',
        list('BioCro:harmonic_energy'),
        list('BioCro:harmonic_oscillator'),
        list(position = 1, velocity = 0),
        list(timestep = 1, mass = 1, spring_constant = 0.5),
        BioCro::default_ode_solvers[['boost_rkck54']]
    )

    # To update: run the code above, and then type the following from an R
    # session running in this directory:
    #
    #  writeLines(new_file, stored_fname)

    expect_equal(
        new_file,
        stored_file
    )
})

test_that('Stored models can be read and run', {
    # Get stored results
    stored_fname <- file.path('..', 'test_data', 'test_model.R')

    expect_silent(
        source(stored_fname) # creates a list called `test_model`
    )

    res <- expect_silent(
        with(test_model, {BioCro::run_biocro(
            initial_values,
            parameters,
            data.frame(time = seq_len(100)),
            direct_modules,
            differential_modules,
            ode_solver
        )})
    )
})
