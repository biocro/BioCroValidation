all_rmse <- function(
    base_model_definition,
    data_driver_pairs,
    data_definitions = list(),
    post_process_function = NULL,
    verbose_startup = TRUE
)
{
    # Make sure all weights in the data-driver pairs are set to 1
    for (i in seq_along(data_driver_pairs)) {
        data_driver_pairs[[i]]$weight <- 1
    }

    # objective_function requires one independent argument, so we will take one
    # from the model
    independent_args <- base_model_definition$parameters[1]

    # Create an objective function with all weights set to 1, no regularization,
    # and no penalty term. Using the `mean` method for normalization will enable
    # us to calculate RMSE later.
    obj_fun <- objective_function(
        base_model_definition,
        data_driver_pairs,
        independent_args,
        list(), # quantity_weights
        data_definitions,
        normalization_method = 'mean',
        stdev_weight_method = 'equal',
        post_process_function = post_process_function,
        verbose_startup = verbose_startup
    )

    # Now call the objective function to get all the error terms
    error_terms <- obj_fun(
        as.numeric(independent_args),
        return_terms = TRUE
    )

    # We only want the terms from the data-driver pairs
    terms_from_data_driver_pairs <- error_terms$terms_from_data_driver_pairs

    # Because the `mean` normalization method was used, each term is actually
    # the mean squared error divided by the number of data-driver pairs
    ndrivers <- length(data_driver_pairs)

    for (i in seq_along(terms_from_data_driver_pairs)) {
        terms_from_data_driver_pairs[[i]]$RMSE_terms <-
            lapply(terms_from_data_driver_pairs[[i]]$least_squares_terms, function(x) {sqrt(x * ndrivers)})
    }

    # Convert to a list of data frames
    rmse_df_list <- lapply(terms_from_data_driver_pairs, function(x) {
        as.data.frame(x$RMSE_terms)
    })

    # Find all unique column names across all data frames
    all_cols <- unique(unlist(lapply(rmse_df_list, colnames)))

    # Map over each data frame to add missing columns filled with NA
    padded_dfs <- lapply(rmse_df_list, function(df) {
        missing_cols <- setdiff(all_cols, colnames(df))
        df[missing_cols] <- NA

        # Ensure the column order is identical
        return(df[, all_cols, drop = FALSE])
    })

    # Combine into a single data frame and return
    rmse_table <- do.call(rbind, padded_dfs)

    rmse_table$data_driver_pair <- rownames(rmse_table)

    rownames(rmse_table) <- NULL

    rmse_table
}
