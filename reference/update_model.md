# Update a BioCro model definition

Following an optimization, it is typically necessary to update the
initial values and/or parameters of a base model definition with new
values determined during the optimization. The `update_model` function
helps to streamline this process. It is expected that this function will
be called after
[`objective_function`](https://biocro.github.io/BioCroValidation/reference/objective_function.md).

## Usage

``` r
update_model(
    base_model_definition,
    independent_args,
    new_ind_arg_values,
    dependent_arg_function = NULL
  )
```

## Arguments

- base_model_definition:

  The same value that was passed to
  [`objective_function`](https://biocro.github.io/BioCroValidation/reference/objective_function.md).

- independent_args:

  The same value that was passed to
  [`objective_function`](https://biocro.github.io/BioCroValidation/reference/objective_function.md).

- new_ind_arg_values:

  A numeric vector representing new values of the independent arguments,
  typically determined by an optimizer.

- dependent_arg_function:

  The same value that was passed to
  [`objective_function`](https://biocro.github.io/BioCroValidation/reference/objective_function.md).

## Value

A list based on `base_model_definition` but with new values of some of
its `initial_values` and `parameters`, as specified by the elements of
`independent_args` and `new_ind_arg_values`.

## Examples

``` r
if (require(BioCro)) {
  # Update the default Soybean-BioCro model with new values of `Leaf` (an
  # initial value) and `alphaStem` (a parameter)
  base_model <- BioCro::soybean

  new_model <- update_model(
    base_model,
    list(Leaf = 1, alphaLeaf = 2), # The values here will not be used
    c(1000, 2000)                  # These are the actual new values
  )

  # Compare the two models
  cat('\n\nComparing initial Leaf values:\n')
  print(base_model$initial_values$Leaf)
  print(new_model$initial_values$Leaf)

  cat('\n\nComparing alphaLeaf values:\n')
  print(base_model$parameters$alphaLeaf)
  print(new_model$parameters$alphaLeaf)
}
#> 
#> 
#> Comparing initial Leaf values:
#> [1] 0.06312
#> [1] 1000
#> 
#> 
#> Comparing alphaLeaf values:
#> [1] 23.36771
#> [1] 2000
```
