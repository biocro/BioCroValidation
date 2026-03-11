# Create a table of lower and upper bounds

During an optimization, it is often necessary to provide lower and upper
bounds for the parameters that are being varied. Typically they are
specified as numeric vectors, which often leads to confusing code, where
the writer and reader must remember which value corresponds to each
argument; for example, "the third element of the lower bound vector is
for alphaLeaf."

The purpose of `bounds_table` is to make the process of specifying
bounds simpler and easier to follow. It is expected that this function
will be called after
[`objective_function`](https://biocro.github.io/BioCroValidation/reference/objective_function.md).

## Usage

``` r
bounds_table(
    independent_args,
    bounds_list,
    initial_ind_arg_values = NULL
  )
```

## Arguments

- independent_args:

  The same value that was passed to
  [`objective_function`](https://biocro.github.io/BioCroValidation/reference/objective_function.md).

- bounds_list:

  A list of named elements, where each element is a numeric vector of
  length 2. The names should correspond to the independent arguments,
  and the values should indicate lower and upper bounds for the
  corresponding parameter (in any order). Any "extra" bounds (that is,
  bounds that do not correspond to any independent argument) will be
  ignored.

- initial_ind_arg_values:

  A numeric vector of initial values for each of the independent
  arguments, supplied in the same order as in `independent_args`.

## Details

The main purpose of this function is to create vectors of lower and
upper bounds, which are returned as the columns of a data frame. For
each independent argument in `independent_args`, the bounds are supplied
via the `bounds_list` input. The syntax is designed so the code calling
this function is easy for a human to parse. (See example below.)

It is also (optionally) possible to provide an initial guess for each
independent argument via the `initial_ind_arg_values` argument. When
provided, these will be checked to make sure they do not lie outside the
bounds; an error will be thrown if any do lie outside the bounds. A
warning will also be thrown if any initial guesses lie on the bounds,
since this can be problematic for some optimizers, such as
[`nmkb`](https://rdrr.io/pkg/dfoptim/man/nmkb.html).

Some optimizers, such as
[`DEoptim`](https://rdrr.io/pkg/DEoptim/man/DEoptim.html), do not
require an initial guess; in this case, there is no strong need to pass
an initial guess to `bounds_table`.

## Value

A data frame with three or four columns: `arg_name`, `lower`, `upper`,
and (optionally) `initial_value`.

The `lower` and `upper` columns are the lower and upper bounds,
determined from `bounds_list`. The `arg_name` column is the argument
name, and the rows of the table are ordered as in `independent_args`.
The `initial_value` column contains initial values, if they were
provided via the `initial_ind_arg_values` argument.

## Examples

``` r
# Make a list of independent arguments; the values are not used for anything
independent_args <- list(
  alphaLeaf = 0,
  alphaRoot = 0,
  alphaStem = 0,
  betaLeaf  = 0,
  betaRoot  = 0,
  betaStem  = 0
)

# Specify bounds and initial guess for each. Note that:
#
# - The bounds will be reordered to follow the same order as the
#   `independent_args`, but the initial guess is assumed to already follow the
#   same order as the `independent_args`.
#
# - The bounds for the two extra parameters are ignored when forming the table.
#
# - The lower and upper bounds can be supplied as (upper, lower)
#   or (lower, upper) pairs.
#
b_ll <- -50 # Lower limit for beta parameters
a_ul <- 50  # Upper limit for alpha parameters

bounds <- bounds_table(
  independent_args,
  list(
    betaStem  = c(0, b_ll),
    betaRoot  = c(0, b_ll),
    betaLeaf  = c(0, b_ll),
    alphaStem = c(0, a_ul),
    alphaRoot = c(0, a_ul),
    alphaLeaf = c(0, a_ul),
    extraPar1 = c(0, 5),
    extraPar2 = c(0, 6)
  ),
  c(1, 1, 1, -1, -1, -1)
)

print(bounds)
#>    arg_name lower upper initial_value
#> 1 alphaLeaf     0    50             1
#> 2 alphaRoot     0    50             1
#> 3 alphaStem     0    50             1
#> 4  betaLeaf   -50     0            -1
#> 5  betaRoot   -50     0            -1
#> 6  betaStem   -50     0            -1

# Now the properly-ordered lower and upper limits can be accessed as follows:
bounds$lower
#> [1]   0   0   0 -50 -50 -50

bounds$lower
#> [1]   0   0   0 -50 -50 -50
```
