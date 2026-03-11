# Changelog

## Changes in BioCroValidation Version 0.3.0 (2026-03-11)

- Fixed typos in the help page for `objective_function`, and in the
  `add_norm` function (defined in `R/objective_function_helpers.R`)

- Allowed user-supplied regularization functions

- Allowed driver-specific initial values and parameters

- Errors that occur while running simulations are now caught so they do
  not prevent an optimization from finishing

- More options for `debug_mode` are now available; the default setting
  (`minimal`) only prints info to the terminal when an issue with a
  simulation occurs

- Improved formatting of output from `write_module`, so that initial
  value and parameter lists are alphabetized, equals signs are aligned,
  and module names are preserved

- Fixed a bug where calling `objective_function` with
  `dependent_arg_function = NULL` and `verbose_mode = TRUE` caused an
  error

## Changes in BioCroValidation Version 0.2.0 (2025-05-23)

- Added 2002 and 2005 SoyFACE biomass and standard deviation data.

- Added several new functions: `objective_function`, `update_model`, and
  `bounds_table`.

- Added two new vignettes: a “Getting Started” article
  (`BioCroValidation.Rmd`) and a user guide illustrating how to perform
  a model parameterization (`parameterizing_soybean_biocro.Rmd`).

## Changes in BioCroValidation Version 0.1.0

- This is the first version of BioCroValidation. At this point, the
  package is in a state of rapid development, and not all changes will
  be described here.

- We are reserving version `1.0.0` for a more stable and complete future
  release; until then, major changes should only increase the minor
  version number.
