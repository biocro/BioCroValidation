# SoyFACE Precipitation Data

Eight years of precipitation data collected at the SoyFACE facility in
Champaign, IL from 2004 - 2011.

## Usage

``` r
soyface_precip
```

## Format

A list of eight named elements, where each element is a data frame with
the following columns:

- `year`: The year

- `doy`: The day of year

- `hour`: The hour of day

- `precip`: The precipitation rate expressed in mm / hr

Each element represents a single year of data, and the name of each
element is the corresponding year.

## Source

Daily precipitation totals were obtained from the supplemental
information of Grey *et al.* 2016
([doi:10.1038/nplants.2016.132](https://doi.org/10.1038/nplants.2016.132)
), which is available at
[doi:10.5061/dryad.g0v62](https://doi.org/10.5061/dryad.g0v62) . Daily
totals were converted to hourly rates by assuming a uniform rate across
each day.

The original data and the processing script are included with the
BioCroValidation package; their locations can be found by typing
`system.file('extdata', 'process_soyface_precip.R', package = 'BioCroValidation')`
or
`system.file('extdata', 'soyface_weather_data', package = 'BioCroValidation')`
in an R session.
