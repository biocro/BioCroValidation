# SoyFACE Soybean Biomass Data

Two years of soybean biomass data collected at the SoyFACE facility in
Champaign, IL during the 2002 and 2005 field seasons.

## Usage

``` r
soyface_biomass
```

## Format

A list of four named elements, where each element is a data frame with
the following columns:

- `DOY`: The day of year

- `Leaf_Mg_per_ha`: Leaf biomass per area expressed in Mg / ha

- `Stem_Mg_per_ha`: Stem biomass per area expressed in Mg / ha

- `Rep_Mg_per_ha`: Pod biomass per area expressed in Mg / ha

- `Seed_Mg_per_ha`: Seed biomass per area expressed in Mg / ha

- `Litter_Mg_per_ha`: Mass of leaf litter accumulated between harvests,
  expressed in Mg / ha

- `CumLitter_Mg_per_ha`: Cumulative leaf litter biomass per aear
  expressed in Mg / ha

The elements named `ambient_2002` and `ambient_2005` represent mean
biomass values measured from plants grown in ambient CO2 conditions
during 2002 and 2005, respectively.

The elements named `ambient_2002_std` and `ambient_2005_std` represent
the standard deviation of biomass values measured from plants grown in
ambient CO2 conditions during 2002 and 2005, respectively.

## Source

The leaf, stem, and pod data was obtained from several files in from the
[Soybean-BioCro GitHub
repository](https://github.com/cropsinsilico/soybean-biocro):

- `Data/SoyFACE_data/2002_ambient_biomass.csv`

- `Data/SoyFACE_data/2005_ambient_biomass.csv`

- `Data/SoyFACE_data/2002_ambient_biomass_std.csv`

- `Data/SoyFACE_data/2005_ambient_biomass_std.csv`

See that repository for more information.

The leaf litter accumulated between harvests was obtained from the
original sources. The cumulative leaf litter was calculated from the
amount accumulated between harvests. The seed mass was estimated as a
fraction of the total pod mass, using proportionality factors determined
from unpublished data collected in Champaign, Illinois during 2021-2022.

These data tables have not been published previously, but were used to
parameterize Soybean-BioCro as used in He *et al.* 2024
([doi:10.1093/insilicoplants/diae009](https://doi.org/10.1093/insilicoplants/diae009)
).
