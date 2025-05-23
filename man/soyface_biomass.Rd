\name{soyface_biomass}

\docType{data}

\alias{soyface_biomass}

\title{SoyFACE Soybean Biomass Data}

\description{
  Two years of soybean biomass data collected at the SoyFACE facility in
  Champaign, IL during the 2002 and 2005 field seasons.
}

\usage{soyface_biomass}

\format{
  A list of four named elements, where each element is a data frame with the
  following columns:
  \itemize{
    \item \code{DOY}: The day of year
    \item \code{Leaf_Mg_per_ha}: Leaf biomass per area expressed in Mg / ha
    \item \code{Stem_Mg_per_ha}: Stem biomass per area expressed in Mg / ha
    \item \code{Rep_Mg_per_ha}: Pod biomass per area expressed in Mg / ha
    \item \code{Seed_Mg_per_ha}: Seed biomass per area expressed in Mg / ha
    \item \code{Litter_Mg_per_ha}: Mass of leaf litter accumulated between
          harvests, expressed in Mg / ha
    \item \code{CumLitter_Mg_per_ha}: Cumulative leaf litter biomass per aear
          expressed in Mg / ha
  }

  The elements named \code{ambient_2002} and \code{ambient_2005} represent mean
  biomass values measured from plants grown in ambient CO2 conditions during
  2002 and 2005, respectively.

  The elements named \code{ambient_2002_std} and \code{ambient_2005_std}
  represent the standard deviation of biomass values measured from plants grown
  in ambient CO2 conditions during 2002 and 2005, respectively.
}

\source{
  The leaf, stem, and pod data was obtained from several files in from the
  \href{https://github.com/cropsinsilico/soybean-biocro}{Soybean-BioCro GitHub repository}:
  \itemize{
    \item \code{Data/SoyFACE_data/2002_ambient_biomass.csv}
    \item \code{Data/SoyFACE_data/2005_ambient_biomass.csv}
    \item \code{Data/SoyFACE_data/2002_ambient_biomass_std.csv}
    \item \code{Data/SoyFACE_data/2005_ambient_biomass_std.csv}
  }

  See that repository for more information.

  The leaf litter accumulated between harvests was obtained from the original
  sources. The cumulative leaf litter was calculated from the amount accumulated
  between harvests. The seed mass was estimated as a fraction of the total pod
  mass, using proportionality factors determined from unpublished data collected
  in Champaign, Illinois during 2021-2022.

  These data tables have not been published previously, but were used to
  parameterize Soybean-BioCro as used in He \emph{et al.} 2024
  (\doi{10.1093/insilicoplants/diae009}).
}

\keyword{datasets}
