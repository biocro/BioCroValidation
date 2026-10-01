\name{soyface_precip}

\docType{data}

\alias{soyface_precip}

\title{SoyFACE Precipitation Data}

\description{
  Eight years of precipitation data collected at the SoyFACE facility in
  Champaign, IL from 2004 - 2011.
}

\usage{soyface_precip}

\format{
  A list of eight named elements, where each element is a data frame with the
  following columns:
  \itemize{
    \item \code{year}: The year
    \item \code{doy}: The day of year
    \item \code{hour}: The hour of day
    \item \code{precip}: The precipitation rate expressed in mm / hr
  }

  Each element represents a single year of data, and the name of each element
  is the corresponding year.
}

\source{
  Daily precipitation totals were obtained from the supplemental information of
  Grey \emph{et al.} 2016 (\doi{10.1038/nplants.2016.132}), which is available
  at \doi{10.5061/dryad.g0v62}. Daily totals were converted to hourly rates by
  assuming a uniform rate across each day.

  The original data and the processing script are included with the
  BioCroValidation package; their locations can be found by typing
  \code{system.file('extdata', 'process_soyface_precip.R', package = 'BioCroValidation')}
  or
  \code{system.file('extdata', 'soyface_weather_data', package = 'BioCroValidation')}
  in an R session.
}

\keyword{datasets}
