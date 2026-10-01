# Load libraries
library(dplyr) # for pipe operator (%>%) and mutate
library(tidyr) # for uncount

# Clear workspace
rm(list = ls())

# Path to raw data CSV
fpath <- file.path('soyface_weather_data', 'soyFACE_weather_data_2004thru2011.csv')

# Load raw data
soyface_precip_raw <- read.csv(fpath)

# Create a list where each element is a data frame representing one year of
# precipitation measurements at hourly values
soyface_precip <- by(
    soyface_precip_raw,
    soyface_precip_raw[['Year']],
    function(x) {
        # Create an hourly version, duplicating the daily value at each hourly
        # time point
        x_hourly <- x %>%
            # Add an hourly sequence per day
            uncount(weights = 24, .id = "hour") %>%
            # Adjust the hour (0 to 23)
            mutate(hour = hour - 1, precip = precip.mm. / 24)

        # Return a subset of columns, using standard BioCro column names
        data.frame(
            year   = x_hourly[['Year']],
            doy    = x_hourly[['DOY']],
            hour   = x_hourly[['hour']],
            precip = x_hourly[['precip']]
        )
    }
)

# Save the results
save(soyface_precip, file = 'soyface_precip.rdata')
