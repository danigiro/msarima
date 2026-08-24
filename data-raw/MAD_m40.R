## Hourly traffic volume on the M-40 ring road, Madrid (MAD_m40)
##
## Raw file: datasetcompleto_3753_0205.csv, half-hourly vehicle counts recorded by
## sensor 3753 between 2024-02-01 00:00 and 2024-05-31 23:30 (5808 observations).
## The two half-hours of each clock hour are summed to obtain hourly volumes.

raw <- read.csv(
  "data-raw/datasetcompleto_3753_0205.csv",
  header = TRUE,
  dec = ".",
  stringsAsFactors = FALSE
)
raw$fecha <- as.POSIXct(raw$fecha, tz = "UTC")
raw <- raw[order(raw$fecha), ]

## the series must be a complete half-hourly grid before aggregating
stopifnot(
  length(unique(raw$id)) == 1,
  !anyNA(raw$volume),
  !anyDuplicated(raw$fecha),
  all(diff(as.numeric(raw$fecha)) == 30 * 60),
  as.integer(format(raw$fecha[1], "%M")) == 0,
  nrow(raw) %% 2 == 0
)

## sum the two half-hours of each hour
hourly <- colSums(matrix(raw$volume, nrow = 2))
stopifnot(length(hourly) == nrow(raw) / 2)

## 2024-02-01 is a Thursday, so observation 1 is the first hour of a Thursday
MAD_m40 <- forecast::msts(hourly, seasonal.periods = c(24, 168))

usethis::use_data(MAD_m40, overwrite = TRUE)
