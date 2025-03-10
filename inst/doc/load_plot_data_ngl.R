## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  fig.width = 10, # Set default plot width (adjust as needed)
  fig.height = 8, # Set default plot height (adjust as needed)
  fig.align = "center" # Center align all plots
)


# knitr::opts_chunk$set(eval = FALSE)

## -----------------------------------------------------------------------------
library(gmwmx2)

## -----------------------------------------------------------------------------
all_stations <- download_all_stations_ngl()
head(all_stations)

## -----------------------------------------------------------------------------
data_1LSU <- download_station_ngl("1LSU")

## -----------------------------------------------------------------------------
attributes(data_1LSU)
head(data_1LSU$df_position)

## -----------------------------------------------------------------------------
head(data_1LSU$df_equipment_software_changes)

## -----------------------------------------------------------------------------
head(data_1LSU$df_earthquakes)

## -----------------------------------------------------------------------------
plot(data_1LSU)
plot(data_1LSU, component = "N")
plot(data_1LSU, component = "E")
plot(data_1LSU, component = "V")

