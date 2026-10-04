## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  fig.width = 10, # Set default plot width (adjust as needed)
  fig.height = 8, # Set default plot height (adjust as needed)
  fig.align = "center" # Center align all plots
)

## -----------------------------------------------------------------------------
library(gmwmx2)

## -----------------------------------------------------------------------------
all_stations <- download_all_stations_ngl()
head(all_stations)

## -----------------------------------------------------------------------------
data_CERN <- download_station_ngl("CERN")

## -----------------------------------------------------------------------------
attributes(data_CERN)
head(data_CERN$df_position)

## ----station-map--------------------------------------------------------------
if (!is.data.frame(data_CERN$df_position) || nrow(data_CERN$df_position) == 0L) {
  message("Skipping map: no position data were downloaded for CERN.")
} else if (!requireNamespace("leaflet", quietly = TRUE)) {
  message("Install the optional 'leaflet' package to display the station map.")
} else {
  latitude <- data_CERN$df_position$nominal_station_latitude
  longitude <- data_CERN$df_position$nominal_station_longitude
  valid <- which(is.finite(latitude) & abs(latitude) <= 90 &
                   is.finite(longitude))

  if (length(valid) == 0L) {
    message("Skipping map: no valid station coordinates are available for CERN.")
  } else {
    latitude <- latitude[valid[1L]]
    longitude <- ((longitude[valid[1L]] + 180) %% 360) - 180

    station_map <- leaflet::leaflet()
    station_map <- leaflet::addTiles(station_map)
    station_map <- leaflet::setView(station_map, lng = longitude,
                                    lat = latitude, zoom = 12)
    leaflet::addMarkers(station_map, lng = longitude, lat = latitude,
                        popup = "CERN GNSS station", label = "CERN")
  }
}

## -----------------------------------------------------------------------------
head(data_CERN$df_equipment_software_changes)

## -----------------------------------------------------------------------------
head(data_CERN$df_earthquakes)

## -----------------------------------------------------------------------------
if (is.data.frame(data_CERN$df_position) && nrow(data_CERN$df_position) > 0L) {
  plot(data_CERN)
  plot(data_CERN, component = "N")
  plot(data_CERN, component = "E")
  plot(data_CERN, component = "V")
} else {
  message("Skipping plots: no position data were downloaded for CERN. Check the download warnings and try again later.")
}

