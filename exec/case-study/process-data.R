devtools::load_all()
library(dplyr)
library(sf)
library(lubridate)
library(raster)
library(stars)

# Load radar data
radar_file = file.path(downloads_dir(), "radar.nc")
data = raster::brick(radar_file)

# Get the coordinates of the data
coords = rasterToPoints(data[[1]], spatial = TRUE) |>
  st_as_sf() |>
  {\(x) cbind(x, st_coordinates(x))}() |>
  dplyr::select(geometry)

# Transform coordinates from [m] to [km] to remove unneccesary zeros
proj = st_crs(coords)[["input"]] |>
  sub(pattern = "units=m", replacement = "units=km") |>
  st_crs()
coords = st_transform(coords, proj)

# Add height information
height_raster = raster::raster(file.path(downloads_dir(), "dem.tif"))
transformed_coords = st_coordinates(st_transform(coords, st_crs(height_raster)))
coords$height = raster::extract(height_raster, transformed_coords)
coords$height = ifelse(is.na(coords$height), 0, coords$height)

# Get the times of the data
times = colnames(data[1, 1]) |>
  sub(pattern = "X", replacement = "") |>
  lubridate::ymd_hms()

# Extract the data into a matrix
data_matrix = raster::as.matrix(data)
colnames(data_matrix) = NULL
data_matrix = t(data_matrix)

mydata = list(
  data = data_matrix,
  coords = coords,
  times = times)
saveRDS(mydata, file.path(downloads_dir(), "radar.rds"))


# ======================================================================
# Plot for the paper
# Be aware that we need to execute the final part of 1-download-data.R,
# where we download DEM data, before we execute the code below
# ======================================================================
radar = readRDS(file.path(downloads_dir(), "radar.rds"))
coords = radar$coords

height_raster = stars::read_stars(file.path(downloads_dir(), "dem.tif"))

# Create a bbox for the area we wish to plot
bbox = (st_bbox(coords) + .8 * c(-50, -50, 50, 50)) |>
  st_as_sfc() |>
  st_transform(st_crs(height_raster))

# Extract all altitudes inside the bbox,
# and transform to lon/lat
tmp = st_crop(height_raster, bbox) |>
  st_transform(crs = st_crs(4326))

# Create a box showing the study area for the case study
bbox2 = st_bbox(coords) |>
  st_as_sfc() |>
  st_transform(st_crs(tmp))

# This is the location of the Rissa radar
rissa = st_point(c(10.203845, 63.690527)) |>
  st_sfc(crs = 4326)

jpeg(file.path(image_dir(), "height-map.jpg"))
plot(tmp, axes = TRUE, reset = FALSE, nbreaks = 25, main = "")
plot(bbox2, lwd = 5, add = TRUE)
plot(rissa, add = TRUE, pch = 17, cex = 2, col = "red")
dev.off()
