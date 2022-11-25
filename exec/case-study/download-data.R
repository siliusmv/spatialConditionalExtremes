devtools::load_all()
library(stars)

# ==============================================================================
# Download radar data.
# This script requires that you have installed the CDO program on your computer.
# See https://code.mpimet.mpg.de/projects/cdo for info and installation of the program.
# ==============================================================================

# Find all valid urls containing radar data
dates = seq(as.Date("2010-01-01"), as.Date("2021-12-31"), by = "1 month")
dates = dates[lubridate::month(dates) %in% c(6, 7, 8)]
urls = unlist(lapply(dates, get_radar_url))

# Arguments to CDO for downloading the data
args = c(
  # This selects an area of size 31x31 close to the radar
  "-selindexbox,561,591,1015,1045",
  # This is the name of the variable that contains hourly mean precipitation estimates
  "-selname,lwe_precipitation_rate")

filename = file.path(downloads_dir(), "radar.nc")

# This will take some time, and may need to be restarted several times if you have
# a bad internet connection
download(urls, args, filename)

# ============================================================================================
# Download DEM
# ============================================================================================
filename = file.path(downloads_dir(), "dem.tif")
url = "https://hoydedata.no/LaserInnsyn/Home/DownloadFile/56"
zipfile = file.path(downloads_dir(), "dem.zip")
download.file(url, zipfile)
dem_dir = file.path(downloads_dir(), "dem")
unzip(zipfile, exdir = dem_dir)
file.remove(zipfile)

files = list.files(dem_dir, full.names = TRUE)
tif_files = grep("*.tif$", files, value = TRUE)

# Combine all the tif_files into one big raster file
stars::st_mosaic(tif_files, dst = filename)
