#############################################################
# Title: Mapping Romalea Grasshopper Collecting Sites and Geographic Provinces
# Creator: Jorge Medina-Duran
# Date: Oct/7/2024
#
# Description:
# This script maps collecting sites of Romalea grasshoppers and assigns samples to biogeographic provinces.
# It performs:
# - Loading and cleaning of collecting data
# - Plotting sampling localities over North America with elevation background
# - Assigning samples to Neotropical and Nearctic provinces
# - Saving province-level assignments
#
# Output:
# - Collecting map (SVG file)
# - Geographic province assignment table (CSV file)
#
# Notes:
# - Change only the `input` and `output_map` paths for reproducibility.
# - Figures and results saved to `/figures/` and `/geographic_provinces/` subdirectories.
#############################################################

# Load libraries
library(dplyr)
library(ggplot2)
library(maps)
library(mapdata)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(raster)
library(geodata)
library(terra)
library(viridis)
library(RColorBrewer)

#### Define paths ----
input <- "YOUR_PATH/input"
output_map <- "YOUR_PATH/output"

#### Load and clean Romalea collecting database ----
romalea.collection <- read.csv(paste0(input, "/databases/Romalea_collecting_data.csv"))

romalea.clean <- romalea.collection %>%
  dplyr::select(host_species, latitude, longitude, TAMUIC_IC_CODE, infected) %>%
  filter(TAMUIC_IC_CODE != "", TAMUIC_IC_CODE != "TAMUIC-IGC-004565") %>%
  unique()

# Categorize by infection status
romalea.clean.infected <- romalea.clean %>% filter(infected == "y")
romalea.clean.not.infected <- romalea.clean %>% filter(infected != "y") %>%
  mutate(infected = "n")

romalea.clean <- bind_rows(romalea.clean.infected, romalea.clean.not.infected) %>%
  arrange(infected)

#### Download and merge elevation data ----
usa_elevation <- elevation_30s(country = "USA", path = tempdir())
mexico_elevation <- elevation_30s(country = "MEX", path = tempdir())
costa_rica_elevation <- elevation_30s(country = "CRI", path = tempdir())
guatemala_elevation <- elevation_30s(country = "GTM", path = tempdir())
honduras_elevation <- elevation_30s(country = "HND", path = tempdir())
el_salvador_elevation <- elevation_30s(country = "SLV", path = tempdir())
nicaragua_elevation <- elevation_30s(country = "NIC", path = tempdir())
belize_elevation <- elevation_30s(country = "BLZ", path = tempdir())
panama_elevation <- elevation_30s(country = "PAN", path = tempdir())

merged.elevation <- mosaic(usa_elevation, mexico_elevation, costa_rica_elevation,
                           guatemala_elevation, honduras_elevation, el_salvador_elevation,
                           nicaragua_elevation, belize_elevation, panama_elevation)

elevation.crop <- crop(merged.elevation, ext(-120, -80, 8, 35))
elevation.df <- as.data.frame(elevation.crop, xy = TRUE)

#### Load base maps ----
north.america <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
neotropic <- st_read(paste0(input, "/shape_files/Provinces/NeotropicMap_Geo.shp"))
neartic <- st_read(paste0(input, "/shape_files/Provinces/A_Neartica_Escalanteetal2021_Geograficas.shp"))

#### Plot collecting sites ----
custom_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FFFF33",
                   "#A65628", "#F781BF", "#FF7F00", "#66C2A5", "#FC8D62",
                   "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F")

collecting.map <- ggplot() +
  geom_raster(data = elevation.df, aes(x = x, y = y, alpha = USA_elv_msk)) +
  scale_alpha_continuous(range = c(0.3, 1), guide = "none") +
  geom_sf(data = north.america, fill = NA, color = "white") +
  geom_sf(data = neotropic, fill = NA, color = "darkred", linewidth = 0.1) +
  geom_sf(data = neartic, fill = NA, color = "darkred", linewidth = 0.1) +
  geom_jitter(data = romalea.clean,
              aes(x = longitude, y = latitude, shape = infected, fill = as.character(host_species)),
              width = 0.3, height = 0, size = 2.8, stroke = 1, color = "black") +
  scale_shape_manual(values = c("y" = 22, "n" = 21)) +
  scale_fill_manual(values = custom_colors) +
  guides(fill = guide_legend(override.aes = list(color = "black", shape = 22, size = 3))) +
  coord_sf(xlim = c(-107, -80), ylim = c(8, 33), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude", fill = "Host species", shape = "Infection status") +
  theme_minimal() +
  theme(legend.position = "right")

dir.create(paste0(output_map, "/figures"), recursive = TRUE, showWarnings = FALSE)

ggsave(paste0(output_map, "/figures/collecting_map.svg"), collecting.map, device = "svg",
       width = 10, height = 8, units = "in", dpi = 600)

#### Assign collecting sites to geographic provinces ----
romalea.sf <- st_as_sf(romalea.clean, coords = c("longitude", "latitude"), crs = 4326)
neotropic_valid <- st_make_valid(neotropic)

# Join Neotropical provinces
romalea.neotropic <- st_join(romalea.sf, neotropic_valid)

# Join Nearctic provinces for N/A samples
romalea.na <- romalea.neotropic[romalea.neotropic$Provincias == "N/A", ]
romalea.neartic <- st_join(romalea.na, neartic)

# Combine filled data
romalea.neartic <- romalea.neartic %>%
  select(-Provincias) %>%
  rename(Provincias = NameProvin)

romalea.filled <- bind_rows(
  romalea.neotropic[romalea.neotropic$Provincias != "N/A", ],
  romalea.neartic
) %>%
  select(-c(Area_km, Area_Km2, Area_Ha)) %>%
  mutate(
    Region = coalesce(Region, Region.y),
    Subregion = coalesce(Subregion, Subregion.y),
    IDSubreg = coalesce(IDSubreg, IDSubreg.y),
    Dominio = coalesce(Dominio, Dominio.y),
    IDDominio = coalesce(IDDominio, IDDominio.y)
  ) %>%
  select(-c(Region.x, Region.y, Subregion.x, Subregion.y, IDSubreg.x, IDSubreg.y,
            Dominio.x, Dominio.y, IDDominio.x, IDDominio.y, IDProvince, IDProv,
            ID_TEE2021, Perimetro, Cita))

# Save output
dir.create(paste0(output_map, "/geographic_provinces"), recursive = TRUE, showWarnings = FALSE)
write.csv(romalea.filled, paste0(output_map, "/geographic_provinces/geographic_provinces.csv"))
