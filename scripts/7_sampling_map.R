### Creator: Jorge Medina
### Date: Oct/7/2024


#### This code maps collecting sites of grasshoppers of the genus Romalea

# Load and install necessary packages
# install.packages(c("dplyr", "ggplot2", "maps", "mapdata", "rnaturalearth",
#                    "rnaturalearthdata", "sf", "raster", "geodata", "terra", "viridis"))

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


# Set working directory
input <- "C:/Users/beto2/Documents/cophylogeny_project/input"

output_map <- "C:/Users/beto2/Documents/cophylogeny_project/output"

# Load database
romalea.collection <- read.csv(paste0(input, "/databases/Romalea_collecting_data.csv"))
colnames(romalea.collection)

# Shapes biogeographic regions
bioregions <- paste0(input, "/shape_files/Provinces/")

# Select working columns
romalea.collection.filtered <- romalea.collection %>% dplyr::select(host_species, latitude,
                                                                    longitude, TAMUIC_IC_CODE,
                                                                    infected)
#View(romalea.collection.filtered)

# remove empty rows
romalea.clean <- romalea.collection.filtered %>% filter(TAMUIC_IC_CODE != "")
romalea.clean <- romalea.clean %>% filter(TAMUIC_IC_CODE != "TAMUIC-IGC-004565")

#View(romalea.clean)

# remove duplicated rows
romalea.clean <- romalea.clean %>% unique()
#View(romalea.clean)

# filter only samples with gregarines
romalea.clean.infected <- romalea.clean %>% filter(infected == "y")

# filter only samples without gregarines
romalea.clean.not.infected <- romalea.clean %>% filter(infected != "y")
romalea.clean.not.infected$infected <- "n"


# # Add an infection status column to each data frame
# romalea.clean.infected$infection_status <- "Parasitized"
# romalea.clean.not.infected$infection_status <- "Not Parasitized"

# Combine the data frames
romalea.clean <- rbind(romalea.clean.infected, romalea.clean.not.infected)

romalea_ordered <- romalea.clean %>%
  arrange(infected)




# Download elevation data
usa_elevation <- elevation_30s(country = "USA", path = tempdir())
mexico_elevation <- elevation_30s(country = "MEX", path = tempdir())
costa_rica_elevation <- elevation_30s(country = "CRI", path = tempdir())
guatemala_elevation <- elevation_30s(country = "GTM", path = tempdir())
honduras_elevation <- elevation_30s(country = "HND", path = tempdir())
el_salvador_elevation <- elevation_30s(country = "SLV", path = tempdir())
nicaragua_elevation <- elevation_30s(country = "NIC", path = tempdir())
belize_elevation <- elevation_30s(country = "BLZ", path = tempdir())
panama_elevation <- elevation_30s(country = "PAN", path = tempdir())


# Merge the downloaded elevation data
merged.elevation <- mosaic(usa_elevation, mexico_elevation, costa_rica_elevation,
                           guatemala_elevation, honduras_elevation,
                           el_salvador_elevation, nicaragua_elevation,
                           belize_elevation, panama_elevation)


# Crop the merged data to the region of interest
elevation.crop <- crop(merged.elevation, ext(-120, -80, 8, 35))


# Convert the cropped elevation data to a data frame
elevation.df <- as.data.frame(elevation.crop, xy = TRUE)

# View the first few rows of the data frame
head(elevation.df)

# Download North America map data
north.america <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")



custom_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FFFF33", 
                   "#A65628", "#F781BF", "#FF7F00", "#66C2A5", "#FC8D62", 
                   "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F")


# Biogeographic regions
neotropic <- st_read(paste0(bioregions, "NeotropicMap_Geo.shp"))
neartic <- st_read(paste0(bioregions, "A_Neartica_Escalanteetal2021_Geograficas.shp"))


collecting.map <- ggplot() +
  # Elevation raster (using transparency)
  geom_raster(data = elevation.df, aes(x = x, y = y, alpha = USA_elv_msk)) +  
  scale_alpha_continuous(range = c(0.3, 1), guide = "none") +  
  
  # Map outline
  geom_sf(data = north.america, fill = NA, color = "white") +
  
  geom_sf(data = neotropic, fill = NA, color = "darkred", linewidth  = 0.1, linetype = "solid") +
  geom_sf(data = neartic, fill = NA, color = "darkred", linewidth  = 0.1, linetype = "solid") +
  
  # Jitter points (inside fill with viridis, black outline)
  geom_jitter(data = romalea.clean, 
              aes(x = longitude, y = latitude, shape = infected, fill = as.character(host_species)), 
              width = 0.3, height = 0, size = 2.8, stroke = 1, color = "black") +  
  
  # Shape customization
  scale_shape_manual(values = c("y" = 22, "n" = 21)) +  
  
  # Use viridis for inside fill of points
  scale_fill_manual(values = custom_colors) +  
  
  # Fix legend to show correct colors **and keep shape**
  guides(fill = guide_legend(override.aes = list(color = "black", shape = 22, size = 3))) +
  
  # Coordinate limits and labels
  coord_sf(xlim = c(-107, -80), ylim = c(8, 33), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude", fill = "Host species", shape = "Infection status") +
  
  # Theme settings
  theme_minimal() +
  theme(legend.position = "right")



#collecting.map
ggsave(paste0(output_map, "/figures/collecting_map.svg"), collecting.map, device = "svg",
       width = 10, height = 8, units = "in", dpi = 600)







### Sample categorization by geographic province
romalea.sf <- st_as_sf(romalea.clean, coords = c("longitude", "latitude"), crs = 4326)

# fix geometries
neotropic_valid <- st_make_valid(neotropic)


# Add information of neotropical provinces
romalea.neotropic <- st_join(romalea.sf, neotropic_valid)

# Add information of neartic provinces
romalea.na <- romalea.neotropic[romalea.neotropic$Provincias == "N/A", ]

romalea.neartic <- st_join(romalea.na, neartic)


### Combine neotropical and neartic information
romalea.neartic <- romalea.neartic %>%
  dplyr::select(-Provincias) %>%
  dplyr::rename(Provincias = NameProvin)


romalea.filled <- bind_rows(
  romalea.neotropic[romalea.neotropic$Provincias != "N/A", ],
  romalea.neartic
)


romalea.filled <- romalea.filled %>%
  dplyr::select(-Area_km, -Area_Km2, -Area_Ha)


romalea.filled <- romalea.filled %>%
  mutate(
    Region     = coalesce(Region, Region.y),
    Subregion  = coalesce(Subregion, Subregion.y),
    IDSubreg   = coalesce(IDSubreg, IDSubreg.y),
    Dominio    = coalesce(Dominio, Dominio.y),
    IDDominio  = coalesce(IDDominio, IDDominio.y),
    IDProvince = coalesce(IDProvince, IDProvince.y),
    IDProv     = coalesce(IDProv, ID_TEE2021)
  ) %>%
  # Drop the redundant columns
  dplyr::select(
    -Region.x, -Region.y,
    -Subregion.x, -Subregion.y,
    -IDSubreg.x, -IDSubreg.y,
    -Dominio.x, -Dominio.y,
    -IDDominio.x, -IDDominio.y,
    -IDProvince.x, -IDProvince.y,
    -IDProvince, -IDDominio,
    -IDSubreg, -IDProv,
    -ID_TEE2021, -Perimetro, -Cita
  )


if(!dir.exists(paste0(output_map, "/geographic_provinces"))){
  dir.create(paste0(output_map, "/geographic_provinces"))
}

write.csv(romalea.filled, paste0(output_map, "/geographic_provinces/geographic_provinces.csv"))

