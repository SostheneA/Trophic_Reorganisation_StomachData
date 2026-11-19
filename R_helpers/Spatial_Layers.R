# Spatial ecoregion layers

library(sf)
library(dplyr)

GASPE_TO <- "Central"

SPATIAL_DIR <- "data/Spatial_data"

sf_use_s2(FALSE)

EAR_sf     <- st_make_valid(st_read(file.path(SPATIAL_DIR, "EAR_map/EAR_map.shp"),   quiet = TRUE))
NAFO_4T_sf <- st_make_valid(st_read(file.path(SPATIAL_DIR, "NAFO/nafo_2014_02.shp"), quiet = TRUE))

z_4t <- NAFO_4T_sf %>% filter(grepl("^4T", level_2)) %>% st_union()

zone_sf <- function(geom_a, geom_b, name, col) {
  g <- st_union(st_geometry(suppressWarnings(st_intersection(geom_a, geom_b))))
  out <- st_sf(x = name, geometry = g)
  names(out)[1] <- col
  out
}

z_4tmn  <- NAFO_4T_sf %>% filter(level_2 %in% c("4TM", "4TN")) %>% st_union()
z_4thlg <- NAFO_4T_sf %>% filter(level_2 %in% c("4TH", "4TL", "4TG")) %>% st_union()

eco_chaleurs <- zone_sf(EAR_sf %>% filter(PID == 50), z_4tmn,  "Chaleur Bay", "EcoZone")
eco_northumb <- zone_sf(EAR_sf %>% filter(PID == 6),  z_4thlg, "Northumberland",    "EcoZone")

nafo_4tgj <- NAFO_4T_sf %>% filter(level_2 %in% c("4TG", "4TJ")) %>% st_union()
eco_4tgj  <- st_sf(EcoZone = "4TGJ clean",
                   geometry = st_union(suppressWarnings(
                     st_difference(nafo_4tgj, st_geometry(eco_northumb)))))

z_nafo_otl <- NAFO_4T_sf %>% filter(level_2 %in% c("4TO", "4TN", "4TL")) %>% st_union()
z_ear_mag  <- EAR_sf %>% filter(regn_nm %in% c("centre", "magdalen_shallows")) %>% st_union()
eco_shediac <- zone_sf(z_ear_mag, z_nafo_otl, "Shediac-Magdalen", "EcoZone") %>%
  st_difference(st_geometry(eco_chaleurs)) %>% st_as_sf()

eco_pure <- NAFO_4T_sf %>%
  filter(level_2 %in% c("4TF", "4TK")) %>%
  mutate(EcoZone = paste(level_2, "Only")) %>%
  select(EcoZone)

new_ecoregions_clean <- bind_rows(eco_chaleurs, eco_northumb, eco_4tgj,
                                  eco_shediac, eco_pure) %>%
  filter(!st_is_empty(.)) %>%
  st_make_valid() %>%
  st_collection_extract("POLYGON")

st_write(new_ecoregions_clean, file.path(SPATIAL_DIR, "new_ecoregions_final.shp"),
         delete_layer = TRUE, quiet = TRUE)
cat("Written:", file.path(SPATIAL_DIR, "new_ecoregions_final.shp"), "\n")

box <- st_sfc(st_polygon(list(matrix(c(-64, 45.5,  -61, 45.5,  -61, 46.47,
                                       -64, 46.47, -64, 45.5),
                                     ncol = 2, byrow = TRUE))), crs = 4326)

v2_chaleurs <- zone_sf(EAR_sf %>% filter(PID == 50), z_4tmn, "Chaleur Bay", "Zone")

z_4thg <- NAFO_4T_sf %>% filter(level_2 %in% c("4TH", "4TG")) %>% st_union()
v2_northcoast <- zone_sf(st_as_sf(z_4thg), box, "North Coast", "Zone")

z_4tjg <- NAFO_4T_sf %>% filter(level_2 %in% c("4TJ", "4TG")) %>% st_union()
v2_4tgj <- st_sf(Zone = "4TGJ coast",
                 geometry = st_union(suppressWarnings(st_difference(z_4tjg, box))))

z_ear_mag2 <- EAR_sf %>%
  filter(regn_nm %in% c("centre", "magdalen_shallows", "northumberland_strait")) %>%
  st_union()
v2_shediac <- zone_sf(z_ear_mag2, z_nafo_otl, "Shediac-Magdalen", "Zone") %>%
  st_difference(st_geometry(v2_chaleurs)) %>% st_as_sf()

v2_pure <- NAFO_4T_sf %>%
  filter(level_2 %in% c("4TF", "4TK")) %>%
  mutate(Zone = paste(level_2, "Only")) %>%
  select(Zone)

new_regions_v2_clean <- bind_rows(v2_chaleurs, v2_northcoast, v2_4tgj,
                                  v2_shediac, v2_pure) %>%
  filter(!st_is_empty(.)) %>%
  st_make_valid() %>%
  st_collection_extract("POLYGON")

st_write(new_regions_v2_clean, file.path(SPATIAL_DIR, "new_regions_v2_final.shp"),
         delete_layer = TRUE, quiet = TRUE)
cat("Written:", file.path(SPATIAL_DIR, "new_regions_v2_final.shp"), "\n")

clip_ear <- function(pid, name) zone_sf(EAR_sf %>% filter(PID == pid), z_4t, name, "Area")

area_chaleurs  <- clip_ear(50, "Chaleur Bay")
area_northumb  <- clip_ear(6,  "Northumberland")
area_central   <- clip_ear(3,  "Central")
area_gaspe     <- clip_ear(1,  GASPE_TO)

area_magdalen  <- clip_ear(5, "Magdalen Shallows") %>%
  st_difference(st_geometry(area_chaleurs)) %>%
  st_difference(st_geometry(area_northumb)) %>%
  st_as_sf()

new_gulf_clean <- bind_rows(area_chaleurs, area_northumb, area_central,
                            area_magdalen, area_gaspe) %>%
  filter(!st_is_empty(.)) %>%
  group_by(Area) %>% summarise(.groups = "drop") %>%
  st_make_valid() %>%
  st_collection_extract("POLYGON")

cat("\nAreas (km2):\n")
print(data.frame(Area = new_gulf_clean$Area,
                 km2 = round(as.numeric(st_area(new_gulf_clean)) / 1e6)))

for (i in seq_len(nrow(new_gulf_clean) - 1)) {
  for (j in (i + 1):nrow(new_gulf_clean)) {
    ov <- suppressWarnings(st_intersection(st_geometry(new_gulf_clean)[i],
                                           st_geometry(new_gulf_clean)[j]))
    a <- if (length(ov)) sum(as.numeric(st_area(ov))) / 1e6 else 0
    if (a > 1) stop(sprintf("Overlap %s x %s: %.0f km2",
                            new_gulf_clean$Area[i], new_gulf_clean$Area[j], a))
  }
}

out_each <- vapply(seq_len(nrow(new_gulf_clean)), function(i) {
  d <- suppressWarnings(st_difference(st_geometry(new_gulf_clean)[i], z_4t))
  if (length(d)) sum(as.numeric(st_area(d))) / 1e6 else 0
}, numeric(1))
stopifnot(all(out_each < 200))

gap <- suppressWarnings(st_difference(z_4t, st_union(new_gulf_clean)))
gap_km2 <- if (length(gap)) sum(as.numeric(st_area(gap))) / 1e6 else 0
cat(sprintf("4T not covered: %.0f km2 (upper estuary only is allowed)\n", gap_km2))
if (gap_km2 > 0) {
  gap_east <- suppressWarnings(st_crop(gap, xmin = -66.5, xmax = -55,
                                       ymin = 40, ymax = 55))
  gap_east_km2 <- if (length(gap_east)) sum(as.numeric(st_area(gap_east))) / 1e6 else 0
  cat(sprintf("  of which east of -66.5 (surveyed sGSL): %.0f km2\n", gap_east_km2))
  stopifnot(gap_east_km2 < 200)
}

st_write(new_gulf_clean, file.path(SPATIAL_DIR, "new_gulf_final.shp"),
         delete_layer = TRUE, quiet = TRUE)
cat("Written:", file.path(SPATIAL_DIR, "new_gulf_final.shp"), "\n")

sf_use_s2(TRUE)
