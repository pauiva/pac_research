####Calculating spatial indicators per individual

# Packages
library(sf)
library(dplyr)
library(lwgeom)  # for st_make_valid
library(tidyr)


# ==== Inputs ====
roads <- st_read("netascore_Stockholm.gpkg", layer = "edge")
roads <- st_transform(roads, 3006)  # SWEREF99 TM

roads <- roads %>%
  mutate(
    index_bike_ft = as.numeric(index_bike_ft),
    index_bike_tf = as.numeric(index_bike_tf),
    bike_index = rowMeans(pick(index_bike_ft, index_bike_tf), na.rm = TRUE)
  )

roads <- roads %>%
  mutate(
    index_cat = ifelse(bike_index >= 0.6, "High", "Low")
  )

# ==== EDIT PATHS ====
neigh_shp   <- "Stadsdel.shp"  # polygons with field 'namn'
target_crs  <- 3006                          # SWEREF99 TM (meters)


if (is.na(st_crs(roads))) st_crs(roads) <- 4326
roads <- st_transform(roads, target_crs)

# Read & prep neighbourhoods
neigh <- st_read(neigh_shp, quiet = TRUE)
if (is.na(st_crs(neigh))) st_crs(neigh) <- 4326
neigh <- st_transform(neigh, target_crs)

roads <- st_crop(roads, st_bbox(neigh))

# ==== INTERSECTION: create per-neighbourhood road segments ====
# Keep only needed road attrs
roads_keep <- roads %>% select(index_cat)

roads_in_neigh <- st_intersection(
  roads_keep,
  neigh %>% select(NAMN)        # <-- neighbourhood name field
) %>%
  st_collection_extract("LINESTRING")

# ==== SEGMENT LENGTHS ====
roads_in_neigh <- roads_in_neigh %>%
  mutate(
    seg_len_m  = as.numeric(st_length(.)),
    seg_len_km = seg_len_m / 1000
  )

# ==== AGGREGATE: length by neighbourhood × category ====
by_nbhd_cat <- roads_in_neigh %>%
  st_drop_geometry() %>%
  group_by(NAMN, index_cat) %>%
  summarise(total_len_km = sum(seg_len_km, na.rm = TRUE), .groups = "drop")


# Wide table + totals + ratio
by_nbhd_wide <- by_nbhd_cat %>%
  pivot_wider(names_from = index_cat, values_from = total_len_km, values_fill = 0) %>%
  mutate(
    total_km       = coalesce(High, 0) + coalesce(Low, 0),
    high_low_ratio = ifelse(coalesce(Low, 0) > 0, High / Low, NA_real_)
  ) %>%
  arrange(NAMN)

neigh_with_stats <- neigh %>%
  left_join(by_nbhd_wide, by = "NAMN")


# ==== SAVE OUTPUTS ====
st_write(roads_in_neigh, "road_segments_by_neighbourhood.gpkg",
         layer = "segments", delete_dsn = TRUE, quiet = TRUE)

st_write(neigh_with_stats, "neighbourhoods_with_high_low_ratio.gpkg",
         layer = "neigh_stats", delete_dsn = TRUE, quiet = TRUE)

write.csv(by_nbhd_cat,  "length_by_neighbourhood_category_long.csv", row.names = FALSE)
write.csv(by_nbhd_wide, "length_by_neighbourhood_category_wide.csv", row.names = FALSE)

cat("Done. Wrote:\n",
    "- road_segments_by_neighbourhood.gpkg (segments)\n",
    "- neighbourhoods_with_high_low_ratio.gpkg (neigh_stats)\n",
    "- length_by_neighbourhood_category_long.csv / _wide.csv\n")




# --- EDIT THESE ---
neigh_with_stats <- st_read("neighbourhoods_with_high_low_ratio.gpkg",
                            layer = "neigh_stats", quiet = TRUE)

individuals <- read.csv("data_pac_12_08_2025.csv")

# Rows with coordinates
pts <- individuals %>%
  filter(!is.na(lon) & !is.na(lat)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# Rows without coordinates
no_coords <- individuals %>%
  filter(is.na(lon) | is.na(lat))

# Combine: sf for valid, plain df for missing
pts_wgs84 <- bind_rows(pts, no_coords)

# 4) Transform both to same CRS
target_crs <- 3006
if (is.na(st_crs(neigh_with_stats))) st_crs(neigh_with_stats) <- 4326
neigh_with_stats <- st_transform(neigh_with_stats, target_crs) %>% st_make_valid()
pts <- st_transform(pts_wgs84, target_crs)

# 5) Spatial join (point-in-polygon)
people_with_stats <- st_join(
  pts,
  neigh_with_stats %>% select(NAMN, High, Low, total_km, high_low_ratio),
  join = st_within,
  left = TRUE
)


# 6) Save results
st_write(people_with_stats, "individuals_with_neigh_stats.gpkg",
         layer = "people_with_stats", delete_dsn = TRUE, quiet = TRUE)
write.csv(st_drop_geometry(people_with_stats),
          "individuals_with_neigh_stats.csv", row.names = FALSE)

cat("✅ Done! Saved:\n- individuals_with_neigh_stats.gpkg\n- individuals_with_neigh_stats.csv\n")




####POIS per neighbourhood
###########################
###########################
##############################

# install.packages(c("osmdata","sf","dplyr","tidyr","lwgeom")) # if needed
library(osmdata)
library(sf)
library(dplyr)
library(tidyr)
library(lwgeom)

# --- Settings ---
target_crs <- 3006      # SWEREF99 TM (meters)
neigh_name <- "NAMN"    # neighbourhood name field in your polygons

# --- 1) Load neighbourhoods ---
neigh <- st_read("neighbourhoods_with_high_low_ratio.gpkg",
                 layer = "neigh_stats", quiet = TRUE)
if (is.na(st_crs(neigh))) st_crs(neigh) <- 4326
neigh_m   <- st_transform(neigh, target_crs) %>% st_make_valid()
neigh_wgs <- st_transform(neigh_m, 4326)  # for OSM bbox
bb <- st_bbox(neigh_wgs)

# --- 2) Run separate OSM queries for the 4 umbrella keys (ALL subtypes) ---
q_shop    <- opq(bb) |> add_osm_feature(key = "shop")
q_amenity <- opq(bb) |> add_osm_feature(key = "amenity")
q_leisure <- opq(bb) |> add_osm_feature(key = "leisure")
q_office  <- opq(bb) |> add_osm_feature(key = "office")

osm_s <- osmdata_sf(q_shop)
osm_a <- osmdata_sf(q_amenity)
osm_l <- osmdata_sf(q_leisure)
osm_o <- osmdata_sf(q_office)

# --- 3) Helper: convert an osmdata list to POINTS in target CRS ---
osm_to_points <- function(osm, target_crs) {
  out <- list()
  if (!is.null(osm$osm_points)        && nrow(osm$osm_points)        > 0) out$pt  <- st_transform(osm$osm_points, target_crs)
  if (!is.null(osm$osm_polygons)      && nrow(osm$osm_polygons)      > 0) out$pg  <- st_point_on_surface(st_transform(osm$osm_polygons, target_crs))
  if (!is.null(osm$osm_multipolygons) && nrow(osm$osm_multipolygons) > 0) out$mpg <- st_point_on_surface(st_transform(osm$osm_multipolygons, target_crs))
  if (!length(out)) return(NULL)
  # standardize columns before binding
  all_cols <- unique(unlist(lapply(out, names)))
  out2 <- lapply(out, function(x){
    miss <- setdiff(all_cols, names(x))
    for (col in miss) x[[col]] <- NA
    x[, all_cols]
  })
  do.call(rbind, out2)
}

pts_s <- osm_to_points(osm_s, target_crs)
pts_a <- osm_to_points(osm_a, target_crs)
pts_l <- osm_to_points(osm_l, target_crs)
pts_o <- osm_to_points(osm_o, target_crs)

# --- 4) Bind everything, standardize cols, keep tidy attrs, dedupe by osm_id ---
poi_list <- Filter(Negate(is.null), list(pts_s, pts_a, pts_l, pts_o))
stopifnot(length(poi_list) > 0)

all_cols <- unique(unlist(lapply(poi_list, names)))
poi_list2 <- lapply(poi_list, function(x){
  miss <- setdiff(all_cols, names(x))
  for (col in miss) x[[col]] <- NA
  x[, all_cols]
})
pois <- do.call(rbind, poi_list2)

keep_cols <- intersect(c("osm_id","shop","amenity","leisure","office","name"), names(pois))
pois <- pois %>% select(any_of(keep_cols), geometry)
if ("osm_id" %in% names(pois)) {
  pois <- pois %>% distinct(osm_id, .keep_all = TRUE)
}

# --- 5) Label by umbrella key (priority if multiple tags present on same feature) ---
pois_labeled <- pois %>%
  mutate(
    poi_key = case_when(
      !is.na(shop)    & shop    != "" ~ "shop",
      !is.na(amenity) & amenity != "" ~ "amenity",
      !is.na(leisure) & leisure != "" ~ "leisure",
      !is.na(office)  & office  != "" ~ "office",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(poi_key))

# --- 6) Collapse to umbrella groups (shops, amenity, leisure, office) ---
pois_grouped <- pois_labeled %>%
  transmute(
    osm_id = if ("osm_id" %in% names(.)) osm_id else NA_character_,
    poi_group = poi_key,   # directly map key -> group
    geometry = geometry
  )

# --- 7) Assign to neighbourhoods, count per group, compute densities ---
poi_counts_grp <- st_join(
  pois_grouped,
  neigh_m %>% dplyr::select(all_of(neigh_name)),
  join = st_within
) %>%
  st_drop_geometry() %>%
  filter(!is.na(.data[[neigh_name]])) %>%
  group_by(.data[[neigh_name]], poi_group) %>%
  summarise(n_pois = n(), .groups = "drop")

neigh_with_density <- neigh_m %>%
  mutate(area_km2 = as.numeric(st_area(.)) / 1e6) %>%
  left_join(
    poi_counts_grp %>%
      tidyr::pivot_wider(names_from = poi_group, values_from = n_pois, values_fill = 0),
    by = setNames(neigh_name, neigh_name)
  ) %>%
  mutate(
    dens_shops    = ifelse(area_km2 > 0, coalesce(shop,    0) / area_km2, NA_real_),
    dens_amenity  = ifelse(area_km2 > 0, coalesce(amenity, 0) / area_km2, NA_real_),
    dens_leisure  = ifelse(area_km2 > 0, coalesce(leisure, 0) / area_km2, NA_real_),
    dens_office   = ifelse(area_km2 > 0, coalesce(office,  0) / area_km2, NA_real_)
  )

# --- 8) Tidy output & save ---
neigh_poi_simple <- neigh_with_density %>%
  dplyr::select(all_of(neigh_name), area_km2,
                shop, amenity, leisure, office,
                dens_shops, dens_amenity, dens_leisure, dens_office,
                geom)

st_write(neigh_poi_simple, "neigh_with_umbrella_poi_density.gpkg",
         layer = "poi_density_simple", delete_dsn = TRUE, quiet = TRUE)
write.csv(st_drop_geometry(neigh_poi_simple),
          "neigh_with_umbrella_poi_density.csv", row.names = FALSE)

cat("Done:\n- neigh_with_umbrella_poi_density.gpkg (layer: poi_density_simple)\n- neigh_with_umbrella_poi_density.csv\n")


#########merge 
individuals_with_dens <- people_with_stats %>%
  left_join(
    neigh_poi_simple %>% st_drop_geometry(),
    by = "NAMN"
  )

# Save
if (inherits(individuals_with_dens, "sf")) {
  st_write(individuals_with_dens, "individuals_with_poi_densities.gpkg",
           layer = "people_with_dens", delete_dsn = TRUE, quiet = TRUE)
}
write.csv(st_drop_geometry(individuals_with_dens),
          "individuals_with_poi_densities.csv", row.names = FALSE)