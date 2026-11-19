# Prey grouping, survey covariates and size classes
# Input: data/prey_taxon_data.rda; data/diet_corr.RData
# Output: data/prey_groups.RData; data/prey_groups_PP.RData; data/Database.rda; data/diet_clean.rda; data/dat_classed.rda

if (!file.exists("data/prey_taxon_data.rda")) stop("Run 01_Extract_and_Clean.R first.")

# Pooled prey-grouping thresholds
local({
project_path <- here::here()

library(data.table)
source("R_helpers/Taxonomy.R", local = TRUE)
load(paste0(project_path, "/data/prey_taxon_data.rda"))

prey_groups <- rbind(species_prey_id, species_prey_noid, fill = TRUE)
data.table::setDT(prey_groups)

is_prey <- !is.na(prey_groups$prey_category) & prey_groups$prey_category == "identified_prey_species"
pc_chr  <- as.character(prey_groups$prey_category)

data.table::set(prey_groups, j = "prey_category_lowest",
                value = ifelse(is_prey & !is.na(prey_groups$verified_name), as.character(prey_groups$verified_name), pc_chr))
data.table::set(prey_groups, j = "tax_level_lowest",
                value = ifelse(is_prey & !is.na(prey_groups$verified_name), as.character(prey_groups$tax_level), "prey_category"))

chosen_ranks <- c("phylum", "subphylum", "class", "subclass", "order", "infraorder",
                  "family", "genus", "species")
valid_cols   <- intersect(ranknfile[ranknfile %in% chosen_ranks], names(prey_groups))
phyla_list   <- unique(prey_groups[kingdom == "Animalia" & !is.na(phylum), phylum])
target_rows  <- prey_groups$phylum %in% phyla_list
thresholds   <- seq(10, 1000, by = 10)

counts <- rank_counts(prey_groups, rev(valid_cols), target_rows)

phylum_tot <- prey_groups[is_prey,
                          .(total_stomach = uniqueN(unlist(list_stomach_id)),
                            total_pred    = uniqueN(unlist(list_predator_id))),
                          by = phylum]
genus_vec  <- as.character(prey_groups$genus)
phylum_vec <- as.character(prey_groups$phylum)

build_resolution <- function(x, n_pred, collapse_genus, col_name, col_level) {
  r <- PreyCategory(prey_groups, chosen_ranks, ranknfile, target_rows = target_rows,
                    n_stomach_threshold = x, n_predator_threshold = n_pred, counts = counts)
  name <- r$name; level <- r$level

  name[!is_prey]  <- pc_chr[!is_prey]
  level[!is_prey] <- "prey_category"

  if (collapse_genus) {
    g   <- unique(name[!is.na(level) & level == "genus"])
    idx <- which(level == "species" & genus_vec %in% g)
    name[idx]  <- genus_vec[idx]
    level[idx] <- "genus"
  }

  bad <- if (n_pred >= 2) phylum_tot[total_stomach < x | total_pred < n_pred, phylum] else
    phylum_tot[total_stomach < x, phylum]
  idx <- which(level == "phylum" & phylum_vec %in% bad)
  name[idx] <- "other_phyla"

  data.table::set(prey_groups, j = col_name,  value = name)
  data.table::set(prey_groups, j = col_level, value = level)
  invisible(NULL)
}

for (x in thresholds) {
  build_resolution(x, n_pred = 2, collapse_genus = TRUE,
                   col_name = paste0("prey_category_", x, "_2"), col_level = paste0("tax_level_", x, "_2"))
  message("[_2] x = ", x, " done")
}

for (x in thresholds) {
  build_resolution(x, n_pred = 1, collapse_genus = FALSE,
                   col_name = paste0("prey_category_", x, "_1"), col_level = paste0("tax_level_", x, "_1"))
  message("[_1] x = ", x, " done")
}

n_cat <- rbindlist(lapply(thresholds, function(x) data.table(
  x = x,
  n_2 = uniqueN(prey_groups[[paste0("prey_category_", x, "_2")]]),
  n_1 = uniqueN(prey_groups[[paste0("prey_category_", x, "_1")]]))))
cat("\nNumber of prey categories by threshold:\n")
print(n_cat[x %in% c(10, 50, 100, 200, 440, 750, 1000)])

build_resolution(440, n_pred = 2, collapse_genus = TRUE,
                 col_name = "prey_category_keep", col_level = "tax_level_keep")

keep_name  <- prey_groups$prey_category_keep
keep_level <- prey_groups$tax_level_keep
is_taxo    <- !is.na(keep_level) & keep_level != "prey_category"
has_child  <- rep(FALSE, nrow(prey_groups))

rank_mat <- lapply(c("phylum", "class", "order", "family", "genus"), function(cc) as.character(prey_groups[[cc]]))
for (cat_name in unique(keep_name[is_taxo])) {
  under <- Reduce(`|`, lapply(rank_mat, function(v) !is.na(v) & v == cat_name)) & is_taxo
  if (length(unique(keep_name[under])) > 1) has_child[!is.na(keep_name) & keep_name == cat_name] <- TRUE
}
data.table::set(prey_groups, j = "has_child", value = has_child)
data.table::set(prey_groups, j = "prey_category_new",
                value = ifelse(has_child, paste0(keep_name, "_others"), keep_name))

cat("\nCategories at x = 440: ", uniqueN(prey_groups$prey_category_keep),
    " (", uniqueN(prey_groups$prey_category_new), " with the sibling rule)\n", sep = "")

save(prey_groups, file = paste0(project_path, "/data/prey_groups.RData"))

})

# Per-predator prey-grouping thresholds
local({
project_path <- here::here()
library(data.table)
source("R_helpers/Taxonomy.R", local = TRUE)
load(paste0(project_path, "/data/prey_groups.RData"))
load(paste0(project_path, "/data/diet_corr.RData"))
data.table::setDT(prey_groups); data.table::setDT(diet_corr)

PP_SUFFIX          <- "_PP"
thresholds_PP      <- seq(5, 250, by = 5)
MAIN_PRED_MIN_STOM <- 100L

chosen_ranks <- c("phylum", "subphylum", "class", "subclass", "order", "infraorder",
                  "family", "genus", "species")
phyla_list   <- unique(prey_groups[kingdom == "Animalia" & !is.na(phylum), phylum])

rank_cols <- intersect(
  c("kingdom","subkingdom","infrakingdom","phylum","subphylum","infraphylum","parvphylum",
    "gigaclass","superclass","class","subclass","infraclass","subterclass","superorder",
    "order","suborder","infraorder","parvorder","section","subsection","superfamily",
    "family","subfamily","tribe","genus","subgenus","species"),
  names(prey_groups))

lookup_cols <- c("prey_species_common_name", "prey_category", "verified_name",
                 "tax_level", rank_cols)
tax_lookup  <- unique(prey_groups[, ..lookup_cols], by = "prey_species_common_name")

main_predators <- diet_corr[!is.na(stomach_id),
                            .(n_sto = uniqueN(stomach_id)),
                            by = predator_species_common_name
][n_sto >= MAIN_PRED_MIN_STOM, predator_species_common_name]
message(length(main_predators), " main predators")

rec_all <- diet_corr[!is.na(stomach_id) & predator_species_common_name %in% main_predators,
                     .(predator = predator_species_common_name, stomach_id, prey_species_common_name)]
rec_all <- tax_lookup[rec_all, on = "prey_species_common_name", nomatch = 0L]

by_cols <- c("prey_species_common_name", "prey_category", "verified_name",
             "tax_level", rank_cols)

prey_groups_PP <- rbindlist(lapply(main_predators, function(p) {

  pg <- rec_all[predator == p,
                .(list_stomach_id  = list(unique(stomach_id)),
                  list_predator_id = list(p)),
                by = by_cols]
  if (!nrow(pg)) return(NULL)

  is_prey     <- !is.na(pg$prey_category) & pg$prey_category == "identified_prey_species"
  pc_chr      <- as.character(pg$prey_category)
  phylum_vec  <- as.character(pg$phylum)
  target_rows <- pg$phylum %in% phyla_list
  valid_cols  <- intersect(ranknfile[ranknfile %in% chosen_ranks], names(pg))
  counts      <- rank_counts(pg, rev(valid_cols), target_rows)
  phylum_tot  <- pg[is_prey, .(total_stomach = uniqueN(unlist(list_stomach_id))), by = phylum]

  for (x in thresholds_PP) {
    r <- PreyCategory(pg, chosen_ranks, ranknfile, target_rows = target_rows,
                      n_stomach_threshold = x, n_predator_threshold = 1, counts = counts)
    name <- r$name; level <- r$level
    name[!is_prey]  <- pc_chr[!is_prey]
    level[!is_prey] <- "prey_category"
    bad <- phylum_tot[total_stomach < x, phylum]
    idx <- which(level == "phylum" & phylum_vec %in% bad)
    name[idx] <- "other_phyla"
    data.table::set(pg, j = paste0("prey_category_", x, PP_SUFFIX), value = name)
    data.table::set(pg, j = paste0("tax_level_", x, PP_SUFFIX),     value = level)
  }

  pg[, predator := p]
  message("[PP] ", p, " (", nrow(pg), " taxa)")
  pg
}), fill = TRUE)

data.table::setcolorder(prey_groups_PP,
                        c("predator", "prey_species_common_name", "prey_category", "tax_level"))

n_pp_cols <- sum(grepl(paste0("^prey_category_\\d+", PP_SUFFIX, "$"), names(prey_groups_PP)))
stopifnot(n_pp_cols == length(thresholds_PP))

cat(sprintf("\nprey_groups_PP: %d rows | %d predators | %d resolution columns\n",
            nrow(prey_groups_PP), uniqueN(prey_groups_PP$predator), n_pp_cols))

save(prey_groups_PP, file = paste0(project_path, "/data/prey_groups_PP.RData"))

})

# Survey covariates, spatial assignment and predator weights
local({
library(data.table)
library(gulf)
library(dplyr)
library(tidyr)
library(lubridate)
library(sf)

load("data/prey_groups.RData")
load("data/diet_corr.RData")

prey_groups_unique <- unique(prey_groups, by = "prey_species_common_name")
Database <- merge(diet_corr[, prey_species_latin_name := NULL], prey_groups_unique,
                  by = "prey_species_common_name", all.x = TRUE)

if (file.exists("data/prey_groups_PP.RData")) {
  load("data/prey_groups_PP.RData")
  pp_cols <- grep("^(prey_category|tax_level)_\\d+_PP$", names(prey_groups_PP), value = TRUE)
  pp_map  <- unique(prey_groups_PP[, c("predator", "prey_species_common_name", pp_cols), with = FALSE])
  if (anyDuplicated(pp_map, by = c("predator", "prey_species_common_name")))
    stop("prey_groups_PP: (predator, prey) key is not unique.")
  n_before <- nrow(Database)
  Database <- merge(Database, pp_map,
                    by.x = c("predator_species_common_name", "prey_species_common_name"),
                    by.y = c("predator", "prey_species_common_name"), all.x = TRUE)
  stopifnot(nrow(Database) == n_before)
  first_pp <- grep("^prey_category_\\d+_PP$", names(Database), value = TRUE)[1]
  n_na <- sum(is.na(Database[[first_pp]]) & !is.na(Database$prey_species_common_name))
  cat("PP family merged:", length(pp_cols) %/% 2, "resolution columns;",
      n_na, "prey rows without PP classification\n")
}

Database[, `:=`(
  somatic_length_cm = as.numeric(somatic_length_cm),
  number_of_prey    = fcoalesce(as.numeric(number_of_prey), 0),
  somatic_wt_g      = fcoalesce(as.numeric(somatic_wt_g), 0),
  year              = cruise_year,
  period            = factor(fcase(cruise_year %in% 2004:2006, "2004-2006",
                                   cruise_year %in% 2018:2019, "2018-2019"),
                             levels = c("2004-2006", "2018-2019"))
)]

lmat_predators <- fread("data/lmat_predators.csv", encoding = "Latin-1")
lmat_predators[lmat == "Unknown" | is.na(lmat), lmat := as.character(lmax)]
setnames(lmat_predators, 1, "predator_species_code")
refs <- lmat_predators[, .(predator_species_code, lmax, lmat)]

Database <- merge(Database, refs, by = "predator_species_code", all.x = TRUE)
Database[, somatic_length_cm_orig := somatic_length_cm]
Database[somatic_length_cm > 2 * lmax, somatic_length_cm := somatic_length_cm / 10]
Database[somatic_length_cm > lmax,     somatic_length_cm := lmax]

correction_log <- unique(
  Database[somatic_length_cm != somatic_length_cm_orig,
           .(stomach_id, predator_species_code, predator_species_common_name,
             original = somatic_length_cm_orig, corrected = somatic_length_cm, lmax,
             action = fcase(
               somatic_length_cm_orig > 2 * lmax & somatic_length_cm == somatic_length_cm_orig / 10, "unit mm->cm",
               somatic_length_cm_orig > 2 * lmax & somatic_length_cm == lmax, "unit mm->cm then cap",
               somatic_length_cm == lmax, "cap at Lmax",
               default = "other"))],
  by = "stomach_id")
print(correction_log[, .N, by = action])
fwrite(correction_log, "data/correction_log_lengths.csv")
Database[, somatic_length_cm_orig := NULL]

sets <- as.data.table(read.card(survey = "rv", sampling = "research", card.type = "set",
                                year = c(2004:2006, 2018:2019), as.data.frame = TRUE))
sets[, `:=`(
  samp.date = ymd(paste(year, month, day, sep = "-")),
  latitude  = (latitude.start + latitude.end) / 2,
  longitude = (longitude.start + longitude.end) / 2,
  depth     = (depth.start + depth.end) / 2
)]

Database <- merge(
  Database,
  sets[, .(stratum, year, month, day, samp.date, set = set.number, vessel.code,
           hour = start.hour, latitude, longitude, depth,
           bottom_temp = bottom.temperature, salinity = bottom.salinity)],
  by = c("year", "set", "vessel.code"), all.x = TRUE)

fix <- Database$samp.date == "2006-09-05" & Database$hour == 19
Database$latitude[fix]  <- mean(Database$latitude[fix],  na.rm = TRUE)
Database$longitude[fix] <- mean(Database$longitude[fix], na.rm = TRUE)

Database <- Database %>%
  mutate(Lat = latitude, Lon = longitude) %>%
  drop_na(Lon, Lat) %>%
  st_as_sf(coords = c("Lon", "Lat"), crs = 4326)

read_layer <- function(path) st_make_valid(st_transform(st_read(path, quiet = TRUE), st_crs(Database)))
new_ecoregions_clean <- read_layer("data/Spatial_data/new_ecoregions_final.shp")
new_regions_v2_clean <- read_layer("data/Spatial_data/new_regions_v2_final.shp")
new_gulf_clean       <- read_layer("data/Spatial_data/new_gulf_final.shp")
NAFO_4T_sf           <- read_layer("data/Spatial_data/NAFO/nafo_2014_02.shp")
EAR_sf               <- read_layer("data/Spatial_data/EAR_map/EAR_map.shp")

sf_use_s2(FALSE)
Database$EcoZone <- new_ecoregions_clean$EcoZone[st_nearest_feature(Database, new_ecoregions_clean)]
Database$Zone    <- new_regions_v2_clean$Zone[st_nearest_feature(Database, new_regions_v2_clean)]
Database$NAFO    <- NAFO_4T_sf$level_2[st_nearest_feature(Database, NAFO_4T_sf)]
Database$Region  <- EAR_sf$regn_nm[st_nearest_feature(Database, EAR_sf)]
Database$Area    <- new_gulf_clean$Area[st_nearest_feature(Database, new_gulf_clean)]
sf_use_s2(TRUE)

Database <- Database %>%
  mutate(EcoZone = replace_na(EcoZone, "4TF Only"),
         Zone    = replace_na(Zone, "4TF Only"),
         NAFO    = replace_na(NAFO, "4TF"),
         is_nutritional_prey = if_else(is_empty == FALSE, 1, 0))

save(Database, file = "data/Database.rda")

invalid_reasons <- c("content weight larger than stomach", "No set number", "Prey weights too large")

retained <- Database %>%
  group_by(predator_species_common_name, predator_species_latin_name) %>%
  mutate(Total_N     = n_distinct(stomach_id),
         N_2004_2006 = n_distinct(stomach_id[period == "2004-2006"]),
         N_2018_2019 = n_distinct(stomach_id[period == "2018-2019"])) %>%
  ungroup() %>%
  filter(Total_N >= 90, N_2004_2006 >= 10, N_2018_2019 >= 10,
         !predator_invalid_reason %in% invalid_reasons)

diet_C <- retained
stomach_fullness <- diet_C %>%
  st_drop_geometry() %>%
  group_by(stomach_id, period, year, predator_species_common_name) %>%
  summarise(has_valid_food = if_else(sum(is_nutritional_prey, na.rm = TRUE) > 0, 1, 0),
            .groups = "drop")
save(diet_C, stomach_fullness, file = "data/fullness.rda")

diet_clean <- retained %>%
  filter(!prey_category_keep %in% c("parasite", "empty", "digested"))
setDT(diet_clean)

# Predator weight estimates
combos <- diet_clean[!is.na(predator_species_code) & !is.na(somatic_length_cm) & !is.na(year),
                     .(n_fish = .N), by = .(predator_species_code, somatic_length_cm, year)]

cache_path <- "data/unique_samples.rda"
key_cols <- c("predator_species_code", "somatic_length_cm", "year")

if (file.exists(cache_path)) {
  load(cache_path)
  setDT(unique_samples)
  stopifnot(all(c(key_cols, "pred_weight_est") %in% names(unique_samples)))
  if (anyDuplicated(unique_samples, by = key_cols)) {
    stop("Duplicate predator-length-year combinations in the weight cache.")
  }
} else {
  unique_samples <- copy(combos[0])
  unique_samples[, pred_weight_est := numeric()]
}

missing <- combos[!unique_samples,
                  on = .(predator_species_code, somatic_length_cm, year)]

if (nrow(missing) > 0L) {
  missing[, pred_weight_est := vapply(seq_len(.N), function(i) {
    estimate_year <- if (predator_species_code[i] == 60 &&
                         year[i] %in% c(2018, 2019)) 2017 else year[i]
    tryCatch({
      w <- gulf::weight(somatic_length_cm[i], species = predator_species_code[i],
                        year = estimate_year) * 1000
      if (length(w)) as.numeric(w[1]) else NA_real_
    }, error = function(e) NA_real_)
  }, numeric(1))]

  unique_samples <- rbindlist(list(unique_samples, missing), fill = TRUE)
  save(unique_samples, file = cache_path)
  message("Predator weights added to cache: ", nrow(missing))
} else if (!file.exists(cache_path)) {
  save(unique_samples, file = cache_path)
}

message("Predator weight combinations in cache: ", nrow(unique_samples))
if (anyNA(unique_samples$pred_weight_est)) {
  warning(sum(is.na(unique_samples$pred_weight_est)),
          " combinations have no weight estimate.", call. = FALSE)
}

diet_clean[unique_samples, on = .(predator_species_code, somatic_length_cm, year),
           pred_weight_est := i.pred_weight_est]

min_wt <- diet_clean[somatic_wt_g > 0, min(somatic_wt_g, na.rm = TRUE)]
diet_clean[somatic_wt_g == 0, somatic_wt_g := min_wt]
diet_clean[, pfi := (somatic_wt_g / pred_weight_est) * 100]

save(diet_clean, file = "data/diet_clean.rda")

})

# Predator size classes
local({
library(data.table)
library(dplyr)
library(tidyr)
library(sf)

load("data/diet_clean.rda")
setDT(diet_clean)

diet_clean[, `:=`(lmat = as.numeric(lmat), lmax = as.numeric(lmax))]
diet_clean[, size_class := fcase(
  lmat == lmax,               "one_class",
  somatic_length_cm <  lmat,  "juvenile",
  somatic_length_cm >= lmat,  "adult"
)]

size_share <- diet_clean %>%
  st_drop_geometry() %>%
  count(predator_species_common_name, period, size_class) %>%
  group_by(predator_species_common_name, period) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(names_from = size_class, values_from = pct, values_fill = 0)
print(size_share, n = Inf)

dat_classed <- diet_clean
save(dat_classed, file = "data/dat_classed.rda")

})
