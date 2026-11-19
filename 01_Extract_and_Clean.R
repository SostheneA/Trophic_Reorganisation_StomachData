# Survey corrections and prey taxonomy
# Input: data/Survey_data/*.csv
# Output: data/diet_corr.RData; data/prey_taxon_data.rda

if (!dir.exists("data/Survey_data")) stop("Missing data/Survey_data")

# Survey records and stomach corrections
local({
project_path <- here::here()

library(data.table)
library(dplyr)

survey_files <- list(
  list(path = "data/Survey_data/4T_Needler_2004.csv", vessel = "N"),
  list(path = "data/Survey_data/4T_Needler_2005.csv", vessel = "N"),
  list(path = "data/Survey_data/4T_TEL_2004.csv", vessel = "T"),
  list(path = "data/Survey_data/4T_TEL_2005.csv", vessel = "T"),
  list(path = "data/Survey_data/4T_TEL_2006.csv", vessel = "T"),
  list(path = "data/Survey_data/4T_TEL_2018.csv", vessel = "T"),
  list(path = "data/Survey_data/4T_TEL_2019.csv", vessel = "T")
)

diet_list <- lapply(survey_files, function(f) {
  dt <- fread(f$path, encoding = "UTF-8")
  dt[, vessel.code := f$vessel]
  return(dt)
})

diet_raw <- rbindlist(diet_list, fill = TRUE)
diet_raw <- diet_raw %>% mutate(year = cruise_year)

setDT(diet_raw)

diet_corr <- copy(diet_raw)

diet_corr$prey_species_latin_name <- ifelse((is.na(diet_corr$prey_species_latin_name) | diet_corr$prey_species_latin_name == "") &
                                             diet_corr$prey_species_common_name  %in%
                                             c("Fish remains", "Eggs - Atlantic herring", "Mud"),
                                           diet_corr$prey_species_common_name,
                                           diet_corr$prey_species_latin_name)

diet_corr$prey_species_latin_name <- ifelse(diet_corr$prey_id == 60840, "Amphipoda", diet_corr$prey_species_latin_name)
diet_corr$prey_species_common_name <- ifelse(diet_corr$prey_id == 60840, "Amphipod unidentified", diet_corr$prey_species_common_name)
diet_corr$prey_species_code <- ifelse(diet_corr$prey_id == 60840, 2800, diet_corr$prey_species_code)

diet_corr$prey_species_latin_name <- ifelse(diet_corr$prey_id == 68342, "Hyaloidea", diet_corr$prey_species_latin_name)
diet_corr$prey_species_common_name <- ifelse(diet_corr$prey_id == 68342, "Hyaloidea", diet_corr$prey_species_common_name)
diet_corr$prey_species_code <- ifelse(diet_corr$prey_id == 68342, NA_integer_, diet_corr$prey_species_code)

diet_corr$predator_species_latin_name <- ifelse(diet_corr$predator_species_latin_name == "Notolepis rissoi",
                                               "Arctozenus risso",
                                               diet_corr$predator_species_latin_name)

diet_corr$predator_species_code <- ifelse(diet_corr$predator_species_code == 727,
                                         712,
                                         diet_corr$predator_species_code)

diet_corr$prey_species_code <- ifelse(diet_corr$prey_species_code == 910,
                                     NA_integer_,
                                     diet_corr$prey_species_code)

diet_corr$prey_species_code <- ifelse(diet_corr$prey_species_code == 2500,
                                     2507,
                                     diet_corr$prey_species_code)

taxa_not_in_all_sp <- c("Nemertea (p.)", "Malacostraca (c.)", "Acteonidae (f.)",
                        "Eusiridae (f.)",
                        "Eunicida (o.)", "Liljeborgiidae (f.)", "Ampithoidae (f.)")

diet_corr$prey_species_latin_name <- ifelse(diet_corr$prey_species_latin_name %in% taxa_not_in_all_sp,
                                           gsub(" .*", "", trimws(diet_corr$prey_species_latin_name)),
                                           diet_corr$prey_species_latin_name)

diet_corr$prey_species_latin_name <- ifelse(diet_corr$prey_species_latin_name %in% c(""),
                                            NA_character_,
                                            diet_corr$prey_species_latin_name)

diet_corr$prey_species_common_name <- ifelse(diet_corr$prey_species_common_name %in% c(""),
                                            NA_character_,
                                            diet_corr$prey_species_common_name)

potential_empty_elements <- c(
  "Parasites,round worms",
  "Cestoda (c.)",
  "Trematoda",
  "Scales - unid.",
  "Scales - unid. fish",
  "Shells - unid. fish",
  "Mollusca sp. empty",
  "Casing - unidentified invertebrate",
  "Mucus",
  "Sand",
  "Stones and rocks",
  "Foreign articles,
                  garbage",
  "Mud",
  "Foreign articles,garbage"
)

emptyness <- diet_corr[,
                       .(
                         prey_list = list({
                           u <- unique(trimws(prey_species_latin_name))
                           u <- u[!is.na(u) & nzchar(u)]
                           u
                         }),
                         codes = paste0(unique(na.omit(prey_species_code)), collapse = ","),
                         n_prey_names = uniqueN(prey_species_latin_name, na.rm = TRUE)
                       ),
                       by = .(
                         stomach_id
                       )
][
  ,
  keep := vapply(prey_list, function(x) length(x) == 0L || all(x %chin% potential_empty_elements), logical(1L))
][keep == T,

][
  ,
  `:=`(
    prey_names = vapply(prey_list, paste, character(1L), collapse = ";"),
    prey_list = NULL,
    keep = NULL
  )
]

empty_stomach_ids <- unique(emptyness$stomach_id)

diet_corr[, is_empty := fifelse(stomach_id %in% empty_stomach_ids, TRUE, FALSE)]
setDT(diet_corr)

problematic_stomach <- diet_corr[!is_empty & is.na(prey_species_common_name),]$stomach_id
problematic_stomach_table1 <- diet_corr[stomach_id %in% c(problematic_stomach),]
problematic_stomach_table <- diet_corr[stomach_id %in% c(problematic_stomach),
                       .(
                         preys = paste0(unique(prey_species_code), collapse = ",")
                       ),
                       by = .(
                         stomach_id
                       )
]

diet_corr[stomach_id %in% problematic_stomach &
            is.na(prey_species_code) &
            is.na(prey_species_common_name) &
            is.na(prey_species_latin_name),

          `:=`(
            prey_species_latin_name  = "no prey information",
            prey_species_common_name = "no prey information",
            prey_species_code        = NA_integer_
          )
]
setDT(diet_corr)

save(diet_corr, file = paste0(project_path, "/data/diet_corr.RData"))

rm(emptyness, problematic_stomach_table, problematic_stomach_table1, taxa_not_in_all_sp)

})

# Prey taxonomy and WoRMS classification
local({
project_path <- here::here()

library(data.table)
library(taxize)
library(worrms)
library(gulf)
library(here)
source("R_helpers/Taxonomy.R", local = TRUE)
load(paste0(project_path, "/data/diet_corr.RData"))
setDT(diet_corr)

species_prey <- diet_corr[,
                          .(
                            n_prey_count_per_prey = sum(number_of_prey, na.rm = TRUE),
                            n_stomach_per_prey = uniqueN(stomach_id, na.rm = TRUE),
                            n_predator_per_prey = uniqueN(predator_species_latin_name, na.rm = TRUE),
                            list_stomach_id       = list(unique(na.omit(stomach_id))),
                            list_predator_id      = list(unique(na.omit(predator_species_code)))
                          ),
                          by = .(
                            code = prey_species_code,
                            prey_species_common_name,
                            prey_species_latin_name
                          )
]

Not_list <- c("")

Inorganic_list <- c("Sand",
                  "Stones and rocks",
                  "Foreign articles,
                  garbage",
                  "Mud",
                  "Foreign articles,garbage"
)

Organic_list <- c("Mucus",
               "Organic debris",
               "Unid remains,digested",
               "Unid fish and invertebrates",
               "Unid. species",
               "Unidentified",
               "no prey information",
               "Odobenidae (f.)"
)

Unid_fish_list <- c("Unid fish and remains",
               "Fish remains",
               "Unid. Finfish",
               "Unid. fish (larvae,juvenile and adults)",
               "Unid fish and eggs")

Unid_invert_list <- c("Unid. marine invertebrata",
              "Ctenophora,coelenterata,porifera",
              "Protochordata sp.",
              "Zooplankton",
              "Invertebrate"
              )

Eggs_larvae_list <- list("Eggs - invertebrate" = c("Invertebrate", "egg"),
                         "Eggs - Invertebrate" = c("Invertebrate", "egg"),
                         "Eggs - unid." = c(NA, "egg"),
                         "Eggs - Unid." = c(NA, "egg"),
                         "Unid. Fish eggs" = c("Teleostei", "egg"),
                         "Eggs - decapoda" = c("Decapoda", "egg"),
                         "Eggs - Decapoda" = c("Decapoda", "egg"),
                         "Hermit crab eggs" = c("Paguridae", "egg"),
                         "Buccinidae sp. eggs" = c("Buccinidae", "egg"),
                         "Eggs - Atlantic herring"= c("Clupea harengus", "egg"),
                         "Larvae - unid." = c(NA, "larvae"),
                         "Larvae - Unid." = c(NA, "larvae"),
                         "Larvae - decapoda (o.)" = c("Decapoda", "larvae"),
                         "Larvae - Decapoda (o.)" = c("Decapoda", "larvae"),
                         "Larvae - crustacea (s.p.)" = c("Crustacea", "larvae"),
                         "Larvae - Crustacea (s.p.)" = c("Crustacea", "larvae"),
                         "Larvae - mollusca (p.)" = c("Mollusca", "larvae"),
                         "Larvae - Mollusca (p.)" = c("Mollusca", "larvae"),
                         "Larvae - polychaeta (c.)" = c("Polychaeta", "larvae"),
                         "Larvae - Polychaeta (c.)" = c("Polychaeta", "larvae"),
                         "Larvae - crab" = c("Pleocyemata", "larvae"),
                         "Larvae - Crab" = c("Pleocyemata", "larvae"),
                         "Larvae - homarus americanus" = c("Homarus americanus",
                                                           "larvae"),
                         "Larvae - Homarus americanus" = c("Homarus americanus",
                                                           "larvae"),
                         "Calappa megalops" = c("Calappa", "megalopae")
)

Parasite_list <- c("Parasites,round worms",
                   "Cestoda (c.)",
                   "Trematoda")

Rest_list <- c("Scales - unid.",
               "Scales - unid. fish",
               "Shells - unid. fish",
               "Mollusca sp. empty",
               "Casing - unidentified invertebrate",
               "Casing - Unidentified Invertebrate",
               "Scales - Unid.",
               "Scales - Unid. fish",
               "Shells - Unid. fish")

Algae_list <- c("Seaweed/Algae/Kelp unidentified",
                "Brown seaweed unidentified",
                "foraminiferans / forams / hole bearers",
                "Brown rockweed unidentified",
                "Red seaweed unidentified",
                "Eel grass"
                )

map_latin <- setNames(vapply(Eggs_larvae_list, `[`, character(1L), 1L),
                      names(Eggs_larvae_list))
map_stage <- setNames(vapply(Eggs_larvae_list, `[`, character(1L), 2L),
                      names(Eggs_larvae_list))

setDT(species_prey)
species_prey[prey_species_latin_name %in% names(Eggs_larvae_list),
             `:=`(
               prey_species_latin_name = map_latin[prey_species_latin_name],
               life_stage              = map_stage[prey_species_latin_name],
               code = NA_integer_
             )]

species_prey[prey_species_latin_name %in% Unid_fish_list,
             `:=`(
               prey_species_latin_name = "Teleostei",
               code = NA_integer_
             )]

species_prey[, prey_category := fcase(
  prey_species_latin_name %in% Not_list, "empty",
  prey_species_latin_name %in% Inorganic_list, "digested",
  prey_species_latin_name %in% Organic_list, "digested",
  prey_species_latin_name %in% Unid_invert_list, "digested",
  (!is.na(life_stage)), "egg_larvae",
  (!is.na(life_stage) & is.na(prey_species_latin_name)), "digested",
  prey_species_latin_name %in% Rest_list, "digested",
  prey_species_latin_name %in% Parasite_list, "parasite",
  default = "identified_prey_species"
)]

species_prey[is.na(prey_species_latin_name) & is.na(prey_species_common_name),
          `:=`(
            prey_category  = "empty")
]

species_prey$prey_category <- ifelse(species_prey$prey_species_common_name %in%  Algae_list,
                                     "plantae_chromista",
                                     species_prey$prey_category)

all.sp <- as.data.table(gulf::species.names)
all.sp[, aphia_id_origin := aphia_id]

all.sp[code == 290, aphia_id_origin := 126417]
all.sp[code == 2699, aphia_id_origin := 1086]
all.sp[code == 2902, aphia_id_origin := 1080]
all.sp[code == 2901, aphia_id_origin := 1080]
all.sp[code == 4998, aphia_id_origin := 325342]
all.sp[code == 5200, aphia_id_origin := 382158]
all.sp[code == 2507, aphia_id_origin := 106670]
all.sp[code == 3199, aphia_id_origin := 883]
all.sp[code == 3101, aphia_id_origin := 883]
all.sp[code == 3102, aphia_id_origin := 883]
all.sp[code == 4310, aphia_id_origin := 105]
all.sp[code == 6700, aphia_id_origin := 123111]

species_prey <- merge(all.sp[, .(code, aphia_id_origin)],
                      species_prey,
                      by = "code",
                      all.y = T
)

species_prey[, prey_size_cat := NA_character_]
species_prey[, prey_size_cat := fifelse(code %in% c(3102, 2902), "small", prey_size_cat)]
species_prey[, prey_size_cat := fifelse(code %in% c(3101, 2901), "large", prey_size_cat)]

aphia_id_origins  <- unique(species_prey[aphia_id_origin > 0, aphia_id_origin])

aphia_id_origins_batches <- split(aphia_id_origins, ceiling(seq_along(aphia_id_origins) / 50))

aphia_ids_verified_list <- lapply(aphia_id_origins_batches, function(batche) {
  tryCatch({
    Sys.sleep(0.2)
    return(wm_record(id = batche))
  }, error = function(e) {
    warning("Error on this batch of ids : ", paste(head(batche), collapse=", "))
    return(NULL)
  })
})

aphia_ids_verified_df <- rbindlist(aphia_ids_verified_list, fill = TRUE)
setDT(aphia_ids_verified_df)
aphia_ids_verified_df[, valid_AphiaID := fifelse(is.na(valid_AphiaID), parentNameUsageID, valid_AphiaID)]
aphia_id_origins_verified <- unique(aphia_ids_verified_df$valid_AphiaID)

if(!exists("class_list_prey_w")){
class_list_prey_w <- taxize::classification(aphia_id_origins_verified, db = "worms", return_id = TRUE)
}
class_tbl_prey  <- rbindlist(class_list_prey_w[!is.na(class_list_prey_w)],
                             fill = TRUE, idcol = "aphia_id")

class_tbl_prey[, rank := tolower(rank)]
class_tbl_prey[rank == "phylum (division)", rank := "phylum"]
class_tbl_prey[rank == "subphylum (subdivision)", rank := "subphylum"]
class_tbl_prey[, rank := factor(rank, levels = tolower(ranknfile), ordered = TRUE)]

class_tbl_prey[, ":="(verified_name = .SD[rank == max(rank), name],
                      tax_level = .SD[, max(rank)]),
               by = aphia_id]

class_tbl_prey <- dcast(class_tbl_prey,
                        aphia_id + verified_name + tax_level ~ rank,
                        value.var = "name")

setDT(class_tbl_prey)

class_tbl_prey <- merge(
  x = class_tbl_prey[, aphia_id := as.integer(aphia_id)],
  y = aphia_ids_verified_df[, .(aphia_id_origin = as.integer(AphiaID), aphia_id = as.integer(valid_AphiaID))],
  by = "aphia_id",
  all.y = TRUE
)

species_prey_aphia <- merge(
  x = species_prey,
  y = class_tbl_prey,
  by = "aphia_id_origin",
  all.y = TRUE
)

aphia_id_origins_rest <- setdiff(aphia_id_origins, unique(class_tbl_prey$aphia_id_origin))
if (length(aphia_id_origins_rest)) message("AphiaIDs without classification: ",
                                           paste(aphia_id_origins_rest, collapse = ", "))

latin_names  <- unique(species_prey[is.na(aphia_id_origin) &
                                                !is.na(prey_species_latin_name) &
                                                !(prey_species_latin_name %in% c("Obsolete", "Parasites,round worms")) &
                                                prey_category %in% c("identified_prey_species", "plantae_chromista", "parasite", "egg_larvae"),
                                              prey_species_latin_name])

latin_names_batches <- split(latin_names, ceiling(seq_along(latin_names) / 50))

latin_names_verified_list <- lapply(latin_names_batches, function(batche) {
  tryCatch({
    Sys.sleep(0.2)
    return(wm_records_names(name = batche))
  }, error = function(e) {
    warning("Error on this batch of names : ", paste(head(batche), collapse=", "))
    return(NULL)
  })
})

latin_names_verified_df <- dplyr::bind_rows(latin_names_verified_list)
setDT(latin_names_verified_df)
latin_names_verified_df[, valid_AphiaID := fifelse(is.na(valid_AphiaID), parentNameUsageID, valid_AphiaID)]
latin_names_id_verified <- unique(latin_names_verified_df$valid_AphiaID)

if(!exists("class_list_prey_w_latin")){
  class_list_prey_w_latin <- taxize::classification(latin_names_id_verified,
                                                    db = "worms",
                                                    return_id = TRUE)
}

class_tbl_prey2  <- rbindlist(class_list_prey_w_latin[!is.na(class_list_prey_w_latin)],
                             fill = TRUE, idcol = "aphia_id")

class_tbl_prey2[, rank := tolower(rank)]
class_tbl_prey2[rank == "phylum (division)", rank := "phylum"]
class_tbl_prey2[rank == "subphylum (subdivision)", rank := "subphylum"]
class_tbl_prey2[, rank := factor(rank, levels = tolower(ranknfile), ordered = TRUE)]

class_tbl_prey2[, ":="(verified_name = .SD[rank == max(rank), name],
                      tax_level = .SD[, max(rank)]),
               by = aphia_id]

class_tbl_prey2 <- dcast(class_tbl_prey2,
                        aphia_id + verified_name + tax_level ~ rank,
                        value.var = "name")

setDT(class_tbl_prey2)

class_tbl_prey2 <- merge(
  x = class_tbl_prey2[, aphia_id := as.integer(aphia_id)],
  y = latin_names_verified_df[, .(prey_species_latin_name = as.character(scientificname),
                                  aphia_id = as.integer(valid_AphiaID))],
  by = "aphia_id",
  all.y = TRUE
)

species_prey_latin <- merge(
  x = subset(species_prey, is.na(aphia_id_origin)),
  y = class_tbl_prey2,
  by = "prey_species_latin_name",
  all.y = T
)

latin_names_rest <- setdiff(latin_names, unique(class_tbl_prey2$prey_species_latin_name))
if (length(latin_names_rest)) message("Latin names without classification: ",
                                      paste(latin_names_rest, collapse = ", "))

species_prey_id <- rbind(species_prey_aphia, species_prey_latin, fill = TRUE)

species_prey_noid_list <- setdiff(unique(species_prey$prey_species_latin_name), unique(species_prey_id$prey_species_latin_name))
species_prey_noid <- species_prey[species_prey$prey_species_latin_name %in% species_prey_noid_list,]

cat("Classified prey names: ", uniqueN(species_prey_id$verified_name),
    " | unclassified: ", nrow(species_prey_noid), "\n", sep = "")

save(species_prey_id, species_prey_noid, file = paste0(project_path, "/data/prey_taxon_data.rda"))

})
