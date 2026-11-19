# Diagnostic labels and saved-result readers

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble); library(stringr)
})

DIAG_LEVELS <- c(
  "Stable Diet", "Emerging Shift",
  "Ghost Shift",
  "Niche Compression/Expansion", "Internal Rebalancing",
  "Niche Restructuring",
  "Partial Diet Shift", "Structural Shift", "Major Shift"
)

INCONCLUSIVE_LAB <- "Inconclusive"

FAMILY_LEVELS <- c("Stability", "Substitution",
                   "Functional change", "Reorganisation")

STATE_TO_FAMILY <- c(
  "Stable Diet"                 = "Stability",
  "Emerging Shift"              = "Stability",
  "Ghost Shift"                 = "Substitution",
  "Niche Compression/Expansion" = "Functional change",
  "Internal Rebalancing"        = "Functional change",
  "Niche Restructuring"         = "Functional change",
  "Partial Diet Shift"          = "Reorganisation",
  "Structural Shift"            = "Reorganisation",
  "Major Shift"                 = "Reorganisation"
)

STATE_LOGIC <- tibble::tribble(
  ~diagnostic,                    ~composition,        ~H_prime,          ~Bs,
  "Stable Diet",                  "ns",                "ns",              "ns",
  "Emerging Shift",               "trend (0.05-0.10)", "ns",              "ns",
  "Ghost Shift",                  "significant",       "ns",              "ns",
  "Niche Compression/Expansion",  "ns",                "ns",              "significant",
  "Internal Rebalancing",         "ns",                "significant",     "ns",
  "Niche Restructuring",          "ns",                "significant",     "significant",
  "Partial Diet Shift",           "significant",       "significant",     "ns",
  "Structural Shift",             "significant",       "ns",              "significant",
  "Major Shift",                  "significant",       "significant",     "significant"
)

family_of <- function(diag) {
  factor(unname(STATE_TO_FAMILY[as.character(diag)]), levels = FAMILY_LEVELS)
}

add_families <- function(df, keys = character(0),
                         diag_col = "diagnostic", val_col = "pct") {
  long <- df %>%
    mutate(.family = family_of(.data[[diag_col]])) %>%
    filter(!is.na(.family)) %>%
    group_by(across(all_of(c(keys, ".family")))) %>%
    summarise(.v = sum(.data[[val_col]], na.rm = TRUE), .groups = "drop")

  wide <- pivot_wider(long, names_from = .family, values_from = .v, values_fill = 0)
  for (m in setdiff(FAMILY_LEVELS, names(wide))) wide[[m]] <- 0
  select(wide, all_of(c(keys, FAMILY_LEVELS)))
}

CONTRAST_LEVELS <- c("P1a", "P1b", "P2", "PTb", "PTa")

CONTRAST_FROM_PERIODS <- c(
  "2004|2006"                = "P1a",
  "2004-2005|2006"           = "P1b",
  "2018|2019"                = "P2",
  "2006|2018"                = "PTb",
  "2004-2006|2018-2019"      = "PTa"
)

CONTRAST_RECODE <- c(P1 = "P1a", P3 = "P1b", P2 = "P2", P4 = "PTb", PT = "PTa")

CONTRAST_LAB_1L <- c(
  P1a = "2004 vs 2006",
  P1b = "2004-2005 vs 2006",
  P2  = "2018 vs 2019",
  PTb = "2006 vs 2018",
  PTa = "2004-2006 vs 2018-2019"
)

CONTRAST_LAB <- c(
  P1a = "2004\nvs 2006",
  P1b = "2004-2005\nvs 2006",
  P2  = "2018\nvs 2019",
  PTb = "2006\nvs 2018",
  PTa = "2004-2006\nvs 2018-2019"
)

contrast_label <- function(x) {
  factor(unname(CONTRAST_LAB_1L[as.character(x)]), levels = unname(CONTRAST_LAB_1L))
}

INTRA_CONTRASTS <- c("P1a", "P1b", "P2")
INTER_CONTRASTS <- c("PTb", "PTa")

contrast_from_periods <- function(p1, p2) {
  key <- paste(p1, p2, sep = "|")
  factor(unname(CONTRAST_FROM_PERIODS[key]), levels = CONTRAST_LEVELS)
}

relabel_contrast <- function(x) {
  factor(unname(CONTRAST_RECODE[as.character(x)]), levels = CONTRAST_LEVELS)
}

SPATIAL_LEVELS_ORD <- c("all_gulf", "ecoregion", "stratum")

LEVEL_LAB <- c(
  all_gulf  = "Gulf-wide (no spatial split)",
  ecoregion = "By ecoregion",
  stratum   = "By survey stratum"
)

LEVEL_SHORT <- c(all_gulf = "Gulf", ecoregion = "Ecoregion", stratum = "Stratum")

LEVEL_SPATIAL_COL <- c(all_gulf = "Area", ecoregion = "Area", stratum = "str")

if (!exists("PREY_FAMILY")) PREY_FAMILY <- "_1"

if (!exists("RDA_DIRS")) {
  RDA_DIRS <- c(all_gulf  = paste0("Sensitivity_all_gulf",  PREY_FAMILY),
                ecoregion = paste0("Sensitivity_ecoregion", PREY_FAMILY),
                stratum   = paste0("Sensitivity_stratum",   PREY_FAMILY))
}

FAMILY_PAL <- c(
  "Stability"         = "#2F4A5A",
  "Substitution"      = "#D96C4A",
  "Functional change" = "#7F9CAB",
  "Reorganisation"    = "#9B2226"
)

FAMILY_SHAPE <- c("Stability" = 16, "Substitution" = 15,
                  "Functional change" = 18, "Reorganisation" = 17)

CURRENCY_LAB <- c(biomass = "Biomass", occurrence = "Occurrence")
CURRENCY_PAL <- c(biomass = "#2c7bb6", occurrence = "#d7191c")

LEVEL_PAL <- c(all_gulf = "#4C6A78", ecoregion = "#C08552", stratum = "#5F7A4F")

decade_bands <- function() {
  list(
    ggplot2::annotate("rect", xmin = 0.5, xmax = 3.5, ymin = -Inf, ymax = Inf,
                      fill = "#D7E6EC", alpha = 0.55),
    ggplot2::annotate("rect", xmin = 3.5, xmax = 5.5, ymin = -Inf, ymax = Inf,
                      fill = "#F6E0D6", alpha = 0.65),
    ggplot2::geom_vline(xintercept = 3.5, linetype = "dashed",
                        colour = "grey60", linewidth = 0.3)
  )
}

PERIOD_PAL <- c("2004-2006" = "#2c7bb6", "2018-2019" = "#d7191c")

AREA_LAB <- c("Baie des Chaleurs" = "Chaleur Bay",
              "Chaleur Bay"       = "Chaleur Bay")
area_label <- function(x) {
  x <- as.character(x)
  ifelse(x %in% names(AREA_LAB), unname(AREA_LAB[x]), x)
}

if (!exists("FIG_FONT")) FIG_FONT <- "sans"

theme_diag <- function(base_size = 11, base_family = FIG_FONT) {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "top",
      legend.title     = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold"),
      axis.text.x      = ggplot2::element_text(size = 8),
      plot.title       = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.caption     = ggplot2::element_text(hjust = 0, size = base_size - 2)
    )
}
if (requireNamespace("ggplot2", quietly = TRUE)) ggplot2::theme_set(theme_diag())

parse_run_name <- function(path) {
  b  <- sub("\\.rda$", "", basename(path))
  rx <- paste0("^(biomass|occurrence)_([A-Za-z0-9_]+)_(",
               paste(SPATIAL_LEVELS_ORD, collapse = "|"), ")_(.+)_vs_(.+)$")
  m  <- regmatches(b, regexec(rx, b))[[1]]
  if (!length(m)) {
    return(list(currency = NA_character_, scenario = NA_character_,
                family = NA_character_, level = NA_character_,
                period_1 = NA_character_, period_2 = NA_character_))
  }
  sc  <- m[3]
  fam <- if (grepl("_", sc)) sub("^[^_]+", "", sc) else ""
  sc  <- sub("_.*$", "", sc)
  list(currency = m[2], scenario = sc, family = fam, level = m[4],
       period_1 = gsub("_", "-", m[5]), period_2 = gsub("_", "-", m[6]))
}

read_run <- function(path) {
  e <- new.env(); load(path, envir = e)
  if (!exists("res", envir = e)) stop("No object 'res' in ", path)
  r <- tibble::as_tibble(e$res$results)
  if (!nrow(r)) return(NULL)

  meta <- parse_run_name(path)

  lvl <- if ("spatial_level" %in% names(r)) as.character(r$spatial_level[1]) else meta$level
  cur <- if ("mode" %in% names(r))          as.character(r$mode[1])          else meta$currency
  p1  <- if ("period_1" %in% names(r))      as.character(r$period_1[1])      else meta$period_1
  p2  <- if ("period_2" %in% names(r))      as.character(r$period_2[1])      else meta$period_2

  ct <- contrast_from_periods(p1, p2)
  if (is.na(ct)) ct <- relabel_contrast(meta$scenario)

  sp_col <- intersect(c("Area", "str"), names(r))[1]
  r$spatial_unit <- if (is.na(sp_col)) NA_character_ else as.character(r[[sp_col]])

  r %>%
    select(-any_of(c("Area", "str"))) %>%
    mutate(currency = cur, contrast = ct, level = lvl,
           period_1 = p1, period_2 = p2, .before = 1)
}

read_all_runs <- function(levels = SPATIAL_LEVELS_ORD, dirs = RDA_DIRS) {
  out <- list()
  for (lv in levels) {
    d <- dirs[[lv]]
    if (is.null(d) || !dir.exists(d)) {
      warning("Folder not found for level '", lv, "': ", d, call. = FALSE)
      next
    }
    fs <- list.files(d, pattern = "\\.rda$", full.names = TRUE)
    if (!length(fs)) warning("No .rda in ", d, call. = FALSE)
    for (f in fs) out[[length(out) + 1]] <- read_run(f)
  }
  res <- bind_rows(out)
  if (!nrow(res)) stop("No results read. Check RDA_DIRS and run 04_Run_Main_Analyses.R first.")
  res %>% mutate(level = factor(level, levels = SPATIAL_LEVELS_ORD))
}

freq_by_resolution <- function(results, keys = c("currency", "contrast", "level")) {
  base <- results %>% filter(diagnostic %in% DIAG_LEVELS)
  if (!nrow(base)) return(base[0, ])

  cnt <- base %>% count(across(all_of(c(keys, "x_threshold"))), diagnostic, name = "n")
  tot <- cnt %>%
    group_by(across(all_of(c(keys, "x_threshold")))) %>%
    summarise(n_tot = sum(n), .groups = "drop")

  grid <- tot %>%
    select(all_of(c(keys, "x_threshold"))) %>%
    distinct() %>%
    tidyr::crossing(diagnostic = DIAG_LEVELS)

  grid %>%
    left_join(cnt, by = c(keys, "x_threshold", "diagnostic")) %>%
    left_join(tot, by = c(keys, "x_threshold")) %>%
    mutate(n = dplyr::coalesce(n, 0L), pct = 100 * n / n_tot) %>%
    group_by(across(all_of(c(keys, "diagnostic")))) %>%
    summarise(pct = mean(pct), .groups = "drop") %>%
    mutate(diagnostic = factor(diagnostic, levels = DIAG_LEVELS))
}

family_by_resolution <- function(results, keys = c("currency", "contrast")) {
  results %>%
    filter(diagnostic %in% DIAG_LEVELS) %>%
    mutate(family = family_of(diagnostic)) %>%
    count(across(all_of(c(keys, "x_threshold"))), family, name = "n") %>%
    group_by(across(all_of(c(keys, "x_threshold")))) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup() %>%
    select(-n) %>%
    tidyr::complete(tidyr::nesting(!!!rlang::syms(c(keys, "x_threshold"))),
                    family = factor(FAMILY_LEVELS, levels = FAMILY_LEVELS),
                    fill = list(pct = 0))
}

spread_by_resolution <- function(results, keys = c("currency", "contrast")) {
  family_by_resolution(results, keys) %>%
    group_by(across(all_of(c(keys, "family")))) %>%
    summarise(pct_mean = mean(pct), pct_sd = sd(pct),
              pct_min = min(pct), pct_max = max(pct),
              n_res = dplyr::n_distinct(x_threshold), .groups = "drop")
}

spread_by_resolution_units <- function(results, keys = c("currency", "contrast", "level"),
                                       unit_col = "spatial_unit") {
  family_by_resolution(results, keys = c(keys, unit_col)) %>%
    group_by(across(all_of(c(keys, "x_threshold", "family")))) %>%
    summarise(pct = mean(pct), .groups = "drop") %>%
    group_by(across(all_of(c(keys, "family")))) %>%
    summarise(pct_mean = mean(pct), pct_sd = sd(pct),
              pct_min = min(pct), pct_max = max(pct),
              n_res = dplyr::n_distinct(x_threshold), .groups = "drop")
}

n_units_by <- function(results, keys = c("currency", "contrast", "level")) {
  results %>%
    distinct(across(all_of(c(keys, "species", "size_class")))) %>%
    count(across(all_of(keys)), name = "n_units")
}

coverage_by <- function(results, keys = c("currency", "contrast", "level")) {
  results %>%
    group_by(across(all_of(keys))) %>%
    summarise(
      n_cells        = dplyr::n(),
      pct_testable   = round(100 * mean(testable), 1),
      pct_reliable   = round(100 * mean(reliable), 1),
      pct_confounded = round(100 * mean(comp_dispersion, na.rm = TRUE), 1),
      .groups = "drop"
    )
}

if (!exists("DIR_TABLES"))  DIR_TABLES  <- "Output_Tables"
if (!exists("DIR_FIGURES")) DIR_FIGURES <- "Output_Figures"
if (!exists("DIR_APPEND"))  DIR_APPEND  <- "Output_Appendices"

for (d in c(DIR_TABLES, DIR_FIGURES, DIR_APPEND)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

load_strata <- function(rds_path = "strata_rv_gulf.rds") {
  shp <- suppressWarnings(
    system.file("extdata/shapefiles/survey.stratum.polygons.shp",
                package = "gulf.spatial"))
  if (nzchar(shp) && file.exists(shp)) {
    s <- sf::read_sf(shp)
    s <- s[s$survey == "rv" & s$region == "gulf" & s$type == "polygon", ]
    s$str <- as.character(s$stratum)
    return(s)
  }
  if (file.exists(rds_path)) return(readRDS(rds_path))
  NULL
}

save_fig <- function(plot, stem, width, height, dpi = 300, dir = DIR_FIGURES) {
  ggplot2::ggsave(file.path(dir, paste0(stem, ".png")), plot,
                  width = width, height = height, dpi = dpi,
                  bg = "white", limitsize = FALSE)
  cat("  wrote ", stem, ".png\n", sep = "")
  invisible(plot)
}

message("Result mappings loaded: ",
        length(DIAG_LEVELS), " states -> ", length(FAMILY_LEVELS), " families, ",
        length(SPATIAL_LEVELS_ORD), " spatial levels.")

suppressPackageStartupMessages({library(dplyr); library(purrr); library(tibble)})

`%||%` <- function(a, b) if (is.null(a)) b else a

.pooled_path <- function(level, fam, tag = "", base = "Sensitivity") {
  file.path("data", base, sprintf("all_runs_%s%s%s.rda", level, fam, tag))
}

.extract <- function(o) {
  if (is.data.frame(o)) return(o)
  if (is.list(o) && !is.null(o$results)) return(o$results)
  NULL
}

load_all_runs <- function(level = "all_gulf", prey_family = "_1", tag = "",
                          base = "Sensitivity",
                          dir_override = NULL, verbose = TRUE) {

  pooled <- .pooled_path(level, prey_family, tag, base)
  res <- NULL

  if (file.exists(pooled)) {
    e <- new.env(); nm <- load(pooled, envir = e)
    res <- bind_rows(map(nm, ~ .extract(get(.x, envir = e))))
    if (verbose) message("Charge depuis ", pooled, " (", length(nm), " objets)")
  } else {
    sens_dir <- dir_override %||% sprintf("%s_%s%s%s", base, level, prey_family, tag)
    if (!dir.exists(sens_dir)) {
      dirs <- grep("^Sensitivity", list.dirs(".", recursive = FALSE,
                                             full.names = FALSE), value = TRUE)
      stop("Ni ", pooled, " ni le dossier ", sens_dir, ".",
           "\nDossiers Sensitivity presents : ",
           if (length(dirs)) paste(dirs, collapse = ", ") else "aucun")
    }
    f <- list.files(sens_dir, pattern = "\\.rda$", full.names = TRUE)
    if (!length(f)) stop("Aucun .rda dans ", sens_dir)
    if (verbose) message("Assemblage depuis ", sens_dir, " (", length(f), " fichiers)")
    res <- bind_rows(map(f, function(ff) {
      e <- new.env(); nm <- load(ff, envir = e)
      bind_rows(map(nm, ~ .extract(get(.x, envir = e))))
    }))
  }

  if (is.null(res) || !nrow(res)) stop("Aucune ligne de resultats recuperee.")

  res$contrast <- contrast_from_periods(res$period_1, res$period_2)
  if (anyNA(res$contrast)) {
    bad <- res %>% filter(is.na(contrast)) %>% distinct(period_1, period_2)
    warning("Paires de periodes non reconnues : ",
            paste(sprintf("%s|%s", bad$period_1, bad$period_2), collapse = ", "))
  }
  res$contrast_type <- ifelse(as.character(res$contrast) %in% INTRA_CONTRASTS,
                              "within", "between")
  res$contrast_label <- contrast_label(res$contrast)

  res$family <- as.character(family_of(res$diagnostic))

  unit <- intersect(c("Area", "str"), names(res))[1]
  res$cell <- paste(res$species, res$size_class, res[[unit]], sep = "|")
  attr(res, "unit_col") <- unit

  if (verbose) {
    message(sprintf("%d lignes | %d cellules | %d seuils | contrastes: %s",
                    nrow(res), dplyr::n_distinct(res$cell),
                    dplyr::n_distinct(res$x_threshold),
                    paste(levels(droplevels(res$contrast)), collapse = ", ")))
  }
  res
}

inventory_runs <- function() {
  pooled <- unlist(lapply(c("Sensitivity", "SensitivityT"), function(b) {
    f <- list.files(file.path("data", b), pattern = "^all_runs_.*\\.rda$")
    if (length(f)) file.path("data", b, f) else character(0)
  }))
  dirs <- grep("^Sensitivity_", list.dirs(".", recursive = FALSE,
                                          full.names = FALSE), value = TRUE)
  dirs <- dirs[!grepl("^Sensitivity_Plot", dirs)]
  list(
    pooled   = tibble(fichier = pooled),
    dossiers = map_dfr(dirs, ~ tibble(dossier = .x,
                                      n_rda = length(list.files(.x, pattern = "\\.rda$"))))
  )
}
