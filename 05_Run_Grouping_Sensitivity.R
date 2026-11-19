# Prey-grouping sensitivity

input_file <- "data/dat_classed.rda"
if (!file.exists(input_file)) stop("Missing input: ", input_file)
if (!file.exists("03_Statistical_Pipeline.R")) stop("Missing 03_Statistical_Pipeline.R")

# Available prey-grouping columns
input_env <- new.env()
load(input_file, envir = input_env)
if (!exists("dat_classed", envir = input_env, inherits = FALSE)) {
  stop("Object 'dat_classed' is missing from ", input_file)
}
prey_columns <- names(input_env$dat_classed)
rm(input_env)

families <- c("_2", "_PP")
spatial_levels <- c("all_gulf", "ecoregion", "stratum")
main_outputs <- file.path("data", "Sensitivity",
                          paste0("all_runs_", spatial_levels, "_1_T.rda"))
if (any(!file.exists(main_outputs))) {
  stop("Missing main-analysis outputs: ",
       paste(main_outputs[!file.exists(main_outputs)], collapse = ", "))
}
available_families <- families[vapply(families, function(fam) {
  any(grepl(paste0("^prey_category_[0-9]+", fam, "$"), prey_columns))
}, logical(1))]

for (fam in setdiff(families, available_families)) {
  warning("Missing prey-grouping columns for ", fam, "; skipping this family.", call. = FALSE)
}
rm(prey_columns)

# Alternative prey-grouping analyses
for (family in available_families) {
  for (level in spatial_levels) {
    local({
      SPATIAL_LEVEL <- level
      PREY_FAMILY <- family
      RUN_TAG <- ""
      N_MIN <- 5
      N_STRICT <- 25
      MIN_SETS <- 3
      R_PERM <- 999
      ALPHA <- 0.05
      BS_TEST <- "welch"
      COMP_RELATIVE <- TRUE
      OCC_BINARY <- FALSE
      ADD_NICHE_PART <- TRUE
      ADD_IND_NICHE <- FALSE
      source("03_Statistical_Pipeline.R", local = TRUE)
      run_all_scenarios(input_file, tag = RUN_TAG)
    })
  }
}

# Cross-scale comparisons
local({
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(sf); library(readr)
})

FAMILIES <- c("_1", available_families)
ALPHA    <- 0.05
OBJ      <- "res_biomass_PT"
OUT_DIR  <- "Sensitivity_comparison"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FAM_LAB <- c("_1" = "Pooled q = 1 (main)", "_2" = "Pooled q = 2", "_PP" = "Per predator")

load("data/dat_classed.rda")
xwalk_raw <- as.data.frame(dat_classed) %>%
  st_drop_geometry() %>%
  mutate(stratum = as.character(stratum)) %>%
  count(Area, stratum, name = "n_rec") %>%
  filter(n_rec > 0)
rm(dat_classed)

xwalk <- xwalk_raw %>%
  group_by(stratum) %>%
  mutate(n_stratum = sum(n_rec), share = n_rec / n_stratum, n_ecoregions = dplyr::n()) %>%
  slice_max(n_rec, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(stratum, eco_majority = Area, share_majority = round(share, 3),
            n_ecoregions, pure = n_ecoregions == 1)

split_strata <- xwalk_raw %>%
  group_by(stratum) %>%
  filter(dplyr::n() > 1) %>%
  mutate(share = round(100 * n_rec / sum(n_rec), 1)) %>%
  ungroup() %>%
  arrange(stratum, desc(n_rec))
write_csv(split_strata, file.path(OUT_DIR, "crosswalk_strata_ecoregions.csv"))

cat("\n--- Strata split across ecoregions ---\n"); print(split_strata, n = Inf)
cat("\n--- Crosswalk coverage ---\n")
print(xwalk %>% summarise(
  n_strata = dplyr::n(), n_pure = sum(pure), n_split = sum(!pure),
  pct_records_pure = round(100 * sum(xwalk_raw$n_rec[xwalk_raw$stratum %in% stratum[pure]]) /
                             sum(xwalk_raw$n_rec), 1)))

tidy_level <- function(res) {
  d <- res$results
  sp_col <- intersect(c("Area", "str"), names(d))[1]
  d$spatial_unit <- if (is.na(sp_col)) NA_character_ else as.character(d[[sp_col]])
  dplyr::select(d, -dplyr::any_of(c("Area", "str")))
}

load_level <- function(level, fam) {
  tag <- if (fam == "_1") "_T" else ""
  f <- paste0("data/Sensitivity/all_runs_", level, fam, tag, ".rda")
  if (!file.exists(f)) return(NULL)
  e <- new.env(); load(f, envir = e)
  if (!exists(OBJ, envir = e)) return(NULL)
  tidy_level(get(OBJ, envir = e))
}

emit <- function(df, stem, title) {
  write_csv(df, file.path(OUT_DIR, paste0(stem, ".csv")))
  cat("\n--- ", title, " ---\n", sep = ""); print(as.data.frame(df), row.names = FALSE)
  invisible(df)
}

compare_family <- function(fam) {
  L <- list(all_gulf = load_level("all_gulf", fam),
            ecoregion = load_level("ecoregion", fam),
            stratum   = load_level("stratum", fam))
  L <- L[!vapply(L, is.null, logical(1))]
  if (length(L) < 2) { message("Family ", fam, ": fewer than two levels found, skipped."); return(NULL) }
  cat("\n\n===================== FAMILY ", fam, " (", FAM_LAB[[fam]], ") =====================\n", sep = "")
  all3 <- bind_rows(L)

  emit(all3 %>% group_by(spatial_level) %>%
         summarise(n_cells = dplyr::n(), n_units = dplyr::n_distinct(spatial_unit),
                   n_resolutions = dplyr::n_distinct(x_threshold),
                   stomachs = sum(n_sto_P1 + n_sto_P2) / dplyr::n_distinct(x_threshold),
                   pct_testable = round(100 * mean(testable), 1),
                   pct_reliable = round(100 * mean(reliable), 1), .groups = "drop"),
       paste0(fam, "_levels_sample"), "Cells, units and stomach totals per level")

  comp <- all3 %>% filter(testable) %>% count(spatial_level, diagnostic) %>%
    group_by(spatial_level) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% ungroup()
  emit(comp %>% select(-n) %>% pivot_wider(names_from = spatial_level, values_from = pct, values_fill = 0),
       paste0(fam, "_levels_composition"), "Diagnostic composition per level (testable cells)")

  emit(all3 %>% filter(testable) %>% group_by(spatial_level) %>%
         summarise(n = dplyr::n(),
                   k_med = median(pmin(n_set_P1, n_set_P2), na.rm = TRUE),
                   BC_med = round(median(BC, na.rm = TRUE), 3),
                   R2_med = round(median(R2_comp, na.rm = TRUE), 3),
                   pct_H_sig = round(100 * mean(p_H < ALPHA, na.rm = TRUE), 1),
                   pct_Bs_sig = round(100 * mean(p_Bs < ALPHA, na.rm = TRUE), 1),
                   pct_comp_sig = round(100 * mean(p_comp < ALPHA, na.rm = TRUE), 1), .groups = "drop"),
       paste0(fam, "_levels_effect_vs_power"), "Effect size vs detection rate, by level")

  emit(all3 %>% filter(testable) %>% group_by(spatial_level) %>%
         summarise(pct_comp_sig = round(100 * mean(p_comp < ALPHA, na.rm = TRUE), 1),
                   pct_disp_sig = round(100 * mean(p_disp < ALPHA, na.rm = TRUE), 1),
                   pct_confounded = round(100 * mean(comp_dispersion, na.rm = TRUE), 1), .groups = "drop"),
       paste0(fam, "_levels_confounded"), "Confounded location/dispersion cells, by level")

  if (all(c("ecoregion", "stratum") %in% names(L))) {
    key <- c("x_threshold", "species", "size_class")
    L3 <- L$stratum %>% filter(testable) %>% left_join(xwalk, by = c("spatial_unit" = "stratum"))
    cat("\n--- Coverage of the ecoregion / stratum comparison ---\n")
    print(L3 %>% summarise(n_stratum_cells = dplyr::n(), n_kept = sum(pure, na.rm = TRUE),
                           pct_kept = round(100 * mean(pure, na.rm = TRUE), 1)))
    L3_pure <- L3 %>% filter(pure)
    L2 <- L$ecoregion %>% filter(testable) %>%
      select(all_of(key), eco_unit = spatial_unit, d_eco = diagnostic)

    emit(inner_join(L2, L3_pure %>% select(all_of(key), eco_unit = eco_majority, d_str = diagnostic),
                    by = c(key, "eco_unit"), relationship = "many-to-many") %>%
           count(d_eco, d_str) %>% group_by(d_eco) %>%
           mutate(pct_within_eco = round(100 * n / sum(n), 1)) %>% ungroup() %>% arrange(d_eco, desc(n)),
         paste0(fam, "_eco_vs_strata"), "Ecoregional verdict vs stratum verdicts (unambiguous strata)")

    emit(inner_join(L2, L3_pure %>% select(all_of(key), eco_unit = eco_majority, d_str = diagnostic),
                    by = c(key, "eco_unit"), relationship = "many-to-many") %>%
           count(across(all_of(c(key, "eco_unit"))), d_eco, d_str) %>%
           group_by(across(all_of(c(key, "eco_unit")))) %>%
           slice_max(n, n = 1, with_ties = FALSE) %>% ungroup() %>%
           summarise(n_eco_cells = dplyr::n(), pct_match = round(100 * mean(d_eco == d_str), 1)),
         paste0(fam, "_eco_vs_strata_match"), "Ecoregion vs modal stratum verdict")
  }

  p <- comp %>%
    mutate(spatial_level = factor(spatial_level, levels = c("all_gulf", "ecoregion", "stratum"))) %>%
    ggplot(aes(x = reorder(diagnostic, pct), y = pct, fill = spatial_level)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.9) +
    coord_flip() +
    scale_fill_brewer(palette = "Set2", name = "Spatial level") +
    labs(title = paste0("Diagnostic composition across spatial levels - ", FAM_LAB[[fam]]),
         subtitle = "Biomass, 2004-2006 vs 2018-2019, testable cells only",
         x = NULL, y = "Within-level frequency (%)") +
    theme_minimal() +
    theme(axis.text.y = element_text(face = "bold"), legend.position = "bottom")
  ggsave(file.path(OUT_DIR, paste0(fam, "_diagnostic_composition_by_level.png")), p,
         width = 10, height = 6, dpi = 300)

  comp %>% mutate(family = fam)
}

comp_all <- bind_rows(lapply(FAMILIES, compare_family))

if (nrow(comp_all)) {
  write_csv(comp_all, file.path(OUT_DIR, "families_levels_composition.csv"))
  p_all <- comp_all %>%
    mutate(spatial_level = factor(spatial_level, levels = c("all_gulf", "ecoregion", "stratum")),
           family = factor(FAM_LAB[family], levels = unname(FAM_LAB))) %>%
    ggplot(aes(x = reorder(diagnostic, pct), y = pct, fill = spatial_level)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.9) +
    coord_flip() +
    facet_wrap(~family, nrow = 1) +
    scale_fill_brewer(palette = "Set2", name = "Spatial level") +
    labs(title = "Diagnostic composition across spatial levels, by prey-grouping family",
         subtitle = "Biomass, 2004-2006 vs 2018-2019, testable cells only",
         x = NULL, y = "Within-level frequency (%)") +
    theme_minimal() +
    theme(axis.text.y = element_text(face = "bold"), legend.position = "bottom",
          strip.text = element_text(face = "bold"))
  ggsave(file.path(OUT_DIR, "families_diagnostic_composition_by_level.png"), p_all,
         width = 4.5 * dplyr::n_distinct(comp_all$family) + 2, height = 6, dpi = 300)
}

cat("\nCross-level comparison written to ", normalizePath(OUT_DIR), "\n", sep = "")

})

# Prey-grouping comparisons
local({
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

PREY_FAMILY <- "_1"
source("R_helpers/Results.R", local = TRUE)

MAIN_PROJECT <- normalizePath(getwd())

ALPHA  <- 0.05
LEVELS <- SPATIAL_LEVELS_ORD
P_COLS <- c("p_comp", "p_H", "p_Bs", "p_disp")
TEST_LAB <- c(p_comp = "Composition (PERMANOVA)", p_H = "Diversity H'",
              p_Bs = "Niche breadth Bs", p_disp = "Dispersion (PERMDISP)")

APPROACHES <- list(
  pooled_q1    = list(root = MAIN_PROJECT, suffix = "_1"),
  pooled_q2    = list(root = MAIN_PROJECT, suffix = "_2"),
  per_predator = list(root = MAIN_PROJECT, suffix = "_PP")
)
APPROACHES <- APPROACHES[vapply(APPROACHES, function(a) a$suffix %in% c("_1", available_families), logical(1))]
APPROACH_LAB <- c(pooled_q1    = "Pooled grouping, q = 1 (main analysis)",
                  pooled_q2    = "Pooled grouping, q = 2",
                  per_predator = "Per-predator grouping")
APPROACH_PAL <- c(pooled_q1 = "#2F4A5A", pooled_q2 = "#7F9CAB", per_predator = "#D96C4A")

SCEN_MAP <- data.table(
  period_1 = c("2004", "2004-2005", "2018", "2006", "2004-2006"),
  period_2 = c("2006", "2006",      "2019", "2018", "2018-2019"),
  contrast = c("P1a",  "P1b",       "P2",   "PTb",  "PTa")
)
FAM_MAP <- data.table(diagnostic = names(STATE_TO_FAMILY),
                      family     = unname(STATE_TO_FAMILY))

save_app <- function(p, stem, w, h) save_fig(p, stem, w, h, dir = DIR_APPEND)

.is_results <- function(x)
  is.data.frame(x) && all(c("x_threshold", "p_comp", "mode") %in% names(x))

.collect_results <- function(x) {
  if (.is_results(x)) return(list(as.data.table(x)))
  if (is.list(x) && !is.data.frame(x))
    return(unlist(lapply(x, .collect_results), recursive = FALSE))
  list()
}

.results_from_rda <- function(f) {
  e <- new.env(); load(f, envir = e)
  unlist(lapply(ls(e), function(o) .collect_results(get(o, envir = e))),
         recursive = FALSE)
}

load_approach <- function(root, suffix, spatial) {
  if (!dir.exists(root)) return(NULL)
  tag <- if (suffix == "_1") "_T" else ""
  pooled <- file.path(root, "data", "Sensitivity",
                      sprintf("all_runs_%s%s%s.rda", spatial, suffix, tag))
  if (file.exists(pooled)) {
    tabs <- .results_from_rda(pooled)
    if (length(tabs)) return(list(dt = rbindlist(tabs, fill = TRUE), src = pooled))
  }
  d <- file.path(root, paste0("Sensitivity_", spatial, suffix, tag))
  if (dir.exists(d)) {
    fs <- list.files(d, pattern = "\\.rda$", full.names = TRUE)
    if (length(fs)) {
      tabs <- unlist(lapply(fs, .results_from_rda), recursive = FALSE)
      if (length(tabs)) return(list(dt = rbindlist(tabs, fill = TRUE),
                                    src = paste0(d, " (", length(fs), " files)")))
    }
  }
  NULL
}

period_bands <- function(n_within, n_total) {
  list(
    annotate("rect", xmin = 0.5, xmax = n_within + 0.5, ymin = -Inf, ymax = Inf,
             fill = "#D7E6EC", alpha = 0.55),
    annotate("rect", xmin = n_within + 0.5, xmax = n_total + 0.5,
             ymin = -Inf, ymax = Inf, fill = "#F6E0D6", alpha = 0.65),
    geom_vline(xintercept = n_within + 0.5, linetype = "dashed",
               colour = "grey60", linewidth = 0.3)
  )
}

for (SPATIAL in LEVELS) {

  cat("\n=============", SPATIAL, "=============\n")

  res_list <- list()
  for (ap in names(APPROACHES)) {
    r <- load_approach(APPROACHES[[ap]]$root, APPROACHES[[ap]]$suffix, SPATIAL)
    if (!is.null(r)) {
      dt <- merge(r$dt[, approach := ap], SCEN_MAP,
                  by = c("period_1", "period_2"), all.x = TRUE)
      res_list[[ap]] <- dt
      xs <- sort(unique(dt$x_threshold))
      message(sprintf("[ok]     %-13s %s | %d rows | X = %d-%d",
                      ap, r$src, nrow(dt), min(xs), max(xs)))
    } else message("[absent] ", ap)
  }
  if (length(res_list) < 2) {
    message("Fewer than two approaches available for ", SPATIAL, " -> skipped")
    next
  }

  both <- rbindlist(res_list, fill = TRUE)
  both[, approach := factor(approach, levels = names(APPROACHES))]
  both[, contrast := factor(contrast, levels = CONTRAST_LEVELS)]

  common_x <- Reduce(intersect, lapply(res_list, function(d) unique(d$x_threshold)))
  sig_x    <- if (length(common_x)) common_x else unique(both$x_threshold)
  message(length(sig_x), " common resolutions used (",
          min(sig_x), "-", max(sig_x), ")")

  long <- melt(both[testable == TRUE & x_threshold %in% sig_x & !is.na(contrast),
                    c("approach", "contrast", "mode", "x_threshold",
                      intersect(P_COLS, names(both))), with = FALSE],
               id.vars = c("approach", "contrast", "mode", "x_threshold"),
               variable.name = "test", value.name = "p")

  pct_by_x <- long[!is.na(p), .(pct_sig = 100 * mean(p < ALPHA), n_cells = .N),
                   by = .(approach, contrast, mode, test, x_threshold)]

  sig_summary <- pct_by_x[, .(pct_mean = round(mean(pct_sig), 1),
                              pct_sd = round(sd(pct_sig), 1),
                              n_thresholds = uniqueN(x_threshold),
                              n_cells_mean = round(mean(n_cells), 1)),
                          by = .(approach, contrast, mode, test)]
  setorder(sig_summary, test, contrast, mode, approach)

  fwrite(sig_summary[, .(approach = APPROACH_LAB[as.character(approach)],
                         contrast = contrast_label(contrast),
                         currency = CURRENCY_LAB[as.character(mode)],
                         test = TEST_LAB[as.character(test)],
                         pct_mean, pct_sd, n_thresholds, n_cells_mean)],
         file.path(DIR_APPEND, paste0("TableB2_grouping_pct_significant_", SPATIAL, ".csv")))

  pB2 <- ggplot(sig_summary, aes(x = contrast, y = pct_mean, fill = approach)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = pmax(pct_mean - pct_sd, 0),
                      ymax = pmin(pct_mean + pct_sd, 100)),
                  position = position_dodge(width = 0.8), width = 0.25) +
    facet_grid(test ~ mode, labeller = labeller(test = TEST_LAB, mode = CURRENCY_LAB)) +
    scale_x_discrete(labels = CONTRAST_LAB) +
    scale_fill_manual(values = APPROACH_PAL, labels = APPROACH_LAB, name = NULL) +
    labs(title = paste0("Share of testable cells significant (p < ", ALPHA,
                        ") by prey-grouping approach, ", LEVEL_SHORT[[SPATIAL]], " scale"),
         subtitle = paste0("Mean +/- SD across ", length(sig_x),
                           " common taxonomic resolutions"),
         x = NULL, y = "Cells significant (%)") +
    theme_diag(base_size = 10)
  save_app(pB2, paste0("FigB2_grouping_pct_significant_", SPATIAL), 10, 11)

  db <- both[testable == TRUE & x_threshold %in% sig_x &
               !is.na(diagnostic) & diagnostic != INCONCLUSIVE_LAB &
               !is.na(contrast)]
  db <- merge(db, FAM_MAP, by = "diagnostic", all.x = TRUE)
  if (anyNA(db$family))
    warning("Diagnostics outside STATE_TO_FAMILY: ",
            paste(unique(db[is.na(family), diagnostic]), collapse = ", "))
  db <- db[!is.na(family)]

  keys <- c("approach", "mode", "contrast", "x_threshold")
  tot  <- db[, .(n_tot = .N), by = keys]
  cnt  <- db[, .(n = .N), by = c(keys, "family")]
  grid <- tot[, CJ(family = FAMILY_LEVELS), by = keys]
  fam_by_x <- merge(grid, cnt, by = c(keys, "family"), all.x = TRUE)
  fam_by_x <- merge(fam_by_x, tot, by = keys)
  fam_by_x[is.na(n), n := 0L]
  fam_by_x[, pct := 100 * n / n_tot]

  fam_summary <- fam_by_x[, .(pct_mean = round(mean(pct), 1),
                              pct_sd   = round(sd(pct), 1),
                              n_thresholds = uniqueN(x_threshold)),
                          by = .(approach, mode, contrast, family)]
  fam_summary[, family := factor(family, levels = FAMILY_LEVELS)]
  setorder(fam_summary, approach, mode, contrast, family)

  out_long <- fam_summary[, .(approach = APPROACH_LAB[as.character(approach)],
                              currency = CURRENCY_LAB[as.character(mode)],
                              contrast = contrast_label(contrast),
                              family, pct_mean, pct_sd, n_thresholds)]
  fwrite(out_long, file.path(DIR_APPEND, paste0("TableB1_grouping_families_", SPATIAL, ".csv")))

  fam_wide <- dcast(out_long, approach + currency + contrast ~ family, value.var = "pct_mean")
  setcolorder(fam_wide, c("approach", "currency", "contrast", FAMILY_LEVELS))
  fwrite(fam_wide, file.path(DIR_APPEND, paste0("TableB1_grouping_families_wide_", SPATIAL, ".csv")))

  dplot <- copy(fam_summary)[, x := as.integer(contrast)]
  present  <- levels(droplevels(dplot$contrast))
  n_within <- sum(INTRA_CONTRASTS %in% present)
  n_total  <- length(CONTRAST_LEVELS)

  pB1 <- ggplot(dplot, aes(x = x, y = pct_mean, colour = approach, group = approach)) +
    period_bands(n_within, n_total) +
    geom_errorbar(aes(ymin = pmax(pct_mean - pct_sd, 0), ymax = pct_mean + pct_sd),
                  width = 0.12, linewidth = 0.35, alpha = 0.8) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = seq_len(n_total), labels = unname(CONTRAST_LAB),
                       limits = c(0.5, n_total + 0.5), expand = c(0, 0)) +
    scale_colour_manual(values = APPROACH_PAL, labels = APPROACH_LAB, name = NULL) +
    facet_grid(family ~ mode, scales = "free_y",
               labeller = labeller(mode = CURRENCY_LAB)) +
    labs(title = paste0("Diagnostic families by prey-grouping approach, ",
                        LEVEL_SHORT[[SPATIAL]], " scale"),
         subtitle = paste0("Mean +/- SD across ", length(sig_x),
                           " common taxonomic resolutions.\nShaded bands: ",
                           "within-period (left) and between-period (right) contrasts."),
         x = NULL, y = "Frequency (%)") +
    theme_diag(base_size = 10)
  save_app(pB1, paste0("FigB1_grouping_families_", SPATIAL), 9, 10)

  message("Written for ", SPATIAL, ": FigB1, FigB2, TableB1 (long + wide), TableB2")
}

cat("\nDone. Appendix outputs in ", normalizePath(DIR_APPEND), "\n", sep = "")

})
