# Derived results and summary figures

if (!file.exists("data/dat_classed.rda")) stop("Missing prepared data")

# Spatial family summaries
local({
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(forcats)
  library(sf); library(readr); library(stringr)
})

PREY_FAMILY <- "_1"
RDA_DIRS <- c(
  all_gulf  = "Sensitivity_all_gulf_1_T",
  ecoregion = "Sensitivity_ecoregion_1_T",
  stratum   = "Sensitivity_stratum_1_T"
)
source("R_helpers/Results.R", local = TRUE)

N_FLAG        <- 5
STRATA_PATH   <- "strata_rv_gulf.rds"
COAST_PATH    <- "gulf.coast.intermediate.csv"
MAP_XLIM      <- c(-66.2, -60.0)
MAP_YLIM      <- c(45.5, 49.2)
MAP_CONTRASTS <- c("P1a", "P1b", "P2", "PTa")

cat("\nReading pipeline results...\n")
res_all <- read_all_runs()

cat("Runs read:\n")
print(res_all %>% count(level, currency, contrast) %>% arrange(level, currency, contrast),
      n = Inf)

fam_gulf <- res_all %>%
  filter(level == "all_gulf") %>%
  freq_by_resolution(keys = c("currency", "contrast")) %>%
  add_families(keys = c("currency", "contrast"))

fam_unit <- res_all %>%
  filter(level != "all_gulf") %>%
  freq_by_resolution(keys = c("currency", "contrast", "level", "spatial_unit")) %>%
  add_families(keys = c("currency", "contrast", "level", "spatial_unit")) %>%
  left_join(
    n_units_by(filter(res_all, level != "all_gulf"),
               keys = c("currency", "contrast", "level", "spatial_unit")),
    by = c("currency", "contrast", "level", "spatial_unit")
  )

fam_level <- bind_rows(
  fam_gulf %>% mutate(level = "all_gulf", .before = 1),
  fam_unit %>%
    group_by(currency, contrast, level) %>%
    summarise(across(all_of(FAMILY_LEVELS), ~mean(.x, na.rm = TRUE)),
              .groups = "drop")
) %>%
  mutate(level = factor(level, levels = SPATIAL_LEVELS_ORD))

long_of <- function(d) {
  d %>%
    pivot_longer(all_of(FAMILY_LEVELS), names_to = "family", values_to = "pct") %>%
    mutate(family   = factor(family, levels = FAMILY_LEVELS),
           currency = factor(currency, levels = names(CURRENCY_LAB)))
}

FAM_PAL_N  <- if (is.null(names(FAMILY_PAL))) setNames(FAMILY_PAL, FAMILY_LEVELS) else FAMILY_PAL
FAMILY_SHAPE_FILL <- setNames(c(21, 22, 23, 24), FAMILY_LEVELS)
X_OFFSET   <- 0.12

fill_of <- function(family, currency)
  ifelse(currency == "biomass", unname(FAM_PAL_N[as.character(family)]), "white")

family_currency_scales <- function() list(
  scale_colour_manual(values = FAM_PAL_N, name = NULL),
  scale_shape_manual(values = FAMILY_SHAPE_FILL, name = NULL),
  scale_fill_identity(),
  scale_alpha_manual(values = c(biomass = 1, occurrence = 1),
                     labels = CURRENCY_LAB, name = NULL),
  guides(
    colour = guide_legend(order = 1,
                          override.aes = list(fill = unname(FAM_PAL_N[FAMILY_LEVELS]))),
    shape  = guide_legend(order = 1,
                          override.aes = list(fill = unname(FAM_PAL_N[FAMILY_LEVELS]))),
    alpha  = guide_legend(order = 2,
                          override.aes = list(shape = 21, colour = "grey25", size = 2.8,
                                              fill = c("grey25", "white"))))
)

cat("\nFigure 5: Gulf-wide regime break\n")

fam_gulf_sd <- spread_by_resolution(filter(res_all, level == "all_gulf"),
                                    keys = c("currency", "contrast")) %>%
  select(currency, contrast, family, sd = pct_sd) %>%
  mutate(currency = factor(currency, levels = names(CURRENCY_LAB)))

d5 <- long_of(fam_gulf) %>%
  filter(!is.na(contrast)) %>%
  left_join(fam_gulf_sd, by = c("currency", "contrast", "family")) %>%
  mutate(x = as.integer(contrast) + ifelse(currency == "biomass", -X_OFFSET, X_OFFSET),
         fill_col = fill_of(family, currency))

fig5 <- ggplot(d5, aes(x, pct, colour = family, shape = family)) +
  decade_bands() +
  annotate("text", x = 2,   y = Inf, label = "WITHIN PERIOD",
           colour = "#2F4A5A", fontface = 2, size = 3.6, vjust = 1.8) +
  annotate("text", x = 4.5, y = Inf, label = "BETWEEN PERIODS",
           colour = "#9B2226", fontface = 2, size = 3.6, vjust = 1.8) +
  geom_errorbar(aes(ymin = pmax(pct - sd, 0), ymax = pmin(pct + sd, 100)),
                width = 0.08, linewidth = 0.4, alpha = 0.7, show.legend = FALSE) +
  geom_point(aes(fill = fill_col, alpha = currency), size = 2.8, stroke = 0.9) +
  scale_x_continuous(breaks = seq_along(CONTRAST_LEVELS),
                     labels = unname(CONTRAST_LAB[CONTRAST_LEVELS]),
                     expand = expansion(add = 0.25)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  family_currency_scales() +
  labs(
    title = "Diet regime break at the Gulf-wide scale",
    subtitle = paste0("Frequency of the four diagnostic families among predator x ",
                      "size-class cells: mean and\nstandard deviation across ~100 prey taxonomic resolutions"),
    x = NULL, y = "Frequency (%)"
  ) +
  theme_diag() +
  theme(legend.box = "horizontal")

save_fig(fig5, "Fig5_regime_break_gulf", width = 9, height = 5.6)

cat("\nFigure 6: by ecoregion\n")

d6 <- long_of(filter(fam_unit, level == "ecoregion")) %>%
  filter(!is.na(contrast)) %>%
  mutate(unit_lab = area_label(spatial_unit))
d6_n <- d6 %>% distinct(currency, unit_lab, contrast, n_units)

fig6 <- ggplot(d6, aes(contrast, pct, fill = family)) +
  geom_col(width = 0.72, colour = "grey25", linewidth = 0.2) +
  geom_text(data = d6_n, aes(contrast, 101, label = paste0("n=", n_units)),
            inherit.aes = FALSE, size = 2.3, vjust = 0, colour = "grey30") +
  scale_x_discrete(labels = CONTRAST_LAB) +
  geom_hline(yintercept = c(25, 50, 75), colour = "white",
             linewidth = 0.3, linetype = "dotted") +
  facet_grid(currency ~ unit_lab,
             labeller = labeller(currency = CURRENCY_LAB)) +
  scale_fill_manual(values = FAMILY_PAL, name = NULL) +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Diet regime break by ecoregion",
    subtitle = paste0("Composition of the four families within each ecoregion. ",
                      "n = predator x size-class units contributing to each bar."),
    x = NULL, y = "Within-ecoregion frequency (%)"
  ) +
  theme_diag() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5),
        panel.grid.major.x = element_blank())

save_fig(fig6, "Fig6_regime_break_ecoregion", width = 11, height = 6.4)

cat("\nFigures 7-8: ecoregion and stratum scales\n")

rb_sd <- spread_by_resolution_units(filter(res_all, level != "all_gulf"),
                                    keys = c("currency", "contrast", "level")) %>%
  transmute(currency = factor(currency, levels = names(CURRENCY_LAB)),
            contrast, level = as.character(level), family, sd = pct_sd)

make_regime_break <- function(lvl, stem, title_txt, subtitle_txt) {
  d <- long_of(filter(fam_level, level == lvl)) %>%
    filter(!is.na(contrast)) %>%
    mutate(level = as.character(level)) %>%
    left_join(rb_sd, by = c("currency", "contrast", "level", "family")) %>%
    mutate(x = as.integer(contrast) + ifelse(currency == "biomass", -X_OFFSET, X_OFFSET),
           fill_col = fill_of(family, currency))

  g <- ggplot(d, aes(x, pct, colour = family, shape = family)) +
    decade_bands() +
    annotate("text", x = 2,   y = Inf, label = "WITHIN PERIOD",
             colour = "#2F4A5A", fontface = 2, size = 3.6, vjust = 1.8) +
    annotate("text", x = 4.5, y = Inf, label = "BETWEEN PERIODS",
             colour = "#9B2226", fontface = 2, size = 3.6, vjust = 1.8) +
    geom_errorbar(aes(ymin = pmax(pct - sd, 0), ymax = pmin(pct + sd, 100)),
                  width = 0.08, linewidth = 0.4, alpha = 0.7, show.legend = FALSE) +
    geom_point(aes(fill = fill_col, alpha = currency), size = 2.8, stroke = 0.9) +
    scale_x_continuous(breaks = seq_along(CONTRAST_LEVELS),
                       labels = unname(CONTRAST_LAB[CONTRAST_LEVELS]),
                       expand = expansion(add = 0.25)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    family_currency_scales() +
    labs(title = title_txt, subtitle = subtitle_txt, x = NULL, y = "Frequency (%)") +
    theme_diag() +
    theme(legend.box = "horizontal")

  save_fig(g, stem, width = 9, height = 5.6)
}

make_regime_break(
  "ecoregion", "Fig7_regime_break_ecoregion",
  "Diet regime break at the ecoregion scale",
  paste0("Family frequencies averaged across the four ecoregions: mean and\n",
         "standard deviation across ~100 prey taxonomic resolutions"))

make_regime_break(
  "stratum", "Fig8_regime_break_stratum",
  "Diet regime break at the stratum scale",
  paste0("Family frequencies averaged across survey strata: mean and\n",
         "standard deviation across ~100 prey taxonomic resolutions"))

cat("\nFigure S3: regime break by ecoregion\n")

unit_sd <- spread_by_resolution(filter(res_all, level != "all_gulf"),
                                keys = c("currency", "contrast", "level", "spatial_unit")) %>%
  transmute(currency = factor(currency, levels = names(CURRENCY_LAB)),
            contrast, level = as.character(level), spatial_unit, family, sd = pct_sd)

make_unit_regime_break <- function(lvl, stem, title_txt, ncol_panels, width, height,
                                   base_size = 10, point_size = 1.8) {
  d <- long_of(filter(fam_unit, level == lvl)) %>%
    filter(!is.na(contrast)) %>%
    mutate(level = as.character(level)) %>%
    left_join(unit_sd, by = c("currency", "contrast", "level", "spatial_unit", "family")) %>%
    mutate(unit_lab = if (lvl == "ecoregion") area_label(spatial_unit) else paste("Stratum", spatial_unit),
           x = as.integer(contrast) + ifelse(currency == "biomass", -X_OFFSET, X_OFFSET),
           fill_col = fill_of(family, currency))
  n_cells <- d %>% group_by(unit_lab) %>% summarise(n = max(n_units, na.rm = TRUE), .groups = "drop")
  d <- d %>% left_join(n_cells, by = "unit_lab") %>%
    mutate(unit_lab = paste0(unit_lab, " (n = ", n, ")"))
  if (lvl == "stratum") d <- d %>% mutate(unit_lab = fct_reorder(unit_lab, as.numeric(spatial_unit)))

  g <- ggplot(d, aes(x, pct, colour = family, shape = family)) +
    decade_bands() +
    geom_errorbar(aes(ymin = pmax(pct - sd, 0), ymax = pmin(pct + sd, 100)),
                  width = 0.1, linewidth = 0.3, alpha = 0.7, show.legend = FALSE) +
    geom_point(aes(fill = fill_col, alpha = currency), size = point_size, stroke = 0.8) +
    facet_wrap(~unit_lab, ncol = ncol_panels) +
    scale_x_continuous(breaks = seq_along(CONTRAST_LEVELS),
                       labels = unname(CONTRAST_LAB[CONTRAST_LEVELS]),
                       expand = expansion(add = 0.25)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
    family_currency_scales() +
    labs(title = title_txt,
         subtitle = paste0("Family frequencies within each unit: mean and SD across ~100 prey resolutions.\n",
                           "Bands: within-period (left) and between-period (right) contrasts; ",
                           "n = predator x size-class cells."),
         x = NULL, y = "Frequency (%)") +
    theme_diag(base_size = base_size) +
    theme(legend.position = "bottom", legend.box = "horizontal",
          axis.text.x = element_text(size = base_size - 3),
          strip.text = element_text(size = base_size - 1))
  save_fig(g, stem, width = width, height = height)
}

if (any(fam_unit$level == "ecoregion"))
  make_unit_regime_break("ecoregion", "FigS3_regime_break_by_ecoregion",
                         "Diet regime break within each ecoregion",
                         ncol_panels = 2, width = 11, height = 8, base_size = 11, point_size = 2.2)

cat("\nFigure 10: cross-scale gradient\n")

grad_sd <- bind_rows(
  spread_by_resolution(filter(res_all, level == "all_gulf"),
                       keys = c("currency", "contrast")) %>%
    mutate(level = "all_gulf"),
  spread_by_resolution_units(filter(res_all, level != "all_gulf"),
                             keys = c("currency", "contrast", "level"))
) %>%
  transmute(currency = factor(currency, levels = names(CURRENCY_LAB)),
            contrast, level = as.character(level), family, sd = pct_sd)

d_grad <- long_of(fam_level) %>%
  filter(contrast == "PTa") %>%
  mutate(level = as.character(level)) %>%
  left_join(grad_sd, by = c("currency", "contrast", "level", "family")) %>%
  mutate(level = factor(level, levels = rev(SPATIAL_LEVELS_ORD)))

fig_grad <- ggplot(d_grad, aes(level, pct, colour = family, group = family)) +
  geom_errorbar(aes(ymin = pmax(pct - sd, 0), ymax = pmin(pct + sd, 100)),
                width = 0.1, linewidth = 0.4, alpha = 0.7, show.legend = FALSE) +
  geom_line(linewidth = 0.9) +
  geom_point(aes(shape = family), size = 3.2) +
  geom_text(aes(label = sprintf("%.1f", pct)), vjust = -1.1,
            size = 2.9, show.legend = FALSE) +
  facet_wrap(~currency, labeller = labeller(currency = CURRENCY_LAB)) +
  scale_colour_manual(values = FAMILY_PAL, name = NULL) +
  scale_shape_manual(values = FAMILY_SHAPE, name = NULL) +
  scale_x_discrete(labels = c(stratum = "Stratum", ecoregion = "Ecoregion",
                              all_gulf = "Gulf-wide")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    title = paste0("Effect of spatial aggregation on the diagnostic composition (",
                   CONTRAST_LAB_1L[["PTa"]], ")"),
    subtitle = paste0("Stratum and ecoregion values are means across units; the ",
                      "Gulf-wide value is a single pooled test.\nBars = SD across resolutions. Aggregation ",
                      "raises apparent Reorganisation and lowers Stability, ",
                      "while Substitution is scale-invariant."),
    x = "Spatial unit, from finest to coarsest", y = "Frequency (%)"
  ) +
  theme_diag()

save_fig(fig_grad, "Fig10_crossscale_gradient", width = 9.5, height = 5.4)

cat("\nFigures 9 and S4: stratum maps\n")

strata <- load_strata(STRATA_PATH)

if (is.null(strata)) {
  message("Stratum geometry not found (neither gulf.spatial nor ", STRATA_PATH,
          "). Maps skipped; every other figure is written.")
} else {
  strata <- sf::st_make_valid(strata)
  if (is.na(sf::st_crs(strata))) sf::st_crs(strata) <- 4326

  lab_xy <- if (all(c("label_x", "label_y") %in% names(strata))) {
    strata %>% sf::st_drop_geometry() %>% transmute(spatial_unit = as.character(str),
                                                    X = label_x, Y = label_y)
  } else {
    cc <- suppressWarnings(sf::st_coordinates(sf::st_centroid(strata)))
    tibble(spatial_unit = as.character(strata$str), X = cc[, 1], Y = cc[, 2])
  }

  coast <- if (file.exists(COAST_PATH)) {
    read_csv(COAST_PATH, show_col_types = FALSE) %>% filter(!is.na(x), !is.na(y))
  } else NULL

  map_dat <- long_of(filter(fam_unit, level == "stratum")) %>%
    filter(!is.na(contrast), contrast %in% MAP_CONTRASTS) %>%
    mutate(low_n = !is.na(n_units) & n_units < N_FLAG)

  strata$spatial_unit <- as.character(strata$str)

  render_map <- function(currencies, facet_formula, title, stem,
                         ncol_panels, nrow_panels, mask_low_n = FALSE) {
    d <- strata %>%
      left_join(filter(map_dat, currency %in% currencies), by = "spatial_unit") %>%
      filter(!is.na(pct)) %>%
      mutate(contrast = droplevels(factor(contrast, levels = CONTRAST_LEVELS)),
             currency = droplevels(factor(currency, levels = names(CURRENCY_LAB))))
    d_col  <- if (mask_low_n) filter(d, !low_n) else d
    d_grey <- if (mask_low_n) filter(d, low_n) else d[0, ]

    lab <- d_col %>%
      sf::st_drop_geometry() %>%
      left_join(lab_xy, by = "spatial_unit") %>%
      transmute(family, contrast, currency, X, Y,
                lab = ifelse(low_n & !mask_low_n, paste0(sprintf("%.0f", pct), "*"),
                             sprintf("%.0f", pct)),
                txt = ifelse(pct < 50, "grey95", "grey10"))

    subtitle <- if (mask_low_n) {
      paste0("Grey strata: fewer than ", N_FLAG, " predator x size-class units ",
             "(value not shown). 2006 vs 2018 is not estimable at this scale.")
    } else {
      paste0("Every stratum shown; * = built on < ", N_FLAG,
             " predator x size-class units (interpret with caution). ",
             "2006 vs 2018 is not estimable at this scale.")
    }

    g <- ggplot()
    if (!is.null(coast)) {
      g <- g + geom_polygon(data = coast, aes(x, y, group = pid),
                            fill = "grey88", colour = NA)
    }
    g <- g +
      geom_sf(data = strata, fill = "grey97", colour = "grey70", linewidth = 0.15) +
      geom_sf(data = d_col, aes(fill = pct), colour = "grey45", linewidth = 0.15)
    if (nrow(d_grey)) {
      g <- g + geom_sf(data = d_grey, fill = "grey78", colour = "grey45", linewidth = 0.15)
    }
    g <- g +
      geom_text(data = lab, aes(X, Y, label = lab, colour = txt),
                size = 2.1, fontface = "bold") +
      scale_colour_identity() +
      facet_grid(facet_formula,
                 labeller = labeller(contrast = CONTRAST_LAB_1L,
                                     currency = CURRENCY_LAB)) +
      scale_fill_viridis_c(option = "viridis", limits = c(0, 100),
                           breaks = seq(0, 100, 25),
                           name = "Within-stratum\nfrequency (%)") +
      scale_x_continuous(breaks = seq(-65, -61, by = 1)) +
      scale_y_continuous(breaks = seq(46, 49, by = 1)) +
      coord_sf(xlim = MAP_XLIM, ylim = MAP_YLIM, expand = FALSE) +
      labs(title = title, subtitle = subtitle) +
      theme_diag(base_size = 11) +
      theme(panel.grid   = element_line(colour = "grey92"),
            panel.spacing.x = unit(0.8, "lines"),
            axis.title   = element_blank(),
            axis.text    = element_text(size = 6),
            strip.text   = element_text(face = "bold"),
            strip.text.y = element_text(angle = 0),
            plot.title   = element_text(face = "bold"))

    save_fig(g, stem,
             width  = 2.55 * ncol_panels + 1.8,
             height = 2.45 * nrow_panels + 1.2,
             dpi = 200)
  }

  n_ct <- length(intersect(MAP_CONTRASTS, levels(droplevels(map_dat$contrast))))
  n_fm <- length(FAMILY_LEVELS)

  render_map("biomass", family ~ contrast,
             "Diagnostic families across survey strata - biomass",
             "Fig9_map_strata_biomass", n_ct, n_fm)

  render_map("occurrence", family ~ contrast,
             "Diagnostic families across survey strata - occurrence",
             "FigS4_map_strata_occurrence", n_ct, n_fm)

  render_map("biomass", family ~ contrast,
             "Diagnostic families across survey strata - biomass, low-n strata masked",
             "FigS5_map_strata_biomass_masked", n_ct, n_fm, mask_low_n = TRUE)

  render_map("occurrence", family ~ contrast,
             "Diagnostic families across survey strata - occurrence, low-n strata masked",
             "FigS6_map_strata_occurrence_masked", n_ct, n_fm, mask_low_n = TRUE)
}

export_csv <- function(x, name) {
  x %>%
    mutate(contrast = contrast_label(contrast)) %>%
    write_csv(file.path(DIR_FIGURES, name))
}
export_csv(fam_gulf,  "data_Fig5_gulf_families.csv")
export_csv(fam_unit,  "data_Fig6_9_unit_families.csv")
export_csv(fam_level, "data_Fig10_crossscale.csv")

cat("\nFigures and their source tables written to ", normalizePath(DIR_FIGURES), "\n", sep = "")

})

# Analytical summary tables
local({
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
})

PREY_FAMILY <- "_1"
RDA_DIRS <- c(
  all_gulf  = "Sensitivity_all_gulf_1_T",
  ecoregion = "Sensitivity_ecoregion_1_T",
  stratum   = "Sensitivity_stratum_1_T"
)
source("R_helpers/Results.R", local = TRUE)

emit <- function(df, stem, caption, digits = 1) {
  write_csv(df, file.path(DIR_TABLES, paste0(stem, ".csv")))
  cat("\n\n===== ", caption, " =====\n", sep = "")
  print(as.data.frame(df), row.names = FALSE, digits = 4)
  cat("  -> ", stem, ".csv\n", sep = "")
  invisible(df)
}

cat("\nReading pipeline results...\n")
res_all <- read_all_runs()

fam_gulf <- res_all %>%
  filter(level == "all_gulf") %>%
  freq_by_resolution(keys = c("currency", "contrast")) %>%
  add_families(keys = c("currency", "contrast"))

fam_unit <- res_all %>%
  filter(level != "all_gulf") %>%
  freq_by_resolution(keys = c("currency", "contrast", "level", "spatial_unit")) %>%
  add_families(keys = c("currency", "contrast", "level", "spatial_unit")) %>%
  left_join(
    n_units_by(filter(res_all, level != "all_gulf"),
               keys = c("currency", "contrast", "level", "spatial_unit")),
    by = c("currency", "contrast", "level", "spatial_unit")
  )

readings <- c(
  "Stable Diet"                 = "Same prey, same proportions, same breadth.",
  "Emerging Shift"              = "Early compositional signal, not yet significant.",
  "Ghost Shift"                 = "Prey identities turn over; diversity and breadth hold.",
  "Niche Compression/Expansion" = "Range of prey used changes on a stable prey list.",
  "Internal Rebalancing"        = "Proportions shift within the same prey set.",
  "Niche Restructuring"         = "Diversity and breadth both change without turnover.",
  "Partial Diet Shift"          = "Turnover accompanied by a change in diversity.",
  "Structural Shift"            = "Turnover accompanied by a change in breadth.",
  "Major Shift"                 = "Turnover with change in both diversity and breadth."
)

table1 <- STATE_LOGIC %>%
  mutate(
    family = as.character(family_of(diagnostic)),
    reading = unname(readings[diagnostic])
  ) %>%
  select(Family = family, State = diagnostic,
         Composition = composition, `Shannon H'` = H_prime, `Niche breadth Bs` = Bs,
         `Ecological reading` = reading) %>%
  arrange(factor(Family, levels = FAMILY_LEVELS),
          factor(State, levels = DIAG_LEVELS))

emit(table1, "Table1_typology",
     "Table 1. The nine diagnostic states, their defining signal combination, and the four families they collapse onto.")

table2 <- fam_gulf %>%
  filter(!is.na(contrast)) %>%
  arrange(currency, contrast) %>%
  mutate(across(all_of(FAMILY_LEVELS), ~round(.x, 1))) %>%
  transmute(Currency = CURRENCY_LAB[as.character(currency)],
            Contrast = CONTRAST_LAB_1L[as.character(contrast)],
            Stability, Substitution,
            `Functional change`, Reorganisation,
            `Turnover total` = round(Substitution + Reorganisation, 1))

emit(table2, "Table2_headline",
     "Table 2. Frequency (%) of the four families by contrast and currency, Gulf-wide level. Turnover total = Substitution + Reorganisation.")

fam_level <- bind_rows(
  fam_gulf %>% mutate(level = "all_gulf", .before = 1),
  fam_unit %>%
    group_by(currency, contrast, level) %>%
    summarise(across(all_of(FAMILY_LEVELS), ~mean(.x, na.rm = TRUE)),
              n_units_total = dplyr::n(), .groups = "drop")
) %>%
  mutate(level = factor(level, levels = SPATIAL_LEVELS_ORD))

table3 <- fam_level %>%
  filter(!is.na(contrast)) %>%
  arrange(currency, contrast, level) %>%
  mutate(across(all_of(FAMILY_LEVELS), ~round(.x, 1))) %>%
  transmute(Currency = CURRENCY_LAB[as.character(currency)],
            Contrast = contrast_label(contrast),
            `Spatial level` = LEVEL_SHORT[as.character(level)],
            `Units` = ifelse(is.na(n_units_total), 1L, n_units_total),
            Stability, Substitution,
            `Functional change`, Reorganisation)

emit(table3, "Table3_crossscale",
     paste0("Table 3. Family frequency (%) by spatial level. Gulf-wide is a single ",
            "pooled test on all trawl sets; ecoregion and stratum values are means ",
            "across units, each unit tested on its own data."))

DATA_PATH <- "data/dat_classed.rda"
LMAT_PATH <- "data/lmat_predators.csv"

if (file.exists(DATA_PATH) && file.exists(LMAT_PATH)) {
  e <- new.env(); load(DATA_PATH, envir = e)
  raw <- get(ls(e)[1], envir = e)
  if (inherits(raw, "sf")) raw <- sf::st_drop_geometry(raw)
  raw <- as.data.frame(raw)

  has_code <- "predator_species_code" %in% names(raw)
  if (!has_code) raw$predator_species_code <- NA_integer_

  lmat_ref <- read_csv(LMAT_PATH, show_col_types = FALSE,
                       locale = locale(encoding = "Latin1")) %>%
    transmute(predator_species_code = as.integer(code),
              predator_species_common_name = Predator,
              `Scientific name` = latin_name,
              `Lmat (cm)` = ifelse(tolower(lmat) %in% c("unknown", "na", ""), NA_character_, lmat),
              `Lmax (cm)` = as.character(lmax),
              Source = coalesce(`Information source`, "")) %>%
    distinct(predator_species_code, .keep_all = TRUE)
  join_key <- if (has_code) "predator_species_code" else "predator_species_common_name"
  lmat_ref <- lmat_ref %>% select(-all_of(setdiff(c("predator_species_code", "predator_species_common_name"), join_key)))

  fmt_range <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) return(NA_character_)
    paste0(format(round(min(x), 1), nsmall = 0), "-", format(round(max(x), 1), nsmall = 0))
  }

  per_period <- raw %>%
    filter(period %in% c("2004-2006", "2018-2019")) %>%
    distinct(predator_species_code, predator_species_common_name, period,
             stomach_id, somatic_length_cm) %>%
    group_by(predator_species_code, predator_species_common_name, period) %>%
    summarise(n_sto = n_distinct(stomach_id),
              range = fmt_range(somatic_length_cm), .groups = "drop")

  tableS1 <- per_period %>%
    pivot_wider(names_from = period, values_from = c(range, n_sto)) %>%
    left_join(lmat_ref, by = join_key) %>%
    transmute(Predator = predator_species_common_name,
              `Scientific name`,
              `Length range 2004-2006 (cm)` = `range_2004-2006`,
              `Length range 2018-2019 (cm)` = `range_2018-2019`,
              `Stomachs 2004-2006` = coalesce(`n_sto_2004-2006`, 0L),
              `Stomachs 2018-2019` = coalesce(`n_sto_2018-2019`, 0L),
              `Lmat (cm)`, `Lmax (cm)`, Source) %>%
    arrange(desc(`Stomachs 2004-2006` + `Stomachs 2018-2019`))

  emit(tableS1, "TableS1_predators",
       paste0("Table S1. The ", nrow(tableS1), " predators retained: length range and number ",
              "of stomachs with prey per period, length at maturity (Lmat) and maximum ",
              "length (Lmax) with their source. Predators without a documented Lmat form ",
              "a single size class."))
} else {
  message("TableS1 skipped: ", DATA_PATH, " or ", LMAT_PATH, " not found.")
}

one_res <- res_all %>%
  filter(level == "all_gulf", currency == "biomass", contrast == "PTa") %>%
  filter(x_threshold == min(x_threshold, na.rm = TRUE))

tableS2 <- one_res %>%
  transmute(Predator = species, `Size class` = size_class,
            `Stomachs 2004-2006` = n_sto_P1, `Stomachs 2018-2019` = n_sto_P2,
            `Sets 2004-2006` = n_set_P1, `Sets 2018-2019` = n_set_P2,
            Testable = testable, Reliable = reliable) %>%
  arrange(Predator, `Size class`)

if (nrow(tableS2)) {
  emit(tableS2, "TableS2_sampling",
       paste0("Table S2. Stomachs and trawl sets per predator and size class, ", CONTRAST_LAB_1L[["PTa"]], " contrast, biomass currency, finest taxonomic resolution."))
} else {
  message("TableS2 skipped: no 2004-2006 vs 2018-2019 biomass run at the Gulf-wide level.")
}

tableS3 <- res_all %>%
  filter(level == "all_gulf") %>%
  freq_by_resolution(keys = c("currency", "contrast")) %>%
  filter(!is.na(contrast)) %>%
  mutate(Family = as.character(family_of(diagnostic)),
         pct = round(pct, 1)) %>%
  mutate(contrast = contrast_label(contrast)) %>%
  select(Currency = currency, Contrast = contrast, Family,
         State = diagnostic, `Frequency (%)` = pct) %>%
  arrange(Currency, Contrast,
          factor(Family, levels = FAMILY_LEVELS),
          factor(State, levels = DIAG_LEVELS))

emit(tableS3, "TableS3_states_gulf",
     "Table S3. Frequency (%) of each of the nine diagnostic states, Gulf-wide level.")

unit_table <- function(lvl, stem, caption) {
  d <- fam_unit %>%
    filter(level == lvl, !is.na(contrast)) %>%
    arrange(currency, contrast, spatial_unit) %>%
    mutate(across(all_of(FAMILY_LEVELS), ~round(.x, 1)),
           flag = ifelse(!is.na(n_units) & n_units < 5, "*", "")) %>%
    transmute(Currency = CURRENCY_LAB[as.character(currency)],
              Contrast = contrast_label(contrast),
              Unit = paste0(area_label(spatial_unit), flag),
              `n units` = n_units,
              Stability, Substitution,
              `Functional change`, Reorganisation)
  if (!nrow(d)) { message(stem, " skipped: level '", lvl, "' not present."); return(invisible(NULL)) }
  emit(d, stem, caption)
}

unit_table("ecoregion", "TableS4_by_ecoregion",
           "Table S4. Family frequency (%) per ecoregion. * = fewer than five predator x size-class units.")

unit_table("stratum", "TableS5_by_stratum",
           "Table S5. Family frequency (%) per survey stratum. * = fewer than five predator x size-class units.")

tableS6 <- coverage_by(res_all, keys = c("level", "currency", "contrast")) %>%
  filter(!is.na(contrast)) %>%
  arrange(level, currency, contrast) %>%
  transmute(`Spatial level` = LEVEL_SHORT[as.character(level)],
            Currency = CURRENCY_LAB[as.character(currency)],
            Contrast = contrast_label(contrast),
            `Cells (all resolutions)` = n_cells,
            `Testable (%)` = pct_testable,
            `Reliable (%)` = pct_reliable,
            `Location/dispersion confounded (%)` = pct_confounded)

emit(tableS6, "TableS6_coverage",
     paste0("Table S6. Coverage and confounding by level and contrast. Percentages ",
            "in Tables 2 to S4 are computed on testable cells only."))

cat("\n\n===== Consistency checks =====\n")

chk_sum <- fam_gulf %>%
  mutate(total = rowSums(across(all_of(FAMILY_LEVELS)))) %>%
  filter(abs(total - 100) > 0.01)
if (nrow(chk_sum)) {
  warning("Family percentages do not sum to 100 in ", nrow(chk_sum),
          " Gulf-wide row(s). Inconclusive cells may not have been excluded.",
          call. = FALSE)
  print(chk_sum)
} else {
  cat("  Family percentages sum to 100 at the Gulf-wide level.\n")
}

n_lv <- dplyr::n_distinct(res_all$level)
cat("  Spatial levels present: ", n_lv, " of 3",
    if (n_lv < 3) "  <- run 04_Run_Main_Analyses.R at all three spatial levels before quoting Table 3" else "", "\n", sep = "")

miss_ct <- setdiff(CONTRAST_LEVELS, unique(as.character(res_all$contrast)))
if (length(miss_ct)) {
  cat("  Contrasts absent from the results: ", paste(miss_ct, collapse = ", "), "\n", sep = "")
}

low_cov <- tableS6 %>% filter(`Testable (%)` < 60)
if (nrow(low_cov)) {
  cat("  Level/contrast combinations below 60 % testable cells:\n")
  print(as.data.frame(low_cov), row.names = FALSE)
}

cat("\nTables written to ", normalizePath(DIR_TABLES), "\n", sep = "")

})

# Sampling and prey-grouping summaries
local({
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(readr)
  library(stringr); library(sf); library(forcats)
})

PREY_FAMILY <- "_1"
RDA_DIRS <- c(
  all_gulf  = "Sensitivity_all_gulf_1_T",
  ecoregion = "Sensitivity_ecoregion_1_T",
  stratum   = "Sensitivity_stratum_1_T"
)
source("R_helpers/Results.R", local = TRUE)

DATA_PATH   <- "data/dat_classed.rda"
STRATA_PATH <- "strata_rv_gulf.rds"
ECO_PATH    <- "data/Spatial_data/new_gulf_final.shp"
MAP_XLIM    <- c(-66.2, -60.0)
MAP_YLIM    <- c(45.5, 49.2)
CRS_MAP     <- "+proj=lcc +lat_1=46 +lat_2=49 +lat_0=47 +lon_0=-63 +datum=WGS84 +units=m"
N_FLAG      <- 5

save_app <- function(p, stem, w, h, dpi = 300) save_fig(p, stem, w, h, dpi, dir = DIR_APPEND)

emit_app <- function(df, stem, caption) {
  if ("contrast" %in% names(df)) df <- dplyr::mutate(df, contrast = contrast_label(contrast))
  write_csv(df, file.path(DIR_APPEND, paste0(stem, ".csv")))
  cat("\n\n===== ", caption, " =====\n", sep = "")
  print(as.data.frame(df), row.names = FALSE, digits = 4)
  cat("  -> ", stem, ".csv\n", sep = "")
  invisible(df)
}

cat("\nReading pipeline results...\n")
res_all <- read_all_runs()

have_raw <- file.exists(DATA_PATH)
if (have_raw) {
  e <- new.env(); load(DATA_PATH, envir = e)
  raw <- get(ls(e)[1], envir = e)
  if (inherits(raw, "sf")) raw <- sf::st_drop_geometry(raw)
  raw <- as.data.frame(raw)
  cat("Raw data: ", nrow(raw), " prey records.\n", sep = "")
} else {
  raw <- NULL
  message("Raw data not found at ", DATA_PATH,
          "; Figures 2 and 3 and the crosswalk are skipped.")
}

if (have_raw && file.exists(ECO_PATH) && all(c("longitude", "latitude") %in% names(raw))) {
  cat("\nFigure 2: trawl sets and ecoregions\n")

  eco <- st_make_valid(st_read(ECO_PATH, quiet = TRUE)) %>% mutate(Area = area_label(Area))
  sets_xy <- raw %>%
    filter(year %in% c(2004:2006, 2018:2019), !is.na(longitude), !is.na(latitude)) %>%
    distinct(year, vessel.code, set, longitude, latitude, Area) %>%
    mutate(Area = area_label(Area))
  panels <- c(sort(unique(as.character(sets_xy$year))), "All years")
  pts <- bind_rows(sets_xy %>% mutate(panel = as.character(year)),
                   sets_xy %>% mutate(panel = "All years")) %>%
    mutate(panel = factor(panel, levels = panels)) %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

  land <- NULL
  if (requireNamespace("rnaturalearth", quietly = TRUE)) {
    land <- tryCatch(rnaturalearth::ne_states(country = c("canada", "united states of america"),
                                              returnclass = "sf"),
                     error = function(e) tryCatch(
                       rnaturalearth::ne_countries(scale = "medium", returnclass = "sf"),
                       error = function(e2) NULL))
  }
  bb <- st_bbox(st_transform(
    st_as_sfc(st_bbox(c(xmin = MAP_XLIM[1], xmax = MAP_XLIM[2],
                        ymin = MAP_YLIM[1], ymax = MAP_YLIM[2]), crs = st_crs(4326))),
    CRS_MAP))

  p2 <- ggplot()
  if (!is.null(land)) p2 <- p2 + geom_sf(data = land, fill = "grey90", colour = "grey70", linewidth = 0.2)
  p2 <- p2 +
    geom_sf(data = eco, aes(fill = Area), alpha = 0.25, colour = "grey30", linewidth = 0.3) +
    geom_sf(data = pts, aes(colour = Area), size = 0.6, alpha = 0.75, show.legend = FALSE) +
    facet_wrap(~panel, ncol = 3) +
    scale_fill_brewer(palette = "Set2", name = NULL) +
    scale_colour_brewer(palette = "Set2") +
    scale_x_continuous(breaks = seq(-65, -61, by = 1)) +
    scale_y_continuous(breaks = seq(46, 49, by = 1)) +
    coord_sf(crs = CRS_MAP, xlim = bb[c("xmin", "xmax")], ylim = bb[c("ymin", "ymax")],
             expand = FALSE) +
    labs(title = "Trawl sets with analysed stomachs, by survey year and pooled",
         x = NULL, y = NULL) +
    theme_diag(base_size = 10) +
    theme(axis.text = element_text(size = 6), panel.grid = element_line(colour = "grey93"),
          panel.spacing.x = unit(1.2, "lines"), legend.position = "bottom")
  if (requireNamespace("ggspatial", quietly = TRUE))
    p2 <- p2 + ggspatial::annotation_scale(
      data = data.frame(panel = factor("All years", levels = panels)),
      location = "br", width_hint = 0.25, text_cex = 0.6)

  save_fig(p2, "Fig2_sets_ecoregions", width = 11, height = 7.6, dir = DIR_FIGURES)
}

cat("\nFigure 4: decision tree\n")
source("R_helpers/Typology.R", local = TRUE)
draw_typology_tree(out_prefix = file.path(DIR_FIGURES, "Fig4_typology_tree"))

if (have_raw) {
  cat("\nFigures S1-S2: sampling\n")

  samp <- raw %>%
    filter(year %in% c(2004, 2005, 2006, 2018, 2019)) %>%
    distinct(stomach_id, year, Area, stratum) %>%
    count(year, Area, stratum, name = "n_stomachs")

  pA <- samp %>%
    count(year, Area, wt = n_stomachs, name = "n") %>%
    ggplot(aes(factor(year), n, fill = Area)) +
    geom_col(width = 0.72, colour = "grey25", linewidth = 0.2) +
    geom_text(aes(label = n), position = position_stack(vjust = 0.5),
              size = 2.6, colour = "white", fontface = "bold") +
    scale_fill_brewer(palette = "Set2", name = NULL) +
    labs(title = "Stomachs analysed per year and ecoregion",
         x = NULL, y = "Stomachs") +
    theme_diag()

  save_app(pA, "FigS1_sampling_by_year", 8, 4.6)

  strata_sf <- load_strata(STRATA_PATH)
  if (!is.null(strata_sf)) {
    strata_sf <- st_make_valid(strata_sf)
    if (is.na(st_crs(strata_sf))) st_crs(strata_sf) <- 4326
    strata_sf$str <- as.character(strata_sf$str)

    d <- strata_sf %>%
      left_join(samp %>% mutate(str = as.character(stratum)) %>%
                  count(str, year, wt = n_stomachs, name = "n"),
                by = "str") %>%
      filter(!is.na(n))

    pB <- ggplot() +
      geom_sf(data = strata_sf, fill = "grey97", colour = "grey75", linewidth = 0.12) +
      geom_sf(data = d, aes(fill = n), colour = "grey45", linewidth = 0.12) +
      facet_wrap(~year, nrow = 1) +
      scale_fill_viridis_c(option = "rocket", direction = -1, trans = "sqrt",
                           name = "Stomachs") +
      coord_sf(xlim = MAP_XLIM, ylim = MAP_YLIM, expand = FALSE) +
      labs(title = "Stomachs analysed per stratum and survey year") +
      theme_diag(base_size = 10) +
      theme(axis.title = element_blank(), axis.text = element_text(size = 5.5),
            strip.text = element_text(face = "bold"),
            panel.grid = element_line(colour = "grey93"))

    save_app(pB, "FigS2_sampling_by_stratum", 14, 4.2, dpi = 200)
  }
}

cat("\nFigure 3A: prey categories bar plot\n")

if (have_raw) {
  prey_cols <- grep(paste0("^prey_category_\\d+", PREY_FAMILY, "$"), names(raw), value = TRUE)
  if (length(prey_cols)) {
    ncat <- tibble(
      col = prey_cols,
      x_threshold = as.numeric(str_extract(prey_cols, "\\d+")),
      n_cat = vapply(prey_cols, function(cc) dplyr::n_distinct(raw[[cc]], na.rm = TRUE),
                     numeric(1))
    ) %>% arrange(x_threshold)

    p3a <- ggplot(ncat, aes(x = x_threshold, y = n_cat)) +
      geom_col(fill = "#2F4A5A", width = 7, alpha = 0.85) +
      scale_x_continuous(
        breaks = seq(0, 1000, by = 50),
        limits = c(-10, 1020)
      ) +
      scale_y_continuous(
        breaks = seq(0, 200, by = 10),
        limits = c(0, 190)
      ) +
      labs(title = "A. Prey categories retained by aggregation threshold",
           subtitle = "Each record is kept at its finest level, or walked up the taxonomic tree until the node reaches x stomachs",
           x = "Threshold x (stomachs)", y = "Distinct prey categories") +
      theme_diag() +
      theme(
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        panel.grid.minor = element_blank()
      )

    save_fig(p3a, "Fig3a_prey_categories", 9, 5, dir = DIR_FIGURES)
    emit_app(ncat %>% select(-col), "TableA0_prey_categories",
             "Prey categories retained at each aggregation threshold.")
  } else {
    message("No prey_category_*", PREY_FAMILY, " columns found; Figure 3A skipped.")
  }
}

sweep_fam <- res_all %>%
  filter(level == "all_gulf", diagnostic %in% DIAG_LEVELS) %>%
  count(currency, contrast, x_threshold, diagnostic, name = "n") %>%
  group_by(currency, contrast, x_threshold) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(family = family_of(diagnostic)) %>%
  filter(!is.na(family), !is.na(contrast)) %>%
  group_by(currency, contrast, x_threshold, family) %>%
  summarise(pct = sum(pct), .groups = "drop")

p3b <- ggplot(sweep_fam, aes(x_threshold, pct, colour = family)) +
  geom_line(linewidth = 0.6) +
  facet_grid(currency ~ contrast,
             labeller = labeller(currency = CURRENCY_LAB,
                                 contrast = CONTRAST_LAB_1L)) +
  scale_colour_manual(values = FAMILY_PAL, name = NULL) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(title = "B. Sensitivity of the diagnostic composition to prey taxonomic resolution",
       subtitle = "Family frequencies across the full sweep; the reported values are the means of these curves",
       x = "Threshold x (stomachs)", y = "Frequency (%)") +
  theme_diag() +
  theme(strip.text = element_text(size = 8))

save_fig(p3b, "Fig3b_sweep_families", 13, 5.6, dir = DIR_FIGURES)

cat("\nFigure 3b: taxonomic resolution sweep - Ecoregions\n")

sweep_eco <- res_all %>%
  filter(level == "ecoregion", diagnostic %in% DIAG_LEVELS, contrast != "2006 vs 2018") %>%
  count(spatial_unit, currency, contrast, x_threshold, diagnostic, name = "n") %>%
  group_by(spatial_unit, currency, contrast, x_threshold) %>%
  mutate(pct_unit = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(family = family_of(diagnostic)) %>%
  filter(!is.na(family)) %>%
  group_by(spatial_unit, currency, contrast, x_threshold, family) %>%
  summarise(pct_unit = sum(pct_unit), .groups = "drop") %>%
  group_by(currency, contrast, x_threshold, family) %>%
  summarise(pct = mean(pct_unit), .groups = "drop")

p3b_eco <- ggplot(sweep_eco, aes(x_threshold, pct, colour = family)) +
  geom_line(linewidth = 0.6) +
  facet_grid(currency ~ contrast,
             labeller = labeller(currency = CURRENCY_LAB, contrast = CONTRAST_LAB_1L)) +
  scale_colour_manual(values = FAMILY_PAL, name = NULL) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(title = "Sensitivity of diagnostic composition to taxonomic resolution — Ecoregions",
       subtitle = "Family frequencies across the full sweep (mean across ecoregions)",
       x = "Threshold x (stomachs)", y = "Frequency (%)") +
  theme_diag() +
  theme(strip.text = element_text(size = 8), legend.position = "top")

save_fig(p3b_eco, "Fig3b_sweep_families_ecoregion", 13, 5.6, dir = DIR_FIGURES)

cat("\nFigure 3b: taxonomic resolution sweep - Survey Strata\n")

sweep_strata <- res_all %>%
  filter(level == "stratum", diagnostic %in% DIAG_LEVELS, contrast != "2006 vs 2018") %>%
  count(spatial_unit, currency, contrast, x_threshold, diagnostic, name = "n") %>%
  group_by(spatial_unit, currency, contrast, x_threshold) %>%
  mutate(pct_unit = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(family = family_of(diagnostic)) %>%
  filter(!is.na(family)) %>%
  group_by(spatial_unit, currency, contrast, x_threshold, family) %>%
  summarise(pct_unit = sum(pct_unit), .groups = "drop") %>%
  group_by(currency, contrast, x_threshold, family) %>%
  summarise(pct = mean(pct_unit), .groups = "drop")

p3b_strata <- ggplot(sweep_strata, aes(x_threshold, pct, colour = family)) +
  geom_line(linewidth = 0.6) +
  facet_grid(currency ~ contrast,
             labeller = labeller(currency = CURRENCY_LAB, contrast = CONTRAST_LAB_1L)) +
  scale_colour_manual(values = FAMILY_PAL, name = NULL) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(title = "Sensitivity of diagnostic composition to taxonomic resolution — Survey Strata",
       subtitle = "Family frequencies across the full sweep (mean across strata)",
       x = "Threshold x (stomachs)", y = "Frequency (%)") +
  theme_diag() +
  theme(strip.text = element_text(size = 8), legend.position = "top")

save_fig(p3b_strata, "Fig3b_sweep_families_stratum", 13, 5.6, dir = DIR_FIGURES)

sweep_spread <- sweep_fam %>%
  group_by(currency, contrast, family) %>%
  summarise(mean_pct = round(mean(pct), 1),
            sd_pct   = round(sd(pct), 1),
            min_pct  = round(min(pct), 1),
            max_pct  = round(max(pct), 1),
            n_resolutions = dplyr::n(), .groups = "drop")

emit_app(sweep_spread, "TableA5_resolution_spread",
         "Family frequency across the taxonomic-resolution sweep: mean, spread and range.")

cat("\nAudits\n")

tableA1 <- res_all %>%
  mutate(k = pmin(n_set_P1, n_set_P2)) %>%
  group_by(level, currency) %>%
  summarise(
    n_cells   = dplyr::n(),
    k_min     = min(k, na.rm = TRUE),
    k_q25     = quantile(k, 0.25, na.rm = TRUE),
    k_median  = median(k, na.rm = TRUE),
    k_q75     = quantile(k, 0.75, na.rm = TRUE),
    pct_k_lt4 = round(100 * mean(k < 4, na.rm = TRUE), 1),
    pct_k_lt6 = round(100 * mean(k < 6, na.rm = TRUE), 1),
    stomachs_per_set = round(median((n_sto_P1 + n_sto_P2) /
                                      pmax(n_set_P1 + n_set_P2, 1), na.rm = TRUE), 1),
    .groups = "drop"
  )

emit_app(tableA1, "TableA1_set_depth",
         "Trawl sets per cell by spatial level. k = min(sets in 2004-2006, sets in 2018-2019).")

tableA2 <- res_all %>%
  filter(testable) %>%
  group_by(level, currency, contrast) %>%
  summarise(
    n            = dplyr::n(),
    k_median     = median(pmin(n_set_P1, n_set_P2), na.rm = TRUE),
    BC_median    = round(median(BC, na.rm = TRUE), 3),
    R2_median    = round(median(R2_comp, na.rm = TRUE), 3),
    pct_H_sig    = round(100 * mean(p_H    < 0.05, na.rm = TRUE), 1),
    pct_Bs_sig   = round(100 * mean(p_Bs   < 0.05, na.rm = TRUE), 1),
    pct_comp_sig = round(100 * mean(p_comp < 0.05, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  filter(!is.na(contrast)) %>%
  arrange(currency, contrast, level)

emit_app(tableA2, "TableA2_power_vs_effect",
         "Effect sizes and detection rates by level. BC and R2 do not depend on significance; the pct_*_sig columns do.")

tableA4 <- res_all %>%
  group_by(level, currency, contrast) %>%
  summarise(
    n_cells          = dplyr::n(),
    n_inconclusive   = sum(!testable),
    n_no_H           = sum(is.na(p_H)),
    n_no_Bs          = sum(is.na(p_Bs)),
    n_no_composition = sum(is.na(p_comp)),
    n_below_min_sets = sum(n_set_P1 < 3 | n_set_P2 < 3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(contrast), n_inconclusive > 0) %>%
  arrange(desc(n_inconclusive))

emit_app(tableA4, "TableA4_untestable",
         "Cells that could not be classified, and which signal was missing.")

if (have_raw && all(c("Area", "stratum") %in% names(raw))) {
  cat("\nPartition crosswalk\n")

  xw <- raw %>%
    mutate(stratum = as.character(stratum)) %>%
    count(Area, stratum, name = "n_rec") %>%
    group_by(stratum) %>%
    mutate(n_stratum = sum(n_rec),
           share = round(100 * n_rec / n_stratum, 1),
           n_ecoregions = dplyr::n()) %>%
    ungroup() %>%
    arrange(desc(n_ecoregions), stratum, desc(n_rec))

  tableA3 <- xw %>%
    filter(n_ecoregions > 1) %>%
    select(Stratum = stratum, Ecoregion = Area,
           `Records` = n_rec, `Share of stratum (%)` = share)

  emit_app(tableA3, "TableA3_partition_crosswalk",
           "Strata split across more than one ecoregion. Ecoregion is assigned per stomach by spatial join; strata are survey polygons.")

  summ <- xw %>%
    distinct(stratum, n_stratum, n_ecoregions) %>%
    summarise(
      n_strata     = dplyr::n(),
      n_split      = sum(n_ecoregions > 1),
      pct_records_in_split = round(
        100 * sum(n_stratum[n_ecoregions > 1]) / sum(n_stratum), 1)
    )
  cat("\nCrosswalk summary:\n"); print(as.data.frame(summ), row.names = FALSE)
}

cat("\nFigure S7: low-sample sensitivity\n")

S5_BAND <- "q95"

res_s7 <- res_all %>%
  mutate(spatial_unit = if_else(level == "all_gulf", "Gulf-wide", spatial_unit))

unit_keys <- c("currency", "contrast", "level", "spatial_unit")

cells_per_unit <- n_units_by(res_s7, keys = unit_keys) %>%
  rename(n_cells = n_units)

fam_ux <- family_by_resolution(res_s7, keys = unit_keys) %>%
  left_join(cells_per_unit, by = unit_keys) %>%
  filter(!is.na(contrast))

sub_lab <- paste0("n cells >= ", N_FLAG)

comp_by_x <- bind_rows(
  fam_ux %>% mutate(subset = "All units"),
  fam_ux %>%
    filter(level != "all_gulf", !is.na(n_cells), n_cells >= N_FLAG) %>%
    mutate(subset = sub_lab)
) %>%
  group_by(level, currency, contrast, subset, family, x_threshold) %>%
  summarise(pct = mean(pct), n_units_kept = n_distinct(spatial_unit),
            .groups = "drop")

sens <- comp_by_x %>%
  group_by(level, currency, contrast, subset, family) %>%
  summarise(pct_mean = mean(pct), pct_sd = sd(pct),
            q_lo = unname(quantile(pct, 0.025)), q_hi = unname(quantile(pct, 0.975)),
            n_res = n_distinct(x_threshold), n_units_kept = min(n_units_kept),
            .groups = "drop") %>%
  mutate(lo = if (S5_BAND == "sd") pmax(pct_mean - pct_sd, 0) else q_lo,
         hi = if (S5_BAND == "sd") pmin(pct_mean + pct_sd, 100) else q_hi,
         level    = factor(level, levels = SPATIAL_LEVELS_ORD),
         subset   = factor(subset, levels = c("All units", sub_lab)),
         family   = factor(family, levels = FAMILY_LEVELS),
         contrast = factor(contrast, levels = CONTRAST_LEVELS),
         x = as.integer(contrast) + if_else(subset == "All units", -0.07, 0.07))

band_lab <- if (S5_BAND == "sd") "mean +/- SD across taxonomic resolutions" else
  "95 % interval (2.5-97.5 %) across taxonomic resolutions"

pS7 <- ggplot(sens, aes(x = x, y = pct_mean, colour = family,
                        shape = subset, group = interaction(family, subset))) +
  geom_errorbar(aes(ymin = lo, ymax = hi, linetype = subset),
                width = 0.1, linewidth = 0.35, alpha = 0.85) +
  geom_line(aes(linetype = subset), linewidth = 0.6) +
  geom_point(size = 2.2, fill = "white") +
  scale_x_continuous(breaks = seq_along(CONTRAST_LEVELS), labels = unname(CONTRAST_LAB),
                     limits = c(0.5, length(CONTRAST_LEVELS) + 0.5), expand = c(0, 0)) +
  facet_grid(currency ~ level,
             labeller = labeller(currency = CURRENCY_LAB, level = LEVEL_SHORT)) +
  scale_colour_manual(values = FAMILY_PAL, name = NULL) +
  scale_shape_manual(values = c(16, 21), name = NULL, drop = FALSE) +
  scale_linetype_manual(values = c("solid", "dashed"), name = NULL, drop = FALSE) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(title = "Sensitivity of the across-unit composites to sparsely sampled units",
       subtitle = paste0("Solid = every unit; dashed = units built on at least ", N_FLAG,
                         " predator x size-class cells. Bars: ", band_lab,
                         ".\nGulf-wide is a single unit and is shown as the reference."),
       x = NULL, y = "Frequency (%)") +
  theme_diag() +
  theme(legend.position = "top")

save_app(pS7, "FigS7_lowN_sensitivity", 13, 6.5)

emit_app(
  sens %>%
    select(level, currency, contrast, subset, family, pct_mean, lo, hi,
           n_res, n_units_kept) %>%
    mutate(across(c(pct_mean, lo, hi), ~round(.x, 1))) %>%
    pivot_wider(names_from = family, values_from = c(pct_mean, lo, hi),
                names_glue = "{family}_{.value}"),
  "TableA6_lowN_sensitivity",
  paste0("Across-unit composites with and without sparsely sampled units; ",
         band_lab, "."))

if (!is.null(raw) && "tax_level_lowest" %in% names(raw)) {
  cat("\nTable A7: identification depth by period\n")

  e2 <- new.env(); sys.source("R_helpers/Taxonomy.R", envir = e2)
  rank_order <- get("ranknfile", envir = e2)

  parent_rank <- function(r) {
    r <- tolower(as.character(r))
    out <- rep(NA_character_, length(r))
    for (p in c("phylum", "class", "order", "family", "genus", "species")) {
      hit <- grepl(paste0("^(super|sub|infra|parv|mega|giga|magn|grand)?", p, "$"), r)
      out[hit & is.na(out)] <- p
    }
    out
  }

  id <- raw %>%
    mutate(period = case_when(year %in% 2004:2006 ~ "2004-2006",
                              year %in% 2018:2019 ~ "2018-2019",
                              TRUE ~ NA_character_)) %>%
    filter(!is.na(period), !is.na(tax_level_lowest)) %>%
    transmute(period,
              rank_raw = tolower(as.character(tax_level_lowest)),
              rank     = coalesce(parent_rank(tax_level_lowest), "life-stage prey category"),
              depth    = match(rank_raw, rank_order))

  rank_lv <- c("phylum", "class", "order", "family", "genus", "species",
               "life-stage prey category")
  i_of <- function(p) match(p, rank_order)

  at_rank <- id %>%
    count(period, rank, name = "n") %>%
    group_by(period) %>%
    mutate(pct_at = round(100 * n / sum(n), 1)) %>%
    ungroup() %>%
    select(-n)

  reach <- id %>%
    filter(!is.na(depth)) %>%
    group_by(period) %>%
    summarise(n_items = n(),
              class   = round(100 * mean(depth >= i_of("class"),   na.rm = TRUE), 1),
              order   = round(100 * mean(depth >= i_of("order"),   na.rm = TRUE), 1),
              family  = round(100 * mean(depth >= i_of("family"),  na.rm = TRUE), 1),
              genus   = round(100 * mean(depth >= i_of("genus"),   na.rm = TRUE), 1),
              species = round(100 * mean(depth >= i_of("species"), na.rm = TRUE), 1),
              .groups = "drop") %>%
    pivot_longer(c(class, order, family, genus, species),
                 names_to = "rank", values_to = "pct_at_or_finer")

  tabA7 <- expand_grid(rank = rank_lv, period = unique(id$period)) %>%
    left_join(at_rank, by = c("period", "rank")) %>%
    left_join(reach %>% select(period, rank, pct_at_or_finer), by = c("period", "rank")) %>%
    mutate(pct_at = coalesce(pct_at, 0)) %>%
    pivot_wider(names_from = period, values_from = c(pct_at, pct_at_or_finer),
                names_glue = "{.value}_{period}") %>%
    mutate(rank = factor(rank, levels = rank_lv)) %>%
    arrange(rank)

  n_items <- id %>% count(period, name = "n_items")
  cat("Prey items: ", paste(n_items$period, n_items$n_items, collapse = " | "), "\n", sep = "")

  emit_app(tabA7, "TableA7_identification_depth",
           paste0("Depth of prey identification by period (% of prey items whose lowest ",
                  "recorded rank is each rank, and % identified to that rank or finer); n = ",
                  paste(n_items$n_items, collapse = " / "), " items"))
} else {
  message("tax_level_lowest not in dat_classed: Table A7 skipped.")
}

cat("\nManifest\n")

manifest <- file.path(DIR_APPEND, "manifest.txt")
con <- file(manifest, "w")

wl <- function(...) writeLines(paste0(...), con)

wl("TROPHIC TRANSITION ANALYSIS - RUN MANIFEST")
wl("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
wl(strrep("-", 78))
wl("")
wl("SPATIAL LEVELS PRESENT")
for (lv in levels(res_all$level)) {
  n <- sum(res_all$level == lv, na.rm = TRUE)
  if (n > 0) wl("  ", lv, ": ", n, " cell-rows across all resolutions and contrasts")
}
wl("")
wl("CONTRASTS PRESENT")
for (ct in CONTRAST_LEVELS) {
  n <- sum(as.character(res_all$contrast) == ct, na.rm = TRUE)
  wl("  ", CONTRAST_LAB_1L[[ct]], ": ", n, " cell-rows")
}
wl("")
wl("TAXONOMIC RESOLUTIONS")
wl("  ", dplyr::n_distinct(res_all$x_threshold), " thresholds, from ",
   min(res_all$x_threshold, na.rm = TRUE), " to ",
   max(res_all$x_threshold, na.rm = TRUE))
wl("")
wl("TEST SETTINGS READ BACK FROM THE RESULTS")
wl("  Niche breadth test: ", paste(unique(res_all$Bs_test), collapse = ", "))
wl("  Prey-grouping family: ",
   if ("prey_family" %in% names(res_all)) paste(unique(res_all$prey_family), collapse = ", ") else PREY_FAMILY)
wl("  Result folders: ", paste(RDA_DIRS, collapse = ", "))
wl("")
wl("OUTPUT FILES")
for (d in c(DIR_FIGURES, DIR_TABLES, DIR_APPEND)) {
  if (!dir.exists(d)) next
  fs <- list.files(d, full.names = TRUE)
  if (!length(fs)) next
  wl("  [", d, "]")
  info <- file.info(fs)
  for (i in seq_along(fs)) {
    wl(sprintf("    %-52s %8.1f KB  %s", basename(fs[i]),
               info$size[i] / 1024, format(info$mtime[i], "%Y-%m-%d %H:%M")))
  }
}
wl("")
wl("SESSION")
writeLines(capture.output(sessionInfo()), con)
close(con)

cat("  -> manifest.txt\n")
cat("\nAppendices written to ", normalizePath(DIR_APPEND), "\n", sep = "")

})
