# Manuscript figures


# Figures 1–7
local({
try(suppressWarnings(Sys.setlocale("LC_ALL", "C.UTF-8")), silent = TRUE)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr)
  library(ggplot2); library(sf); library(forcats)
})

PREY_FAMILY <- "_1"
RDA_DIRS <- c(
  all_gulf = "Sensitivity_all_gulf_1_T",
  ecoregion = "Sensitivity_ecoregion_1_T",
  stratum = "Sensitivity_stratum_1_T"
)
source("R_helpers/Results.R", local = TRUE)
sf_use_s2(FALSE)

OUT <- "Output_Paper_Figures"
dir.create(file.path(OUT, "data"), recursive = TRUE, showWarnings = FALSE)

FAM   <- c("Stability", "Substitution", "Functional change", "Reorganization")
PAL   <- setNames(unname(FAMILY_PAL), FAM)
SHAPE <- setNames(unname(FAMILY_SHAPE), FAM)
SHAPE_FILL <- setNames(c(21, 22, 23, 24), FAM)
X_OFF <- 0.12
CT    <- CONTRAST_LEVELS
CL1   <- CONTRAST_LAB_1L
CL2   <- CONTRAST_LAB
LV    <- c("all_gulf", "ecoregion", "stratum")
LL    <- c(all_gulf = "sGSL-wide", ecoregion = "Ecoregion", stratum = "Stratum")
ML    <- c(biomass = "Biomass", occurrence = "Occurrence")
TYPE_COL <- c(within = "#4C8CB5", between = "#C8553D")
TYPE_LAB <- c(within = "Within periods", between = "Between periods")
ECO_FULL <- c(Central = "Central sGSL", "Chaleur Bay" = "Chaleur Bay",
              "Magdalen Shallows" = "Magdalen Shallows",
              Northumberland = "Northumberland Strait")
MAP_CT   <- c("P1a", "P1b", "P2", "PTa")
CT_SENS  <- c("P1a", "P1b", "P2", "PTa")
N_FLAG   <- 5

clean_fam <- function(d) {
  names(d) <- sub("Reorganisation", "Reorganization", names(d), fixed = TRUE)
  for (v in intersect(c("family", "modal_family"), names(d)))
    d[[v]] <- sub("Reorganisation", "Reorganization", d[[v]], fixed = TRUE)
  d
}
# Robustness summaries
robustness_files <- c(
  "M0_family_frequencies_reference.csv",
  "M3_bootstrap_CI_family_pct.csv",
  "M4_alpha_trend_sensitivity.csv",
  "M5_cell_classification_stability.csv"
)
robustness_tags <- c("", "_T")
robustness_complete <- vapply(robustness_tags, function(tag) {
  all(vapply(LV, function(level) {
    root <- sprintf("Robustness_%s_1%s", level, tag)
    all(file.exists(file.path(root, robustness_files)))
  }, logical(1)))
}, logical(1))

if (!any(robustness_complete)) {
  missing <- unlist(lapply(robustness_tags, function(tag) {
    unlist(lapply(LV, function(level) {
      root <- sprintf("Robustness_%s_1%s", level, tag)
      paths <- file.path(root, robustness_files)
      paths[!file.exists(paths)]
    }))
  }))
  stop("Incomplete robustness summaries. Missing files:\n  ",
       paste(unique(missing), collapse = "\n  "), call. = FALSE)
}
robustness_tag <- robustness_tags[which(robustness_complete)[1L]]
message("Robustness summaries: Robustness_<level>_1", robustness_tag)

rr <- function(level, file) {
  path <- file.path(sprintf("Robustness_%s_1%s", level, robustness_tag), file)
  read_csv(path, show_col_types = FALSE) %>%
    clean_fam() %>% mutate(level = level)
}

theme_p <- function(base = 10) {
  theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom", legend.title = element_blank(),
          axis.text = element_text(colour = "#222222"),
          plot.margin = margin(5.5, 12, 5.5, 5.5))
}
theme_set(theme_p())

save_p <- function(p, id, w, h, data = NULL) {
  ggsave(file.path(OUT, paste0(id, ".png")), p, width = w, height = h,
         dpi = 300, bg = "white", limitsize = FALSE)
  if (!is.null(data)) write_csv(data, file.path(OUT, "data", paste0(id, ".csv")))
  cat("  ", id, "\n", sep = "")
}

# Figure inputs
required_files <- c(
  "data/dat_classed.rda",
  "data/Spatial_data/new_gulf_final.shp",
  "Output_Figures/data_Fig6_9_unit_families.csv",
  "Output_Appendices/TableA0_prey_categories.csv",
  "Output_Appendices/TableA7_identification_depth.csv",
  "NullCalibration_all_gulf_1/null_fpr_by_sample_size.csv",
  file.path("Robustness_all_gulf_1_BAL",
            c("balanced_family_per_replicate.csv", "balanced_vs_unbalanced.csv")),
  file.path("Robustness_ecoregion_1_BAL",
            c("balanced_family_per_replicate.csv", "balanced_vs_unbalanced.csv"))
)
required_dirs <- c(RDA_DIRS, "Sensitivity_all_gulf_2", "Sensitivity_all_gulf_PP")
missing_files <- required_files[!file.exists(required_files)]
missing_dirs <- required_dirs[!dir.exists(required_dirs)]
if (length(missing_files) || length(missing_dirs)) {
  missing <- c(missing_files, paste0(missing_dirs, "/"))
  stop("Missing figure inputs:\n  ", paste(missing, collapse = "\n  "), call. = FALSE)
}

cat("Reading saved results...\n")

ref <- map_dfr(LV, rr, file = "M0_family_frequencies_reference.csv") %>%
  mutate(family = factor(family, FAM), contrast = factor(contrast, CT))

cells <- map_dfr(LV, rr, file = "M5_cell_classification_stability.csv")
n_class <- cells %>%
  mutate(spatial_unit = sub(".*\\|", "", cell)) %>%
  count(level, contrast, mode, spatial_unit, name = "n_classifiable")

units <- read_csv("Output_Figures/data_Fig6_9_unit_families.csv", show_col_types = FALSE) %>%
  clean_fam() %>%
  mutate(contrast = names(CL1)[match(contrast, CL1)], mode = currency) %>%
  left_join(n_class, by = c("level", "contrast", "mode", "spatial_unit")) %>%
  mutate(n_classifiable = replace_na(n_classifiable, 0L))

eco_order <- names(ECO_FULL)
str_order <- as.character(sort(as.integer(unique(units$spatial_unit[units$level == "stratum"]))))
unit_order <- c(eco_order, str_order)

res_all <- read_all_runs()

regime_plot <- function(d, cts = CT, labels = TRUE, base = 10, pt = 2.6) {
  nb <- sum(cts %in% c("PTb", "PTa"))
  nw <- length(cts) - nb
  d <- d %>% filter(contrast %in% cts) %>%
    mutate(x = match(as.character(contrast), cts) + ifelse(mode == "biomass", -X_OFF, X_OFF),
           family = factor(family, FAM),
           fill_col = ifelse(mode == "biomass", unname(PAL[as.character(family)]), "white"))
  g <- ggplot(d, aes(x, pct_mean, colour = family, shape = family)) +
    annotate("rect", xmin = 0.5, xmax = nw + 0.5, ymin = -Inf, ymax = Inf, fill = "#D7E6EC", alpha = 0.55) +
    annotate("rect", xmin = nw + 0.5, xmax = nw + nb + 0.5, ymin = -Inf, ymax = Inf, fill = "#F6E0D6", alpha = 0.65) +
    geom_vline(xintercept = nw + 0.5, linetype = "dashed", colour = "grey55", linewidth = 0.3)
  if (labels) g <- g +
    annotate("text", x = (nw + 1) / 2, y = Inf, label = "WITHIN PERIOD", colour = "#2F4A5A",
             fontface = 2, size = 3.6, vjust = 1.8) +
    annotate("text", x = nw + (nb + 1) / 2, y = Inf,
             label = if (nb > 1) "BETWEEN PERIODS" else "BETWEEN\nPERIODS", colour = "#9B2226",
             fontface = 2, size = 3.6, vjust = 1.3, lineheight = 0.9)
  g +
    geom_errorbar(aes(ymin = pmax(pct_mean - sd_across_x, 0), ymax = pmin(pct_mean + sd_across_x, 100)),
                  width = 0.08, linewidth = 0.4, alpha = 0.7, show.legend = FALSE) +
    geom_point(aes(fill = fill_col, alpha = mode), size = pt, stroke = 0.9) +
    scale_x_continuous(breaks = seq_along(cts), labels = unname(CL2[cts]), expand = expansion(add = 0.25)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_colour_manual(values = PAL) + scale_shape_manual(values = SHAPE_FILL) +
    scale_fill_identity() +
    scale_alpha_manual(values = c(biomass = 1, occurrence = 1), labels = ML) +
    guides(colour = guide_legend(order = 1, override.aes = list(fill = unname(PAL))),
           shape  = guide_legend(order = 1, override.aes = list(fill = unname(PAL))),
           alpha  = guide_legend(order = 2, override.aes = list(shape = 21, colour = "grey25", size = 2.8,
                                                                fill = c("grey25", "white")))) +
    labs(x = NULL, y = "Family frequency (%)") +
    theme_p(base) + theme(legend.position = "top", legend.box = "horizontal")
}

units_plot <- function(ct) {
  d <- crossing(mode = names(ML), spatial_unit = unit_order) %>%
    left_join(filter(units, contrast == ct) %>% select(mode, spatial_unit, all_of(FAM), n_classifiable),
              by = c("mode", "spatial_unit")) %>%
    mutate(n_classifiable = replace_na(n_classifiable, 0L),
           y = length(unit_order) + 2 - match(spatial_unit, unit_order) -
             if_else(spatial_unit %in% str_order, 1.3, 0),
           Currency = factor(ML[mode], levels = ML))
  yl <- d %>% filter(mode == "biomass") %>%
    mutate(name = coalesce(unname(ECO_FULL[spatial_unit]), spatial_unit),
           label = paste0(name, " (", n_classifiable, ")",
                          if_else(n_classifiable > 0 & n_classifiable < N_FLAG, " *", "")))
  long <- d %>% pivot_longer(all_of(FAM), names_to = "family", values_to = "pct") %>%
    mutate(family = factor(family, FAM)) %>%
    group_by(mode, spatial_unit) %>% arrange(family, .by_group = TRUE) %>%
    mutate(xmax = cumsum(replace_na(pct, 0)), xmin = xmax - replace_na(pct, 0)) %>% ungroup()
  miss <- filter(d, n_classifiable == 0)
  ggplot() +
    geom_rect(data = miss, aes(xmin = 0, xmax = 100, ymin = y - .36, ymax = y + .36),
              fill = "#EFEFEF", colour = "#CCCCCC") +
    geom_text(data = miss, aes(x = 50, y = y, label = "No classifiable cell"),
              colour = "#777777", size = 2.3) +
    geom_rect(data = filter(long, n_classifiable > 0),
              aes(xmin = xmin, xmax = xmax, ymin = y - .36, ymax = y + .36, fill = family),
              colour = "white", linewidth = .12) +
    geom_hline(yintercept = length(unit_order) - 2.65, colour = "#AAAAAA", linewidth = .3) +
    facet_wrap(~Currency, nrow = 1) +
    scale_fill_manual(values = PAL, drop = FALSE) +
    scale_y_continuous(breaks = yl$y, labels = yl$label, expand = expansion(add = .7)) +
    scale_x_continuous(breaks = seq(0, 100, 25), expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 100)) +
    labs(x = "Classifiable cells (%)", y = NULL) +
    theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 8.2))
}

sweep_data <- function(lvl) {
  keys <- if (lvl == "all_gulf") c("currency", "contrast") else c("currency", "contrast", "spatial_unit")
  res_all %>% filter(level == lvl) %>%
    family_by_resolution(keys = keys) %>%
    group_by(currency, contrast, x_threshold, family) %>%
    summarise(pct = mean(pct), .groups = "drop") %>%
    filter(!is.na(contrast), lvl == "all_gulf" | contrast != "PTb") %>%
    mutate(family = factor(sub("Reorganisation", "Reorganization", family), FAM),
           Contrast = factor(CL1[as.character(contrast)], levels = unname(CL1)),
           Currency = factor(ML[currency], levels = ML))
}
sweep_plot <- function(d) {
  ggplot(d, aes(x_threshold, pct, colour = family)) +
    geom_line(linewidth = 0.5) +
    facet_grid(Currency ~ Contrast) +
    scale_colour_manual(values = PAL) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
    scale_x_continuous(breaks = c(0, 250, 500, 750, 1000)) +
    labs(x = "Aggregation threshold (stomachs)", y = "Family frequency (%)") +
    theme_p(9) + theme(axis.text.x = element_text(size = 6.5))
}

load_strata_geo <- function() {
  shp <- system.file("extdata/shapefiles/survey.stratum.polygons.shp", package = "gulf.spatial")
  if (!nzchar(shp)) shp <- "geo/survey.stratum.polygons.shp"
  if (!file.exists(shp)) return(NULL)
  s <- st_read(shp, quiet = TRUE)
  s <- s[s$survey == "rv" & s$region == "gulf" & s$type == "polygon", ]
  s$spatial_unit <- as.character(s$label)
  st_make_valid(s)
}
load_coast <- function() {
  f <- system.file("extdata/gulf.coast.intermediate.csv", package = "gulf.spatial")
  if (!nzchar(f)) f <- "geo/gulf.coast.intermediate.csv"
  if (!file.exists(f)) return(NULL)
  read_csv(f, show_col_types = FALSE) %>% filter(!is.na(x), !is.na(y))
}
strata_geo <- load_strata_geo()
coast <- load_coast()

map_plot <- function(cur, mask_low_n = FALSE) {
  if (is.null(strata_geo)) return(NULL)
  lab_xy <- strata_geo %>% st_drop_geometry() %>%
    transmute(spatial_unit, X = label_x, Y = label_y)
  md <- units %>% filter(level == "stratum", contrast %in% MAP_CT, mode == cur) %>%
    pivot_longer(all_of(FAM), names_to = "family", values_to = "pct") %>%
    mutate(family = factor(family, FAM),
           contrast = factor(CL1[contrast], levels = CL1[MAP_CT]),
           low_n = n_classifiable < N_FLAG) %>%
    filter(n_classifiable > 0)
  d <- strata_geo %>% select(spatial_unit) %>% inner_join(md, by = "spatial_unit")
  d_col  <- if (mask_low_n) filter(d, !low_n) else d
  d_grey <- if (mask_low_n) filter(d, low_n) else d[0, ]
  lab <- d_col %>% st_drop_geometry() %>% left_join(lab_xy, by = "spatial_unit") %>%
    mutate(lab = paste0(sprintf("%.0f", pct), ifelse(low_n & !mask_low_n, "*", "")),
           txt = ifelse(pct < 50, "grey95", "grey10"))
  g <- ggplot()
  if (!is.null(coast)) g <- g + geom_polygon(data = coast, aes(x, y, group = pid), fill = "grey88")
  g <- g +
    geom_sf(data = strata_geo, fill = "grey97", colour = "grey70", linewidth = 0.15) +
    geom_sf(data = d_col, aes(fill = pct), colour = "grey45", linewidth = 0.15)
  if (nrow(d_grey)) g <- g + geom_sf(data = d_grey, fill = "grey78", colour = "grey45", linewidth = 0.15)
  g +
    geom_text(data = lab, aes(X, Y, label = lab, colour = txt), size = 2.1, fontface = "bold") +
    scale_colour_identity() +
    facet_grid(family ~ contrast) +
    scale_fill_viridis_c(limits = c(0, 100), breaks = seq(0, 100, 25),
                         name = "Within-stratum frequency (%)") +
    scale_x_continuous(breaks = seq(-65, -61, 1)) +
    scale_y_continuous(breaks = seq(46, 49, 1)) +
    coord_sf(xlim = c(-66.2, -60.0), ylim = c(45.5, 49.2), expand = FALSE) +
    theme_p(10) +
    theme(panel.grid = element_line(colour = "grey92"), axis.title = element_blank(),
          axis.text = element_text(size = 6), strip.text.y = element_text(angle = 0),
          legend.position = "top", legend.title = element_text(size = 9),
          legend.key.width = unit(1.2, "cm"))
}

cat("\nManuscript figures\n")

CRS_MAP     <- "+proj=lcc +lat_1=46 +lat_2=49 +lat_0=47 +lon_0=-63 +datum=WGS84 +units=m"
STRATA_KEEP <- as.character(c(401:403, 415:439))
STRATA_COL  <- "#1F4E99"
ECO_COL     <- c("Central sGSL" = "#66C2A5", "Chaleur Bay" = "#FC8D62",
                 "Magdalen Shallows" = "#808080", "Northumberland Strait" = "#E78AC3")
LAB_NUDGE   <- tibble(spatial_unit = c("438", "439"), dx = c(-0.15, 0.12), dy = c(0.08, -0.08))

raw_env <- new.env(); load("data/dat_classed.rda", envir = raw_env)
raw <- get(ls(raw_env)[1], envir = raw_env)
if (inherits(raw, "sf")) raw <- st_drop_geometry(raw)
raw <- as.data.frame(raw) %>%
  filter(year %in% c(2004:2006, 2018:2019)) %>%
  mutate(Area = unname(ECO_FULL[area_label(Area)]),
         period = if_else(year <= 2006, "2004-2006", "2018-2019"))
eco <- st_make_valid(st_read("data/Spatial_data/new_gulf_final.shp", quiet = TRUE)) %>%
  mutate(Area = unname(ECO_FULL[area_label(Area)]))
if (is.na(st_crs(eco))) st_crs(eco) <- 4326

if (is.null(strata_geo)) {
  message("Strata geometry absent: Figure 1 skipped.")
} else {
  str1 <- strata_geo
  if (is.na(st_crs(str1))) st_crs(str1) <- 4326
  str1 <- str1 %>% filter(spatial_unit %in% STRATA_KEEP) %>%
    group_by(spatial_unit) %>% summarise(.groups = "drop")
  miss <- setdiff(STRATA_KEEP, str1$spatial_unit)
  if (length(miss)) warning("Figure 1: strata absent from the layer: ", paste(miss, collapse = ", "))

  str_lab <- st_point_on_surface(str1)
  xy  <- st_coordinates(str_lab)
  nud <- left_join(tibble(spatial_unit = str_lab$spatial_unit), LAB_NUDGE, by = "spatial_unit") %>%
    mutate(dx = coalesce(dx, 0), dy = coalesce(dy, 0))
  st_geometry(str_lab) <- st_sfc(lapply(seq_len(nrow(xy)), function(i)
    st_point(c(xy[i, 1] + nud$dx[i], xy[i, 2] + nud$dy[i]))), crs = st_crs(str_lab))

  per <- raw %>% group_by(period) %>%
    summarise(st = n_distinct(stomach_id), se = n_distinct(paste(year, vessel.code, set)), .groups = "drop") %>%
    mutate(lab = sprintf("%s: %s stomachs, %s sets", period,
                         trimws(format(st, big.mark = ",")), trimws(format(se, big.mark = ","))))
  pts1 <- raw %>% filter(!is.na(longitude), !is.na(latitude)) %>%
    distinct(period, year, vessel.code, set, longitude, latitude) %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

  bb <- st_bbox(st_transform(st_as_sfc(st_bbox(c(xmin = -66.2, xmax = -60, ymin = 45.5, ymax = 49.2),
                                               crs = st_crs(4326))), CRS_MAP))
  land <- if (requireNamespace("rnaturalearth", quietly = TRUE))
    tryCatch(rnaturalearth::ne_states(country = c("canada", "united states of america"), returnclass = "sf"),
             error = function(e) NULL) else NULL

  p <- ggplot()
  if (!is.null(land)) p <- p + geom_sf(data = land, fill = "grey90", colour = "grey70", linewidth = 0.2)
  p <- p +
    geom_sf(data = eco, aes(fill = Area), alpha = 0.18, colour = NA) +
    geom_sf(data = str1, fill = NA, colour = STRATA_COL, linewidth = 0.3) +
    geom_sf(data = eco, aes(colour = Area), fill = NA, linewidth = 0.8) +
    geom_sf(data = pts1, colour = "black", size = 0.5, alpha = 0.7) +
    geom_sf_text(data = str_lab, aes(label = spatial_unit), size = 2.3,
                 colour = STRATA_COL, fontface = "bold") +
    facet_wrap(~period, ncol = 2, labeller = labeller(period = setNames(per$lab, per$period))) +
    scale_fill_manual(values = ECO_COL, limits = names(ECO_COL)) +
    scale_colour_manual(values = ECO_COL, limits = names(ECO_COL), guide = "none") +
    scale_x_continuous(breaks = seq(-65, -61, 1)) + scale_y_continuous(breaks = seq(46, 49, 1)) +
    coord_sf(crs = CRS_MAP, xlim = bb[c("xmin", "xmax")], ylim = bb[c("ymin", "ymax")], expand = FALSE) +
    labs(x = NULL, y = NULL) +
    guides(fill = guide_legend(override.aes = list(alpha = 0.35))) +
    theme_p(10) + theme(axis.text = element_text(size = 6), panel.grid = element_line(colour = "grey93"),
                        panel.spacing.x = unit(1.2, "lines"), strip.text = element_text(size = 9))
  if (requireNamespace("ggspatial", quietly = TRUE))
    p <- p + ggspatial::annotation_scale(
      data = data.frame(period = "2018-2019"), location = "br", width_hint = 0.12, text_cex = 0.45,
      height = unit(0.12, "cm"), bar_cols = c("grey20", "white"), line_width = 0.3,
      pad_x = unit(0.15, "cm"), pad_y = unit(0.15, "cm"))
  save_p(p, "Figure_1", 11, 5.4, per)
}

ncat <- read_csv("Output_Appendices/TableA0_prey_categories.csv", show_col_types = FALSE)
p <- ggplot(ncat, aes(x_threshold, n_cat)) +
  geom_col(fill = "#4E6470", width = 7) +
  scale_x_continuous(breaks = seq(0, 1000, 100)) +
  scale_y_continuous(breaks = seq(0, 200, 20), expand = expansion(mult = c(0, 0.03))) +
  labs(x = "Aggregation threshold x (stomachs)", y = "Prey categories retained") +
  theme(panel.grid.major.x = element_blank())
save_p(p, "Figure_2", 8, 4.2, ncat)

plot_tree <- function() {
  txt <- c(Stability = "white", Substitution = "white", "Functional change" = "#14232B", Reorganization = "white")
  sm <- function(x, y, fam, state, sig) list(
    annotate("rect", xmin = x - 10.2, xmax = x + 10.2, ymin = y - 6.6, ymax = y + 6.6,
             fill = PAL[[fam]], colour = "#1F2D3D", linewidth = .35),
    annotate("text", x = x, y = y + 3.6, label = fam, fontface = "bold", size = 3.05, colour = txt[[fam]]),
    annotate("text", x = x, y = y + 0.2, label = state, size = 2.55, colour = txt[[fam]], lineheight = .9),
    annotate("text", x = x, y = y - 4.0, label = sig, size = 2.7, colour = txt[[fam]], fontface = "italic"))
  arr <- function(x, y, xend, yend, colour = "#4A5A66")
    annotate("segment", x = x, y = y, xend = xend, yend = yend, colour = colour, linewidth = .5,
             arrow = arrow(length = unit(2.2, "mm"), type = "closed"))
  ggplot() + xlim(0, 100) + ylim(-5, 65) + theme_void() +
    annotate("rect", xmin = 20, xmax = 76, ymin = 56.2, ymax = 63.8, fill = "white", colour = "#1F2D3D", linewidth = .5) +
    annotate("text", x = 48, y = 61.4, label = "Classifiable cell?", fontface = "bold", size = 3.7) +
    annotate("text", x = 48, y = 58.2, label = "At least 5 stomachs with prey per period; all three tests computable",
             size = 2.6, colour = "#3C4B57") +
    annotate("rect", xmin = 83, xmax = 99, ymin = 56.2, ymax = 63.8, fill = "#E3E7EA", colour = "#8A9BA8", linewidth = .4) +
    annotate("text", x = 91, y = 61.0, label = "Inconclusive", fontface = "bold", size = 3.0, colour = "#3C4B57") +
    annotate("text", x = 91, y = 58.4, label = "excluded from\npercentages", size = 2.5, colour = "#3C4B57", lineheight = .9) +
    arr(76, 60, 82.6, 60) + annotate("text", x = 79.3, y = 62.1, label = "no", size = 2.7, colour = "#3C4B57") +
    arr(48, 56.2, 48, 50.6) + annotate("text", x = 50.4, y = 53.4, label = "yes", size = 2.7, colour = "#3C4B57") +
    annotate("rect", xmin = 20, xmax = 76, ymin = 43, ymax = 50.6, fill = "white", colour = "#1F2D3D", linewidth = .5) +
    annotate("text", x = 48, y = 48.2, label = "Composition test (C)", fontface = "bold", size = 3.7) +
    annotate("text", x = 48, y = 45.1, label = "PERMANOVA on Bray-Curtis dissimilarities; trawl set as replicate",
             size = 2.6, colour = "#3C4B57") +
    arr(34, 43, 24, 39.2, "#2E7D32") +
    annotate("text", x = 13.5, y = 41.2, label = "C not significant", size = 3, colour = "#2E7D32", fontface = "bold") +
    arr(62, 43, 76, 39.2, "#B3261E") +
    annotate("text", x = 87, y = 41.2, label = "C significant (p < 0.05)", size = 3, colour = "#B3261E", fontface = "bold") +
    annotate("rect", xmin = 1.5, xmax = 48.5, ymin = -1.5, ymax = 39, fill = "#F5F7F5", colour = "#2E7D32", linewidth = .5) +
    annotate("rect", xmin = 51.5, xmax = 98.5, ymin = -1.5, ymax = 39, fill = "#FBF3F2", colour = "#B3261E", linewidth = .5) +
    annotate("text", x = 25, y = 36.6, label = "Composition not significantly different", fontface = "bold", size = 3.1, colour = "#2E7D32") +
    annotate("text", x = 75, y = 36.6, label = "Composition significantly different", fontface = "bold", size = 3.1, colour = "#B3261E") +
    sm(13.5, 26, "Stability", "Stable Diet;\nEmerging Shift*", "H′ ns · Bs ns") +
    sm(36.5, 26, "Functional change", "Internal Rebalancing", "H′ sig · Bs ns") +
    sm(13.5, 10, "Functional change", "Niche Compression/\nExpansion", "H′ ns · Bs sig") +
    sm(36.5, 10, "Functional change", "Niche Restructuring", "H′ sig · Bs sig") +
    sm(63.5, 26, "Substitution", "Ghost Shift", "H′ ns · Bs ns") +
    sm(86.5, 26, "Reorganization", "Partial Diet Shift", "H′ sig · Bs ns") +
    sm(63.5, 10, "Reorganization", "Structural Shift", "H′ ns · Bs sig") +
    sm(86.5, 10, "Reorganization", "Major Shift", "H′ sig · Bs sig") +
    annotate("text", x = 1.5, y = -4.2, hjust = 0, size = 2.45, colour = "#3C4B57", lineheight = .95,
             label = paste0("* Emerging Shift: composition trend (0.05 ≤ p < 0.10) with H′ and Bs non-significant; ",
                            "classified as Stability.\nsig, p < 0.05; ns, p ≥ 0.05."))
}
save_p(plot_tree(), "Figure_3", 9.5, 6.4)

d4 <- filter(ref, level == "all_gulf")
save_p(regime_plot(d4), "Figure_4", 9, 5, d4)

save_p(units_plot("P1b"), "Figure_5", 9.5, 9.4, filter(units, contrast == "P1b"))
save_p(units_plot("PTa"), "Figure_6", 9.5, 9.4, filter(units, contrast == "PTa"))

cat("\nSupplementary figures\n")

RANKS <- c("class", "order", "family", "genus", "species")
idd <- read_csv("Output_Appendices/TableA7_identification_depth.csv", show_col_types = FALSE) %>%
  filter(rank %in% RANKS) %>%
  select(rank, starts_with("pct_at_or_finer")) %>%
  pivot_longer(-rank, names_to = "period", values_to = "pct", names_prefix = "pct_at_or_finer_") %>%
  mutate(rank = factor(rank, RANKS, c("Class", "Order", "Family", "Genus", "Species")),
         period = factor(period, c("2004-2006", "2018-2019")))
p <- ggplot(idd, aes(rank, pct, colour = period, group = period)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.4) +
  geom_text(aes(label = sprintf("%.1f", pct),
                vjust = ifelse(period == "2004-2006", ifelse(as.integer(rank) <= 2, 1.9, -1.1),
                               ifelse(as.integer(rank) <= 2, -1.1, 1.9))),
            size = 3, show.legend = FALSE) +
  scale_colour_manual(values = c("2004-2006" = "#2F4A5A", "2018-2019" = "#D96C4A")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(x = "Taxonomic rank (identified at this rank or finer)", y = "Ranked prey records (%)")
save_p(p, "Figure_S1", 6.5, 4.2, idd)

dS3 <- filter(ref, level == "ecoregion")
save_p(regime_plot(dS3), "Figure_S2", 9, 5, dS3)

eco_sd <- res_all %>% filter(level == "ecoregion") %>%
  spread_by_resolution(keys = c("currency", "contrast", "spatial_unit")) %>%
  clean_fam() %>%
  transmute(mode = currency, contrast = factor(contrast, CT), spatial_unit,
            family = factor(family, FAM), pct_mean, sd_across_x = pct_sd)
n_eco <- n_class %>% filter(level == "ecoregion") %>% group_by(spatial_unit) %>%
  summarise(n = max(n_classifiable), .groups = "drop")
dS4 <- eco_sd %>% left_join(n_eco, by = "spatial_unit") %>%
  mutate(unit = factor(paste0(ECO_FULL[spatial_unit], " (n = ", n, ")")))
p <- regime_plot(dS4, labels = FALSE, base = 10, pt = 2) + facet_wrap(~unit, ncol = 2) +
  theme(axis.text.x = element_text(size = 7))
save_p(p, "Figure_S3", 11, 8, dS4)

save_p(units_plot("P1a"), "Figure_S4", 9.5, 9.4, filter(units, contrast == "P1a"))
save_p(units_plot("P2"),  "Figure_S5", 9.5, 9.4, filter(units, contrast == "P2"))
save_p(units_plot("PTb"), "Figure_S6", 9.5, 9.4, filter(units, contrast == "PTb"))

dS8 <- filter(ref, level == "stratum")
save_p(regime_plot(dS8), "Figure_S7", 9, 5, dS8)

p <- map_plot("biomass");    if (!is.null(p)) save_p(p, "Figure_S8", 12, 8.8) else message("Strata geometry absent: maps skipped.")
p <- map_plot("occurrence"); if (!is.null(p)) save_p(p, "Figure_S9", 12, 8.8)

dS10 <- ref %>% filter(contrast == "PTa") %>%
  mutate(Scale = factor(LL[level], levels = LL), Currency = factor(ML[mode], levels = ML))
p <- ggplot(dS10, aes(Scale, pct_mean, colour = family, shape = family, group = family)) +
  geom_errorbar(aes(ymin = pmax(pct_mean - sd_across_x, 0), ymax = pct_mean + sd_across_x),
                width = 0.1, linewidth = 0.4, alpha = 0.7) +
  geom_line(linewidth = 0.9) + geom_point(size = 3.2) +
  facet_wrap(~Currency) +
  scale_colour_manual(values = PAL) + scale_shape_manual(values = SHAPE) +
  scale_x_discrete(limits = rev(levels(dS10$Scale))) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(x = "Spatial unit, from finest to coarsest", y = "Family frequency (%)") +
  theme(legend.position = "top")
save_p(p, "Figure_S10", 9.5, 5.4, dS10)

for (k in seq_along(LV)) {
  d <- sweep_data(LV[k])
  save_p(sweep_plot(d), paste0("Figure_S", 10 + k), if (k == 1) 10 else 8.4, 4.8, d)
}

RULES <- c("_1" = "Pooled (primary)", "_2" = "Pooled cross-predator", "_PP" = "Per-predator")
RULE_COL <- c("Pooled (primary)" = "#2F4A5A", "Pooled cross-predator" = "#7F9CAB", "Per-predator" = "#D96C4A")
gr <- bind_rows(lapply(names(RULES), function(pf) {
  r <- suppressWarnings(read_all_runs(levels = "all_gulf",
                                      dirs = c(all_gulf = paste0("Sensitivity_all_gulf",
                                                                if (pf == "_1") "_1_T" else pf))))
  family_by_resolution(r, keys = c("currency", "contrast")) %>%
    group_by(currency, contrast, family) %>%
    summarise(pct_mean = mean(pct), pct_sd = sd(pct), n_thresholds = n_distinct(x_threshold),
              .groups = "drop") %>%
    mutate(Rule = RULES[[pf]])
})) %>%
  clean_fam() %>%
  filter(!is.na(contrast)) %>%
  mutate(Rule = factor(Rule, levels = unname(RULES)), family = factor(family, FAM),
         x = match(as.character(contrast), CT), Currency = factor(ML[currency], levels = ML))
pd <- position_dodge(width = 0.5)
p <- ggplot(gr, aes(x, pct_mean, colour = Rule, group = Rule)) +
  annotate("rect", xmin = 0.5, xmax = 3.5, ymin = -Inf, ymax = Inf, fill = "#D7E6EC", alpha = 0.55) +
  annotate("rect", xmin = 3.5, xmax = 5.5, ymin = -Inf, ymax = Inf, fill = "#F6E0D6", alpha = 0.65) +
  geom_vline(xintercept = 3.5, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_errorbar(aes(ymin = pmax(pct_mean - pct_sd, 0), ymax = pct_mean + pct_sd),
                width = 0.15, linewidth = 0.4, position = pd) +
  geom_point(size = 2.2, position = pd) +
  facet_grid(family ~ Currency, scales = "free_y") +
  scale_colour_manual(values = RULE_COL) +
  scale_x_continuous(breaks = 1:5, labels = unname(CL2[CT]), expand = expansion(add = 0.25)) +
  labs(x = NULL, y = "Family frequency (%)") +
  theme(legend.position = "top", axis.text.x = element_text(size = 7.5))
save_p(p, "Figure_7", 9, 9.5, gr)

dS16 <- cells %>% filter(contrast == "PTa") %>%
  mutate(Scale = factor(LL[level], levels = LL), Currency = factor(ML[mode], levels = ML))
p <- ggplot(dS16, aes(support)) +
  stat_ecdf(geom = "step", colour = "#2F4A5A", linewidth = .8) +
  geom_vline(xintercept = .8, linetype = 2, colour = "#D96C4A") +
  facet_grid(Currency ~ Scale) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
  labs(x = "Modal-family support across 100 resolutions", y = "Cumulative fraction of cells")
save_p(p, "Figure_S14", 9, 5.2)

dS18 <- map_dfr(LV, rr, file = "M4_alpha_trend_sensitivity.csv") %>%
  filter(contrast == "PTa", abs(trend - .1) < 1e-9) %>%
  mutate(Scale = factor(LL[level], levels = LL), Currency = factor(ML[mode], levels = ML),
         family = factor(family, FAM))
p <- ggplot(dS18, aes(alpha, pct_mean, colour = family)) +
  geom_line() + geom_point(size = 1.4) +
  facet_grid(Currency ~ Scale) +
  scale_colour_manual(values = PAL) +
  labs(x = "Significance threshold", y = "Family frequency (%)")
save_p(p, "Figure_S15", 9, 5.2, dS18)

boot <- map_dfr(LV, rr, file = "M3_bootstrap_CI_family_pct.csv") %>%
  filter(family %in% c("Stability", "Substitution")) %>%
  mutate(Scale = factor(LL[level], levels = LL), Currency = factor(ML[mode], levels = ML),
         panel = factor(paste(Scale, Currency, family, sep = " | "),
                        levels = unlist(lapply(LL, function(sc) unlist(lapply(ML, function(m)
                          paste(sc, m, c("Stability", "Substitution"), sep = " | ")))))),
         Contrast = factor(CL1[contrast], levels = unname(CL1)))
env <- boot %>% group_by(panel, contrast_type) %>%
  summarise(lo = min(ci_lo), hi = max(ci_hi), first = min(as.integer(droplevels(Contrast))),
            .groups = "drop")
nct <- boot %>% group_by(panel) %>% summarise(n = n_distinct(contrast), .groups = "drop")
env <- env %>% left_join(nct, by = "panel") %>%
  mutate(xmin = ifelse(contrast_type == "within", 0.6, 3.6),
         xmax = ifelse(contrast_type == "within", 3.4, n + 0.4))
gaps <- env %>% select(panel, contrast_type, lo, hi) %>%
  pivot_wider(names_from = contrast_type, values_from = c(lo, hi)) %>%
  mutate(family = sub(".* \\| ", "", panel),
         gap = ifelse(family == "Stability", lo_within - hi_between, lo_between - hi_within),
         lab = sprintf("gap %+.1f", gap))
p <- ggplot(boot) +
  geom_pointrange(aes(Contrast, pct_mean, ymin = ci_lo, ymax = ci_hi, colour = contrast_type),
                  size = 0.25, linewidth = 0.5) +
  geom_rect(data = env, aes(xmin = xmin, xmax = xmax, ymin = lo, ymax = hi, fill = contrast_type), alpha = 0.18) +
  geom_text(data = gaps, aes(x = 2, y = 97, label = lab), size = 2.7, colour = "grey20") +
  facet_wrap(~panel, ncol = 4, scales = "free_x") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  scale_colour_manual(values = TYPE_COL, labels = TYPE_LAB, breaks = names(TYPE_LAB)) +
  scale_fill_manual(values = TYPE_COL, labels = TYPE_LAB, breaks = names(TYPE_LAB)) +
  labs(x = NULL, y = "Family frequency (%)") +
  theme_p(9) + theme(axis.text.x = element_text(size = 6.3, angle = 40, hjust = 1),
                     strip.text = element_text(size = 7.5))
save_p(p, "Figure_S16", 10, 8, gaps)

dS20 <- read_csv("NullCalibration_all_gulf_1/null_fpr_by_sample_size.csv", show_col_types = FALSE) %>%
  pivot_longer(-c(mode, n_bin, n), names_to = "Test", values_to = "Rate") %>%
  mutate(Test = recode(Test, "Composition (PERMANOVA)" = "Composition",
                       "Diversite (H')" = "Diversity", "Amplitude (Bs)" = "Niche breadth"),
         Currency = factor(ML[mode], levels = ML), n_bin = factor(n_bin, levels = unique(n_bin)))
p <- ggplot(dS20, aes(n_bin, Rate, colour = Test, group = Test)) +
  geom_hline(yintercept = 5, linetype = 2, colour = "grey50") +
  geom_line() + geom_point() +
  facet_wrap(~Currency, nrow = 1) +
  scale_colour_manual(values = c(Composition = "#397D3B", Diversity = "#2864A1", "Niche breadth" = "#C45C3B")) +
  labs(x = "Sets per period (smaller count)", y = "Null rejection rate (%)")
save_p(p, "Figure_S17", 9, 4.2, dS20)

balanced_plot <- function(dir, drop = character(0)) {
  rep_d <- read_csv(file.path(dir, "balanced_family_per_replicate.csv"), show_col_types = FALSE) %>%
    clean_fam() %>% filter(family %in% c("Stability", "Substitution"), !contrast %in% drop) %>%
    mutate(contrast = factor(contrast, setdiff(CT, drop)), Currency = factor(ML[mode], levels = ML))
  unb <- read_csv(file.path(dir, "balanced_vs_unbalanced.csv"), show_col_types = FALSE) %>%
    clean_fam() %>% filter(family %in% c("Stability", "Substitution"), !contrast %in% drop) %>%
    mutate(contrast = factor(contrast, setdiff(CT, drop)), Currency = factor(ML[mode], levels = ML))
  list(data = unb, plot = ggplot(rep_d, aes(contrast, pct, fill = contrast_type)) +
         geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.3) +
         geom_point(position = position_jitter(width = .1, height = 0, seed = 416), size = 0.9, alpha = .7) +
         geom_point(data = unb, aes(contrast, pct_unbalanced), inherit.aes = FALSE, shape = 4, size = 2.4, stroke = 0.9) +
         facet_grid(Currency ~ family) +
         scale_fill_manual(values = c(within = "#8DB1C4", between = "#E5A08C"),
                           labels = TYPE_LAB, breaks = names(TYPE_LAB)) +
         scale_x_discrete(labels = CL2) +
         scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
         labs(x = NULL, y = "Family frequency (%)") +
         theme(axis.text.x = element_text(size = 7.5)))
}
b <- balanced_plot("Robustness_all_gulf_1_BAL");  save_p(b$plot, "Figure_S18", 8, 5.6, b$data)
b <- balanced_plot("Robustness_ecoregion_1_BAL"); save_p(b$plot, "Figure_S19", 8, 5.6, b$data)

cat("\nDone: ", length(list.files(OUT, pattern = "^Figure_")), " figures in ", OUT, "\n", sep = "")

})
