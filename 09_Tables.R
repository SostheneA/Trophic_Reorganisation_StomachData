# Main manuscript and supplementary tables
# Input: data/dat_classed.rda; robustness summaries
# Output: Output_Paper_Tables/

if (!file.exists("data/dat_classed.rda")) stop("Missing prepared data")
local({
  suppressPackageStartupMessages({library(dplyr);library(sf);library(tidyr);library(readr)})
  e <- new.env(parent = emptyenv())
  load("data/dat_classed.rda", envir = e)
  raw <- if (inherits(e$dat_classed, "sf")) sf::st_drop_geometry(e$dat_classed) else e$dat_classed

  # Predator characteristics and sample sizes
  s1 <- tibble::tribble(
    ~predator, ~sci, ~len_p1, ~len_p2, ~lmat, ~lmax, ~source,
    "Atlantic cod", "Gadus morhua", "5-98", "4-83", "37", "101", "Swain et al. 2019",
    "American plaice", "Hippoglossoides platessoides", "5-50", "3-49", "30", "54", "Ricard et al. 2016",
    "Yellowtail flounder", "Myzopsetta ferruginea", "13.6-34", "12-31", "25", "36", "Surette and Swain 2016",
    "Thorny skate", "Amblyraja radiata", "8.7-56", "10-56", "53", "56", "Templeman 1987",
    "Capelin", "Mallotus villosus", "8-18", "2-18", "15", "18", "Froese and Pauly 2025",
    "Daubed shanny", "Leptoclinus maculatus", "7-15", "6-15", "13", "18", "Froese and Pauly 2025",
    "Moustache sculpin", "Triglops murrayi", "3-16", "4-17", "-", "17", "Froese and Pauly 2025",
    "Longhorn sculpin", "Myoxocephalus octodecemspinosus", "8-37", "4-37", "24", "37", "Beacham 1983",
    "Winter flounder", "Pseudopleuronectes americanus", "9.2-45", "8-47", "25", "47", "Scharf et al. 2000",
    "White hake", "Urophycis tenuis", "5-75", "7-68", "36", "75", "Turcotte et al. 2025",
    "Atlantic herring", "Clupea harengus", "3.3-36", "9.5-36", "26", "36", "Beaudry-Sylvestre et al. 2022",
    "Alligatorfish", "Aspidophoroides monopterygius", "5-19", "4-19", "10", "19", "Froese and Pauly 2025",
    "Atlantic hookear sculpin", "Artediellus atlanticus", "3.5-9", "3-8", "-", "9", "Froese and Pauly 2025",
    "Atlantic poacher", "Leptagonus decagonus", "6.7-24", "2-24", "-", "24", "Froese and Pauly 2025",
    "Witch flounder", "Glyptocephalus cynoglossus", "15-46", "20-51", "30", "52", "Swain et al. 2012",
    "Dusky snailfish", "Liparis gibbus", "3-31.6", "2-32", "-", "32", "Froese and Pauly 2025",
    "Spatulate sculpin", "Icelus spatula", "3.9-14", "3-14", "-", "14", "Froese and Pauly 2025",
    "Laval's eelpout", "Lycodes lavalaei", "8-67", "7-57", "-", "68", "Froese and Pauly 2025",
    "Greenland halibut", "Reinhardtius hippoglossoides", "6.1-62", "8-62", "40", "73", "Benoît et al. 2026",
    "Arctic staghorn sculpin", "Gymnocanthus tricuspis", "3.4-24", "5-25", "6", "27", "Froese and Pauly 2025",
    "Stout eelblenny", "Anisarchus medius", "6.4-17", "7-17", "-", "31", "Froese and Pauly 2025",
    "Alewife", "Alosa pseudoharengus", "10-29", "8-29", "22", "29", "Breau and Gibson 2024",
    "Fourline snakeblenny", "Eumesogrammus praecisus", "6-24", "8-24", "-", "24", "Froese and Pauly 2025",
    "Rainbow smelt", "Osmerus mordax", "6.8-24", "4-23", "11", "25", "Scott and Scott 1988",
    "Shorthorn sculpin", "Myoxocephalus scorpius", "6-43", "3-43", "16", "44", "Froese and Pauly 2025",
    "Redfish unidentified", "Sebastes spp.", "6.1-47", "3.1-44", "19", "48", "Brûlé et al. 2024",
    "Atlantic spiny lumpsucker", "Eumicrotremus spinosus", "2.7-10", "2-10", "-", "10", "Froese and Pauly 2025",
    "Smooth skate", "Malacoraja senta", "9-55", "13-51", "54", "68", "Sulikowski et al. 2009",
    "Snakeblenny", "Lumpenus lampretaeformis", "6.9-43", "12-35", "20", "43", "Froese and Pauly 2025",
    "Arctic alligatorfish", "Aspidophoroides olrikii", "4.5-9", "5-8", "-", "9", "Froese and Pauly 2025",
    "Fourbeard rockling", "Enchelyopus cimbrius", "6.6-28", "8-26", "25", "30", "Froese and Pauly 2025",
    "Marlin-spike grenadier", "Nezumia bairdii", "5.5-32.9", "7.8-31", "-", "33", "Froese and Pauly 2025",
    "Three-spined stickleback", "Gasterosteus aculeatus", "2.5-8", "3-8", "5", "8", "Froese and Pauly 2025",
    "Atlantic mackerel", "Scomber scombrus", "14.5-31", "9-33.5", "26", "45", "Lehoux et al. 2025",
    "Greenland cod", "Gadus macrocephalus", "9-50", "23-53", "46", "58", "Froese and Pauly 2025",
    "Arctic hookear sculpin", "Artediellus uncinatus", "4.1-7.8", "3.2-8", "-", "12", "Froese and Pauly 2025",
    "Sea raven", "Hemitripterus americanus", "5-51", "5-41", "-", "51", "Froese and Pauly 2025",
    "Vahl's eelpout", "Lycodes vahlii", "6.4-44", "16-44", "-", "46", "Froese and Pauly 2025",
    "Northern sand lance", "Ammodytes dubius", "9-24", "7-23", "-", "24", "Froese and Pauly 2025",
    "Windowpane flounder", "Scophthalmus aquosus", "12.7-24.7", "11-64", "22", "64", "Froese and Pauly 2025"
  )

  key <- function(x) tolower(trimws(gsub("[’`]", "'", as.character(x))))
  fmt <- function(x) trimws(format(x, big.mark = ","))

  cand <- names(raw)[vapply(raw, function(x) is.character(x) || is.factor(x), logical(1))]
  hits <- sapply(cand, function(cn) {
    u <- unique(key(raw[[cn]]))
    c(common = sum(key(s1$predator) %in% u), sci = sum(key(s1$sci) %in% u))
  })
  print(head(t(hits)[order(-colSums(hits)), , drop = FALSE], 5))
  if (max(hits) == 0) stop("No column of `raw` matches the Table S1 predator names.")

  best     <- which(hits == max(hits), arr.ind = TRUE)[1, ]
  PRED_COL <- colnames(hits)[best["col"]]
  BY       <- rownames(hits)[best["row"]]
  message("Predator column: ", PRED_COL, " (matched on ", BY, " names, ",
          max(hits), " of ", nrow(s1), ")")

  n_tab <- raw %>%
    mutate(k = key(.data[[PRED_COL]])) %>%
    group_by(k, period) %>%
    summarise(n = n_distinct(stomach_id),
              sets = n_distinct(paste(year, vessel.code, set)), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = period, values_from = c(n, sets), values_fill = 0)

  tab_s1 <- s1 %>%
    mutate(k = key(if (BY == "common") predator else sci)) %>%
    left_join(n_tab, by = "k")

  miss  <- tab_s1$predator[is.na(tab_s1$`n_2004-2006`)]
  extra <- setdiff(n_tab$k, tab_s1$k)
  if (length(miss))  warning("Table S1: no stomachs matched for: ", paste(miss, collapse = ", "))
  if (length(extra)) message("In the data but not in Table S1: ", paste(extra, collapse = ", "))
  tab_s1 <- tab_s1 %>% mutate(n_total = `n_2004-2006` + `n_2018-2019`)
  low <- tab_s1 %>% filter(n_total < 90 | `n_2004-2006` < 10 | `n_2018-2019` < 10)
  if (nrow(low)) warning("Below the 90 / 10 threshold: ", paste(low$predator, collapse = ", "))

  out_s1 <- tab_s1 %>%
    transmute(Predator = predator, `Scientific name` = sci,
              `Length range (n) 2004-2006` = sprintf("%s (%s)", len_p1, fmt(`n_2004-2006`)),
              `Length range (n) 2018-2019` = sprintf("%s (%s)", len_p2, fmt(`n_2018-2019`)),
              `Lmat (cm)` = lmat, `Lmax (cm)` = lmax, Source = source) %>%
    bind_rows(tibble(Predator = "Total",
                     `Length range (n) 2004-2006` = sprintf("(%s)", fmt(sum(tab_s1$`n_2004-2006`, na.rm = TRUE))),
                     `Length range (n) 2018-2019` = sprintf("(%s)", fmt(sum(tab_s1$`n_2018-2019`, na.rm = TRUE)))))

  dir.create("Output_Paper_Tables", recursive = TRUE, showWarnings = FALSE)
  write.csv(out_s1, "Output_Paper_Tables/Table_S1.csv", row.names = FALSE, na = "")
  write.csv(tab_s1 %>% select(predator, starts_with("n_"), starts_with("sets_")),
            "Output_Paper_Tables/Table_S1_counts.csv", row.names = FALSE)
  print(out_s1, n = Inf)

  library(flextable)

  w <- tab_s1 %>%
    transmute(Predator = predator, sci,
              len_p1, n1 = fmt(`n_2004-2006`),
              len_p2, n2 = fmt(`n_2018-2019`),
              lmat, lmax, source) %>%
    bind_rows(tibble(Predator = "Total", sci = "",
                     len_p1 = "", n1 = fmt(sum(tab_s1$`n_2004-2006`, na.rm = TRUE)),
                     len_p2 = "", n2 = fmt(sum(tab_s1$`n_2018-2019`, na.rm = TRUE)),
                     lmat = "", lmax = "", source = "")) %>%
    mutate(p1 = "", p2 = "")

  ft <- flextable(w, col_keys = c("Predator", "sci", "p1", "p2", "lmat", "lmax", "source")) %>%
    compose(j = "p1", value = as_paragraph(len_p1, " (", as_b(n1), ")")) %>%
    compose(j = "p2", value = as_paragraph(len_p2, " (", as_b(n2), ")")) %>%
    set_header_labels(sci = "Scientific name",
                      p1 = "Length range (n)\n2004-2006", p2 = "Length range (n)\n2018-2019",
                      lmat = "Lmat\n(cm)", lmax = "Lmax\n(cm)", source = "Source") %>%
    italic(j = "sci", part = "body") %>%
    bold(part = "header") %>%
    bold(i = nrow(w), part = "body") %>%
    font(fontname = "Arial", part = "all") %>%
    fontsize(size = 9, part = "all") %>%
    align(j = c("p1", "p2", "lmat", "lmax"), align = "center", part = "all") %>%
    autofit()

  save_as_docx(ft, path = "Output_Paper_Tables/Table_S1.docx",
               pr_section = officer::prop_section(
                 page_size = officer::page_size(orient = "landscape")))
  message("Saved: ", normalizePath("Output_Paper_Tables/Table_S1.docx"))
})

# Publication tables
local({
try(suppressWarnings(Sys.setlocale("LC_ALL", "C.UTF-8")), silent = TRUE)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr)
  library(flextable); library(officer)
})

PREY_FAMILY <- "_1"
source("R_helpers/Results.R", local = TRUE)

OUT <- "Output_Paper_Tables"
dir.create(file.path(OUT, "csv"), recursive = TRUE, showWarnings = FALSE)

FAM <- c("Stability", "Substitution", "Functional change", "Reorganization")
CT  <- c("P1a", "P1b", "P2", "PTb", "PTa")
CTL <- c(P1a = "2004 / 2006", P1b = "2004-05 / 2006", P2 = "2018 / 2019",
         PTb = "2006 / 2018", PTa = "2004-06 / 2018-19")
LV  <- c("all_gulf", "ecoregion", "stratum")
LL  <- c(all_gulf = "sGSL-wide", ecoregion = "Ecoregion", stratum = "Stratum")
ML  <- c(biomass = "Biomass", occurrence = "Occurrence")
ECO <- c(Central = "Central sGSL", "Chaleur Bay" = "Chaleur Bay",
         "Magdalen Shallows" = "Magdalen Shallows", Northumberland = "Northumberland Strait")
f1  <- function(x, d = 1) formatC(x, format = "f", digits = d)
msd <- function(m, s) paste0(f1(m), " (", f1(s), ")")

clean_fam <- function(d) {
  names(d) <- sub("Reorganisation", "Reorganization", names(d), fixed = TRUE)
  for (v in intersect(c("family", "modal_family"), names(d)))
    d[[v]] <- sub("Reorganisation", "Reorganization", d[[v]], fixed = TRUE)
  d
}
rd <- function(f) clean_fam(read_csv(f, show_col_types = FALSE))
rr <- function(level, file) rd(sprintf("Robustness_%s_1_T/%s", level, file)) %>% mutate(level = level)

norm_ct <- function(x) {
  x <- gsub("–", "-", x, fixed = TRUE)
  key <- c("2004 vs 2006" = "P1a", "2004-2005 vs 2006" = "P1b", "2018 vs 2019" = "P2",
           "2006 vs 2018" = "PTb", "2004-2006 vs 2018-2019" = "PTa")
  unname(key[x])
}

ft_style <- function(d) {
  flextable(d) %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 9, part = "all") %>%
    bold(part = "header") %>%
    bg(bg = "#E9EEF1", part = "header") %>%
    border_remove() %>%
    hline_top(border = fp_border(width = 1), part = "header") %>%
    hline_bottom(border = fp_border(width = 0.75), part = "header") %>%
    hline_bottom(border = fp_border(width = 1), part = "body") %>%
    align(align = "center", part = "all") %>%
    align(j = 1, align = "left", part = "all") %>%
    padding(padding = 2.5, part = "all") %>%
    autofit() %>%
    fit_to_width(max_width = 6.5)
}

TABLES <- list()
add_table <- function(label, caption, d) {
  TABLES[[length(TABLES) + 1]] <<- list(label = label, caption = caption, data = d)
  write_csv(d, file.path(OUT, "csv", paste0(gsub(" ", "_", label), ".csv")), na = "")
  cat("  ", label, " (", nrow(d), " rows)\n", sep = "")
}

ref <- map_dfr(LV, rr, file = "M0_family_frequencies_reference.csv") %>%
  mutate(family = factor(family, FAM), contrast = factor(contrast, CT))

cells <- map_dfr(LV, rr, file = "M5_cell_classification_stability.csv")
n_class <- cells %>% mutate(spatial_unit = sub(".*\\|", "", cell)) %>%
  count(level, contrast, mode, spatial_unit, name = "n")

units <- rd("Output_Figures/data_Fig6_9_unit_families.csv") %>%
  mutate(contrast = norm_ct(contrast), mode = currency) %>%
  left_join(n_class, by = c("level", "contrast", "mode", "spatial_unit")) %>%
  mutate(n = replace_na(n, 0L))

freq_wide <- function(d) d %>%
  arrange(mode, contrast, family) %>%
  mutate(value = msd(pct_mean, sd_across_x)) %>%
  transmute(Currency = ML[mode], Contrast = CTL[as.character(contrast)], family, value) %>%
  pivot_wider(names_from = family, values_from = value) %>%
  select(Currency, Contrast, all_of(FAM))

cat("\nTables\n")

add_table("Table 1",
  paste("Diagnostic typology. Each predator × size-class cell is scored on three tests, prey composition",
        "(C; PERMANOVA), diversity (H′; jackknifed Shannon) and niche breadth (Bs; standardized Levins),",
        "as significant (sig, p < 0.05) or non-significant (ns). A composition trend (0.05 ≤ p < 0.10) is",
        "distinguished only when H′ and Bs are both non-significant; otherwise the cell is classified on the",
        "significant signal."),
  tribble(
    ~`Diagnostic state`,           ~C,              ~`H′`, ~Bs,   ~Family,
    "Stable Diet",                 "ns (p ≥ 0.10)", "ns",  "ns",  "Stability",
    "Emerging Shift",              "trend",         "ns",  "ns",  "Stability",
    "Ghost Shift",                 "sig",           "ns",  "ns",  "Substitution",
    "Internal Rebalancing",        "ns",            "sig", "ns",  "Functional change",
    "Niche Compression/Expansion", "ns",            "ns",  "sig", "Functional change",
    "Niche Restructuring",         "ns",            "sig", "sig", "Functional change",
    "Partial Diet Shift",          "sig",           "sig", "ns",  "Reorganization",
    "Structural Shift",            "sig",           "ns",  "sig", "Reorganization",
    "Major Shift",                 "sig",           "sig", "sig", "Reorganization"))

add_table("Table 2",
  paste("sGSL-wide family frequencies (%). Values are means across 100 prey resolutions; parentheses give",
        "standard deviations across resolutions. Frequencies exclude Inconclusive cells. Absent families",
        "contribute zero. Independently rounded percentages need not sum to exactly 100%."),
  freq_wide(filter(ref, level == "all_gulf")))

add_table("Table 3",
  paste("Family frequencies (%) at the ecoregion and stratum scales (the values plotted in Figures S3 and S8).",
        "Percentages are averaged across classifiable spatial units and then across resolutions. Parentheses",
        "give standard deviations across resolutions. Spatial-unit means are unweighted."),
  bind_rows(lapply(c("ecoregion", "stratum"), function(l)
    freq_wide(filter(ref, level == l)) %>% mutate(Scale = LL[[l]], .before = 1))))

unit_rows <- function(lv, m) {
  ord <- if (lv == "ecoregion") names(ECO) else
    as.character(sort(as.integer(unique(units$spatial_unit[units$level == "stratum"]))))
  units %>% filter(level == lv, contrast == "PTa", mode == m) %>%
    mutate(rank = match(spatial_unit, ord)) %>% arrange(rank) %>%
    transmute(Unit = coalesce(unname(ECO[spatial_unit]), spatial_unit), n,
              across(all_of(FAM), ~ f1(.x)))
}
add_table("Table 4",
  paste("Full between-period (2004-2006 versus 2018-2019) family frequencies (%) by ecoregion, the values",
        "shown in the top panels of Figure 6. Values are means across 100 prey resolutions; n is the number",
        "of classifiable predator \u00d7 size-class cells. A, biomass; B, occurrence."),
  bind_rows(mutate(unit_rows("ecoregion", "biomass"),    Currency = "A) Biomass",    .before = 1),
            mutate(unit_rows("ecoregion", "occurrence"), Currency = "B) Occurrence", .before = 1)))

add_table("Table S1",
  paste("Predator taxa retained in the analysis, with the length at maturity (Lmat) used to separate juvenile",
        "and adult size classes (Scharf et al., 2000, and regional sources). Taxa without a published regional",
        "Lmat were analysed as a single size class."),
  rd("Output_Tables/TableS1_predators.csv"))

s2 <- rd("Output_Tables/TableS2_sampling.csv") %>% select(-any_of(c("Testable", "Reliable")))
s2 <- bind_rows(mutate(s2, across(everything(), as.character)),
                tibble(Predator = "Total", `Size class` = "",
                       `Stomachs 2004-2006` = format(sum(s2$`Stomachs 2004-2006`), big.mark = ","),
                       `Stomachs 2018-2019` = format(sum(s2$`Stomachs 2018-2019`), big.mark = ","),
                       `Sets 2004-2006` = "", `Sets 2018-2019` = ""))
add_table("Table S2",
  paste("Predator × size-class cells classifiable in the full between-period contrast (2004-2006 versus",
        "2018-2019), with the numbers of non-empty stomachs and trawl sets in each sampling window. Stomach",
        "totals are summed over cells; trawl sets are shared among predators, so they are not summed."),
  s2)

cap_unit <- function(cur) paste0(
  "Full between-period (2004-2006 versus 2018-2019) family frequencies (%) in ", cur,
  " by stratum, the values shown in the bottom panels of Figure 6. Values are means across 100 prey",
  " resolutions; n is the number of classifiable predator \u00d7 size-class cells; * marks strata with fewer",
  " than five cells.")
strat_tab <- function(m) unit_rows("stratum", m) %>% mutate(Unit = ifelse(n < 5, paste0(Unit, "*"), Unit))
add_table("Table S3", cap_unit("biomass"),    strat_tab("biomass"))
add_table("Table S4", cap_unit("occurrence"), strat_tab("occurrence"))

# Rebuild Table A5 from the saved primary runs when its appendix CSV is absent.
a5_file <- "Output_Appendix_rebuilt/TableA5_resolution_range.csv"
if (!file.exists(a5_file)) {
  message("Rebuilding missing appendix table: ", a5_file)
  primary <- load_all_runs("all_gulf", "_1", tag = "_T", verbose = FALSE)
  unit_col <- attr(primary, "unit_col")
  if (is.null(unit_col) || is.na(unit_col) || !unit_col %in% names(primary)) {
    stop("Cannot identify the spatial-unit column in saved primary runs.")
  }
  a5 <- primary %>%
    filter(testable, !is.na(family)) %>%
    mutate(family = sub("Reorganisation", "Reorganization", as.character(family), fixed = TRUE)) %>%
    count(contrast, contrast_type, mode, x_threshold,
          !!rlang::sym(unit_col), family, name = "n") %>%
    group_by(contrast, contrast_type, mode, x_threshold, !!rlang::sym(unit_col)) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup() %>%
    tidyr::complete(
      tidyr::nesting(contrast, contrast_type, mode, x_threshold,
                     !!rlang::sym(unit_col)),
      family = FAM, fill = list(n = 0L, pct = 0)) %>%
    group_by(contrast, contrast_type, mode, x_threshold, family) %>%
    summarise(pct = mean(pct), .groups = "drop") %>%
    group_by(contrast, contrast_type, mode, family) %>%
    summarise(Mean = mean(pct), Minimum = min(pct), Maximum = max(pct),
              SD = sd(pct), .groups = "drop")
  contrast_labels <- c(
    P1a = "2004 vs 2006", P1b = "2004–2005 vs 2006",
    P2 = "2018 vs 2019", PTb = "2006 vs 2018",
    PTa = "2004–2006 vs 2018–2019")
  mode_labels <- c(biomass = "Biomass", occurrence = "Occurrence")
  a5_out <- a5 %>%
    mutate(
      Contrast = unname(contrast_labels[as.character(contrast)]),
      Currency = unname(mode_labels[as.character(mode)]),
      cell = sprintf("%s [%s–%s]", f1(Mean), f1(Minimum), f1(Maximum))) %>%
    select(Contrast, Currency, family, cell) %>%
    pivot_wider(names_from = family, values_from = cell) %>%
    select(Contrast, Currency, all_of(FAM)) %>%
    arrange(Currency, Contrast)
  dir.create(dirname(a5_file), recursive = TRUE, showWarnings = FALSE)
  write_csv(a5_out, a5_file)
  message("Saved: ", a5_file)
}

add_table("Table S5",
  paste("sGSL-wide family frequencies across 100 primary resolutions. Entries show mean [minimum-maximum].",
        "These ranges describe resolution sensitivity, not confidence intervals."),
  rd("Output_Appendix_rebuilt/TableA5_resolution_range.csv") %>%
    mutate(Contrast = CTL[norm_ct(Contrast)],
           across(all_of(FAM), ~ gsub("–", "-", .x, fixed = TRUE))))

RULES <- c("_PP" = "Per-predator", "_1" = "Pooled (primary)", "_2" = "Pooled cross-predator")
grp <- bind_rows(lapply(LV, function(lv) bind_rows(lapply(names(RULES), function(pf) {
  r <- suppressWarnings(read_all_runs(levels = lv, dirs = setNames(paste0("Sensitivity_", lv, pf, if (pf == "_1") "_T" else ""), lv)))
  keys <- if (lv == "all_gulf") c("currency", "contrast") else c("currency", "contrast", "spatial_unit")
  family_by_resolution(r, keys = keys) %>%
    group_by(currency, contrast, x_threshold, family) %>% summarise(pct = mean(pct), .groups = "drop") %>%
    group_by(currency, contrast, family) %>% summarise(pct_mean = mean(pct), .groups = "drop") %>%
    mutate(level = lv, Rule = RULES[[pf]])
})))) %>% clean_fam()
W3 <- c("P1a", "P1b", "P2")
grp_tab <- function(lv) grp %>% filter(level == lv) %>%
  group_by(Currency = ML[currency], Rule) %>%
  summarise(
    `Stability within`    = paste0(f1(min(pct_mean[family == "Stability"    & contrast %in% W3])), "-",
                                   f1(max(pct_mean[family == "Stability"    & contrast %in% W3]))),
    `Substitution within` = paste0(f1(min(pct_mean[family == "Substitution" & contrast %in% W3])), "-",
                                   f1(max(pct_mean[family == "Substitution" & contrast %in% W3]))),
    `Stability between`    = f1(pct_mean[family == "Stability"    & contrast == "PTa"]),
    `Substitution between` = f1(pct_mean[family == "Substitution" & contrast == "PTa"]),
    `Leading family` = family[contrast == "PTa"][which.max(pct_mean[contrast == "PTa"])],
    .groups = "drop") %>%
  mutate(Rule = factor(Rule, unname(RULES))) %>% arrange(Currency, Rule) %>% mutate(Rule = as.character(Rule))
cap_grp <- function(scale) paste0(
  "Comparison of the three grouping rules at the ", scale, " scale, across all resolutions of each rule",
  " (100 thresholds for the pooled rules, 50 for the per-predator rule). Within-period entries are ranges",
  " of the three within-period contrast means; between-period entries refer to the full contrast.")
add_table("Table S6",  cap_grp("sGSL-wide"), grp_tab("all_gulf"))
add_table("Table S7",  cap_grp("ecoregion"), grp_tab("ecoregion"))
add_table("Table S8", cap_grp("stratum"),   grp_tab("stratum"))

add_table("Table S9",
  paste("Classification stability across prey resolutions, by scale, contrast and currency; 2006 versus",
        "2018 is not shown at the ecoregion and stratum scales."),
  cells %>% filter(level == "all_gulf" | contrast != "PTb") %>%
    group_by(level, contrast, mode) %>%
    summarise(Cells = n(), `Median support` = f1(median(support), 2),
              `Support ≥ 0.80 (%)` = f1(100 * mean(support >= 0.8), 1),
              `Mean families per cell` = f1(mean(n_families), 2), .groups = "drop") %>%
    mutate(level = factor(level, LV), contrast = factor(contrast, CT)) %>%
    arrange(level, mode, contrast) %>%
    transmute(Scale = LL[as.character(level)], Contrast = CTL[as.character(contrast)], Currency = ML[mode],
              Cells, `Median support`, `Support ≥ 0.80 (%)`, `Mean families per cell`))

cal_dir <- if (file.exists("Calibrated_all_gulf_1_T/calibrated_vs_nominal.csv") &&
               file.exists("Calibrated_all_gulf_1_T/calibrated_alphas.csv")) {
  "Calibrated_all_gulf_1_T"
} else {
  "Calibrated_all_gulf_1"
}
bh  <- rr("all_gulf", "M2_FDR_family_frequencies.csv") %>% filter(contrast == "PTa")
cal <- rd(file.path(cal_dir, "calibrated_vs_nominal.csv")) %>% filter(contrast == "PTa")
add_table("Table S10",
  paste("sGSL-wide full between-period family frequencies (%) under nominal, BH-adjusted and null-calibrated",
        "thresholds. Values are means across 100 resolutions; the two adjustments were evaluated separately."),
  bind_rows(transmute(bh,  mode, family, Thresholds = "Nominal",     pct = pct_raw),
            transmute(bh,  mode, family, Thresholds = "BH adjusted", pct = pct_bh),
            transmute(cal, mode, family, Thresholds = "Calibrated",  pct = pct_calibrated)) %>%
    mutate(Thresholds = factor(Thresholds, c("Nominal", "BH adjusted", "Calibrated")),
           family = factor(family, FAM), pct = f1(pct)) %>%
    arrange(mode, Thresholds) %>%
    pivot_wider(names_from = family, values_from = pct) %>%
    transmute(Currency = ML[mode], Thresholds = as.character(Thresholds), !!!syms(FAM)))

fpr <- rd("NullCalibration_all_gulf_1/null_false_positive_rates.csv")
alp <- rd(file.path(cal_dir, "calibrated_alphas.csv"))
add_table("Table S11",
  paste("Null rejection rates at α = 0.05 and calibrated significance thresholds. Period labels were",
        "randomized within cells: 200 randomizations at each of ten prey resolutions."),
  fpr %>% left_join(alp, by = "mode") %>%
    transmute(Currency = ML[mode], `C rate (%)` = f1(fpr_comp, 2), `H′ rate (%)` = f1(fpr_H, 2),
              `Bs rate (%)` = f1(fpr_Bs, 2), `α C` = f1(alpha_C, 4), `α H′` = f1(alpha_H, 4),
              `α Bs` = f1(alpha_Bs, 4)))

add_table("Table S12",
  paste("sGSL-wide family frequencies (%) with equalized sampling effort. Entries show the mean and, in",
        "parentheses, the standard deviation across ten resampling replicates, each summarized over ten prey",
        "resolutions (50-950 in increments of 100)."),
  rd("Robustness_all_gulf_1_BAL/balanced_family_frequencies.csv") %>%
    filter(family %in% c("Stability", "Substitution")) %>%
    mutate(contrast = factor(contrast, CT), value = msd(mean_pct, sd_pct)) %>%
    arrange(mode, contrast) %>%
    transmute(Currency = ML[mode], Contrast = CTL[as.character(contrast)], family, value) %>%
    pivot_wider(names_from = family, values_from = value))

read_resample_design <- function(level) {
  name <- sprintf("resample_design_%s_1.csv", level)
  csv_paths <- c(name, file.path("data/Sensitivity", name),
                 file.path("Sensitivity", name))

  for (path in csv_paths) {
    if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0) next
    design <- read_csv(path, show_col_types = FALSE)
    if (!is.data.frame(design) || !"keep" %in% names(design) || !nrow(design))
      stop("Invalid saved sampling design: ", path)
    design$keep <- as.logical(design$keep)
    if (anyNA(design$keep)) stop("Invalid retention flags in ", path)
    return(design)
  }

  checkpoints <- unlist(lapply(c("data/Sensitivity", "Sensitivity"), function(folder)
    file.path(folder, sprintf("all_runs_%s_1_%s.rda", level,
                              c("BALCONS", "BAL")))))
  for (path in checkpoints) {
    if (!file.exists(path)) next
    e <- new.env(parent = baseenv())
    objects <- tryCatch(load(path, envir = e), error = function(err) {
      warning("Cannot read ", path, ": ", conditionMessage(err), call. = FALSE)
      character()
    })
    for (object_name in objects) {
      object <- e[[object_name]]
      design <- if (is.list(object) && is.data.frame(object$inventory)) {
        object$inventory
      } else if (is.data.frame(object)) {
        attr(object, "design")
      } else {
        NULL
      }
      if (!is.data.frame(design) || !"keep" %in% names(design) || !nrow(design)) next
      design$keep <- as.logical(design$keep)
      if (anyNA(design$keep)) stop("Invalid retention flags in ", path)
      write_csv(design, name)
      message("Restored sampling design from ", path, " -> ", name)
      return(design)
    }
  }
  stop("Missing sampling design for ", level,
       ". Neither resample_design_*_1.csv nor saved BAL/BALCONS checkpoints contain the design.")
}

add_table("Table S13",
  paste("Number of inventory cells retained and excluded when equalizing sampling effort, by spatial scale.",
        "These counts are not the classifiable-cell counts used in the main text."),
  map_dfr(LV, function(l) {
    d <- read_resample_design(l)
    tibble(Scale = LL[[l]], `Inventory cells` = nrow(d), Retained = sum(d$keep), Excluded = sum(!d$keep))
  }))

add_table("Table S14",
  paste("Separation in percentage points after restricting contrasts to the common range of period-specific",
        "set counts. Both between-period contrasts contribute. Positive values indicate separation."),
  map_dfr(LV, function(l) rd(sprintf("Effort_%s_1/C_separation_verdict.csv", l)) %>% mutate(level = l)) %>%
    transmute(Scale = LL[level], Currency = ML[mode], Family = family,
              `Unrestricted gap` = f1(gap_all), `Common-effort gap` = f1(gap_common)))

FIG_DIR <- "Output_Paper_Figures"
FIG_CAP <- c(
  S1 = "Trawl sets from which predator stomachs were analysed, by survey year and pooled over years, in the four ecoregions. Labels give the number of non-empty stomachs analysed in each ecoregion; panel titles give the yearly totals of stomachs and trawl sets. Totals cover all retained stomachs and therefore differ slightly from Table S2, which counts only cells classifiable in the full between-period contrast.",
  S2 = "Prey identification depth by period: percentage of ranked prey records identified at each taxonomic rank or finer, in 2004-2006 and 2018-2019.",
  S3 = "Ecoregion-scale family frequencies across five contrasts. Points show the mean family frequency across classifiable ecoregions (all four ecoregions contribute to every contrast) and bars show the standard deviation of that mean across 100 prey resolutions.",
  S4 = "Family frequencies within each ecoregion across five contrasts (n = largest number of classifiable predator \u00d7 size-class cells across contrasts and currencies). Points show means and bars show standard deviations across 100 prey resolutions (solid lines, biomass; dashed lines, occurrence). Shading as in Figure 4.",
  S5 = "Family frequencies for 2004 vs 2006 by ecoregion (top) and stratum (bottom), in biomass and occurrence. Layout and colours as in Figure 5. Values are means across 100 resolutions; parentheses count classifiable predator \u00d7 size-class cells and * marks units with fewer than five. Grey rows have no classifiable cell.",
  S6 = "Family frequencies for 2018 vs 2019 by ecoregion (top) and stratum (bottom), in biomass and occurrence. Layout and colours as in Figure 5. Values are means across 100 resolutions; parentheses count classifiable predator \u00d7 size-class cells and * marks units with fewer than five. Grey rows have no classifiable cell.",
  S7 = "Family frequencies for 2006 vs 2018 by ecoregion (top) and stratum (bottom), in biomass and occurrence. Layout and colours as in Figure 5. Values are means across 100 resolutions; parentheses count classifiable predator \u00d7 size-class cells and * marks units with fewer than five. Grey rows have no classifiable cell.",
  S8 = "Stratum-scale family frequencies across five contrasts. Points show the mean family frequency across classifiable strata and bars show the standard deviation of that mean across 100 prey resolutions.",
  S9 = "Family frequencies mapped on the survey strata, biomass, for four of the five contrasts (2006 versus 2018 is not mapped). Values are within-stratum frequencies (%) averaged across 100 prey resolutions; * marks strata built on fewer than five classifiable predator \u00d7 size-class cells.",
  S10 = "As Figure S9, occurrence.",
  S11 = "Full between-period (2004-2006 versus 2018-2019) family frequencies at three spatial scales under the primary grouping rule. Points show means and bars show standard deviations across 100 prey resolutions; ecoregion and stratum values are means across spatial units.",
  S12 = "sGSL-wide frequencies across 100 prey aggregation thresholds (10-1,000 stomachs, increments of 10). Families absent at a threshold are shown at zero.",
  S13 = "As Figure S12, ecoregion scale, without 2006 versus 2018. At each threshold, family frequencies are computed within each ecoregion and averaged across classifiable ecoregions.",
  S14 = "As Figure S12, stratum scale, without 2006 versus 2018. At each threshold, family frequencies are computed within each stratum and averaged across classifiable strata.",
  S15 = "Cumulative distribution of modal-family support in the full between-period contrast. Each curve shows, for a given scale and currency, the fraction of cells (y) with modal-family support at or below the value on the x axis; a curve that stays low until 0.80 and then rises steeply indicates that most cells keep the same family across resolutions. The vertical line marks support of 0.80.",
  S16 = "Full between-period family frequencies under alternative significance thresholds. Thresholds shown are \u03b1 = 0.01, 0.025, 0.05 and 0.075; the upper bound for flagging a composition trend is fixed at 0.10.",
  S17 = "Bootstrap 95% intervals for Stability and Substitution, by scale and currency (2,000 cell-bootstrap replicates). Points show the mean and bars the 95% interval of each contrast; shaded boxes show the within-period envelope (three contrasts) and the between-period envelope (two contrasts). The gap is the distance between the two envelopes in percentage points; positive values indicate separation, negative values overlap.",
  S18 = "Null rejection rates by the smaller period-specific set count. The dashed line marks the nominal 5% rate.",
  S19 = "Stability and Substitution under equalized effort, sGSL-wide. Points are the ten resampling replicates; boxes show medians and interquartile ranges. Whiskers extend to observations within 1.5 interquartile ranges. Each replicate is summarized over ten prey resolutions. Crosses show the frequencies without effort balancing.",
  S20 = "As Figure S19, ecoregion scale. Equalizing effort retained 68 of 200 inventory cells at this scale (Table S13), so these frequencies describe a small subset of the ecoregion cells."
)
FIG_SECTION <- c(S1 = "Data", S3 = "Results", S12 = "Sensitivity", S17 = "Robustness")

new_doc <- function() read_docx() %>%
  body_set_default_section(prop_section(
    page_size = page_size(orient = "portrait", width = 8.5, height = 11),
    page_margins = page_mar(top = 0.8, bottom = 0.8, left = 1, right = 1)))
TNR <- function(size, bold = FALSE) fp_text(bold = bold, font.family = "Times New Roman", font.size = size)
cap_par <- function(label, caption) fpar(ftext(paste0(label, ". "), TNR(10, TRUE)), ftext(caption, TNR(10)),
                                         fp_p = fp_par(padding.bottom = 4))
heading <- function(doc, txt) body_add_fpar(doc, fpar(ftext(txt, TNR(13, TRUE))))
marker  <- function(doc, txt) body_add_fpar(doc, fpar(ftext(txt, TNR(11, TRUE))))
add_tab <- function(doc, t) {
  doc <- body_add_fpar(doc, cap_par(t$label, t$caption))
  body_add_flextable(doc, ft_style(t$data))
}
add_fig <- function(doc, id) {
  f <- file.path(FIG_DIR, paste0("Figure_", id, ".png"))
  if (!file.exists(f)) { warning("missing ", f); return(doc) }
  info <- magick::image_info(magick::image_read(f))
  w <- 6.5; h <- min(w * info$height / info$width, 8.6); w <- h * info$width / info$height
  doc <- body_add_img(doc, f, width = w, height = h)
  body_add_fpar(doc, cap_par(paste("Figure", id), FIG_CAP[[id]]))
}
is_main <- function(t) !grepl("^Table S", t$label)

doc <- heading(new_doc(), "Tables")
main <- Filter(is_main, TABLES)
for (k in seq_along(main)) {
  doc <- add_tab(doc, main[[k]])
  if (k < length(main)) doc <- body_add_break(doc)
}
print(doc, target = file.path(OUT, "Main_Tables.docx"))

SECTION <- c("Table S1" = "Data", "Table S3" = "Results", "Table S5" = "Sensitivity", "Table S10" = "Robustness")
doc <- heading(new_doc(), "Supplementary figures")
ids <- paste0("S", 1:20)
for (k in seq_along(ids)) {
  if (ids[k] %in% names(FIG_SECTION)) doc <- marker(doc, FIG_SECTION[[ids[k]]])
  doc <- add_fig(doc, ids[k])
  doc <- body_add_break(doc)
}
doc <- heading(doc, "Supplementary tables")
supp <- Filter(Negate(is_main), TABLES)
for (k in seq_along(supp)) {
  t <- supp[[k]]
  if (t$label %in% names(SECTION)) doc <- marker(doc, SECTION[[t$label]])
  doc <- add_tab(doc, t)
  if (k < length(supp)) doc <- body_add_break(doc)
}
print(doc, target = file.path(OUT, "Supplementary_materials.docx"))
cat("\nWrote Main_Tables.docx (", length(main), " tables) and Supplementary_materials.docx (",
    length(ids), " figures, ", length(supp), " tables) in ", OUT, "\n", sep = "")

})
