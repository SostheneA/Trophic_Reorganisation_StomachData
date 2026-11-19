# Post-hoc tests and classification robustness
# Input: data/Sensitivity/all_runs_<level>_1_T.rda
# Output: Robustness_<level>_1_T/

SPATIAL_LEVELS <- c("all_gulf", "ecoregion", "stratum")

# Input files
input_files <- file.path(
  "data", "Sensitivity",
  paste0("all_runs_", SPATIAL_LEVELS, "_1_T.rda")
)
missing <- input_files[!file.exists(input_files)]
if (length(missing)) {
  stop("Missing saved results: ", paste(missing, collapse = ", "))
}
if (!file.exists("R_helpers/Results.R")) {
  stop("Missing helper: R_helpers/Results.R")
}

for (level in SPATIAL_LEVELS) {
  message("Post-hoc robustness: ", level)
  local({
    SPATIAL_LEVEL <- level
    PREY_FAMILY <- "_1"
    RUN_TAG <- "_T"
    RUN_BASE <- "Sensitivity"

    # Effect sizes, multiple testing and bootstrap intervals
    suppressPackageStartupMessages({
      library(dplyr); library(tidyr); library(purrr); library(readr)
      library(ggplot2); library(tibble)
    })

    if (!exists("PREY_FAMILY"))   PREY_FAMILY   <- "_1"
    if (!exists("SPATIAL_LEVEL")) SPATIAL_LEVEL <- "all_gulf"
    if (!exists("RUN_TAG"))       RUN_TAG       <- ""
    if (!exists("RUN_BASE"))      RUN_BASE      <- "Sensitivity"

    ALPHA <- 0.05; TREND_UPPER <- 0.10
    B_BOOT <- 2000; SEED <- 20260911

    MARGIN_H  <- 0.20
    MARGIN_BS <- 0.05

    OUT <- sprintf("Robustness_%s%s%s%s", SPATIAL_LEVEL, PREY_FAMILY, RUN_TAG,
                   if (RUN_BASE == "Sensitivity") "" else "_T")
    dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

    source("R_helpers/Results.R", local = TRUE)

    runs <- load_all_runs(SPATIAL_LEVEL, PREY_FAMILY, tag = RUN_TAG, base = RUN_BASE)
    UNIT <- attr(runs, "unit_col")

    runs <- runs %>% filter(testable)
    cls  <- runs %>% filter(!is.na(family))

    reclassify <- function(d, alpha = ALPHA, trend = TREND_UPPER,
                           pc = "p_comp", ph = "p_H", pb = "p_Bs") {
      H  <- !is.na(d[[ph]]) & d[[ph]] < alpha
      B  <- !is.na(d[[pb]]) & d[[pb]] < alpha
      Cs <- !is.na(d[[pc]]) & d[[pc]] < alpha
      Ct <- !is.na(d[[pc]]) & d[[pc]] < trend
      dg <- dplyr::case_when(
        !H & !B & !Cs & !Ct ~ "Stable Diet",
        !H & !B & !Cs &  Ct ~ "Emerging Shift",
        !H & !B &  Cs       ~ "Ghost Shift",
        H & !B & !Cs       ~ "Internal Rebalancing",
        !H &  B & !Cs       ~ "Niche Compression/Expansion",
        H &  B & !Cs       ~ "Niche Restructuring",
        H & !B &  Cs       ~ "Partial Diet Shift",
        !H &  B &  Cs       ~ "Structural Shift",
        H &  B &  Cs       ~ "Major Shift",
        .default = INCONCLUSIVE_LAB
      )
      ok <- !is.na(d[[pc]]) & !is.na(d[[ph]]) & !is.na(d[[pb]])
      dg[!ok] <- INCONCLUSIVE_LAB
      d$diagnostic_alt <- dg
      d$family_alt <- as.character(family_of(dg))
      d
    }

    pct_fam <- function(d, fam = "family") {
      step1 <- d %>%
        filter(!is.na(.data[[fam]])) %>%
        count(contrast, contrast_type, mode, x_threshold, !!rlang::sym(UNIT),
              family = .data[[fam]], name = "n") %>%
        group_by(contrast, contrast_type, mode, x_threshold, !!rlang::sym(UNIT)) %>%
        mutate(pct = 100 * n / sum(n)) %>%
        ungroup() %>%
        tidyr::complete(
          tidyr::nesting(contrast, contrast_type, mode, x_threshold,
                         !!rlang::sym(UNIT)),
          family = FAMILY_LEVELS, fill = list(n = 0L, pct = 0))
      step2 <- step1 %>%
        group_by(contrast, contrast_type, mode, x_threshold, family) %>%
        summarise(pct = mean(pct), .groups = "drop")
      step2 %>%
        group_by(contrast, contrast_type, mode, family) %>%
        summarise(sd_across_x = sd(pct), pct_mean = mean(pct), .groups = "drop") %>%
        select(contrast, contrast_type, mode, family, pct_mean, sd_across_x)
    }

    ref <- pct_fam(cls)
    write_csv(ref, file.path(OUT, "M0_family_frequencies_reference.csv"))

    chk <- ref %>% group_by(contrast, mode) %>%
      summarise(total = sum(pct_mean), .groups = "drop") %>%
      filter(abs(total - 100) > 0.01)
    if (nrow(chk)) { print(as.data.frame(chk))
      stop("Family frequencies do not sum to 100%.") }
    message("Family frequency totals verified.")
    print(as.data.frame(ref %>% filter(family %in% c("Stability", "Substitution"))))

    cls <- cls %>%
      mutate(se_delta_H = if ("t_H" %in% names(.)) abs(delta_H) / pmax(abs(t_H), 1e-9)
             else NA_real_,
             se_delta_Bs = if ("t_Bs" %in% names(.)) abs(delta_Bs) / pmax(abs(t_Bs), 1e-9)
             else NA_real_)

    HAS_BS_SE <- any(is.finite(cls$se_delta_Bs))
    if (!HAS_BS_SE) {
      message("Bs standard errors are unavailable; TOST is limited to H'.")
    }

    eff_summary <- cls %>%
      group_by(contrast, contrast_type, mode, family) %>%
      summarise(n = n(),
                dH_med  = median(abs(delta_H),  na.rm = TRUE),
                dH_q90  = quantile(abs(delta_H),  .90, na.rm = TRUE),
                dBs_med = median(abs(delta_Bs), na.rm = TRUE),
                dBs_q90 = quantile(abs(delta_Bs), .90, na.rm = TRUE),
                BC_med  = median(BC, na.rm = TRUE),
                R2_med  = median(R2_comp, na.rm = TRUE), .groups = "drop")
    write_csv(eff_summary, file.path(OUT, "M1_effect_sizes_by_family.csv"))

    tost <- cls %>%
      filter(family == "Substitution") %>%
      mutate(lo_H = delta_H - 1.645 * se_delta_H,
             hi_H = delta_H + 1.645 * se_delta_H,
             equiv_H = lo_H > -MARGIN_H & hi_H < MARGIN_H,
             lo_Bs = if (HAS_BS_SE) delta_Bs - 1.645 * se_delta_Bs else NA_real_,
             hi_Bs = if (HAS_BS_SE) delta_Bs + 1.645 * se_delta_Bs else NA_real_,
             equiv_Bs = if (HAS_BS_SE) lo_Bs > -MARGIN_BS & hi_Bs < MARGIN_BS else NA,
             equiv_both = if (HAS_BS_SE) equiv_H & equiv_Bs else NA)

    tost_summary <- tost %>%
      group_by(contrast, contrast_type, mode) %>%
      summarise(n_substitution = n(),
                pct_equiv_H  = 100 * mean(equiv_H,  na.rm = TRUE),
                pct_equiv_Bs = if (HAS_BS_SE) 100 * mean(equiv_Bs, na.rm = TRUE) else NA_real_,
                pct_equiv_both = if (HAS_BS_SE) 100 * mean(equiv_both, na.rm = TRUE) else NA_real_,
                .groups = "drop")
    write_csv(tost_summary, file.path(OUT, "M1_TOST_substitution.csv"))
    print(as.data.frame(tost_summary))

    tost_curve <- tibble(margin_H = seq(0.05, 0.60, by = 0.05)) %>%
      mutate(pct_equiv = map_dbl(margin_H, ~ 100 * mean(
        tost$lo_H > -.x & tost$hi_H < .x, na.rm = TRUE)))
    write_csv(tost_curve, file.path(OUT, "M1b_TOST_margin_sensitivity.csv"))

    p_eff <- cls %>%
      filter(contrast_type == "between") %>%
      transmute(mode, family = factor(family, levels = FAMILY_LEVELS),
                `|delta H'| (nats)` = abs(delta_H), `|delta Bs|` = abs(delta_Bs)) %>%
      pivot_longer(c(`|delta H'| (nats)`, `|delta Bs|`),
                   names_to = "index", values_to = "value") %>%
      ggplot(aes(family, value, fill = family)) +
      geom_violin(scale = "width", alpha = .6, colour = NA) +
      geom_boxplot(width = .13, outlier.size = .3) +
      facet_grid(index ~ mode, scales = "free_y") +
      labs(x = NULL, y = "Absolute effect size",
           title = "Univariate effect sizes by family (between-period contrasts)") +
      theme_bw(base_size = 10) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
    ggsave(file.path(OUT, "M1_effect_sizes.png"), p_eff, width = 9, height = 6, dpi = 300)

    runs_fdr <- runs %>%
      group_by(contrast, mode, x_threshold) %>%
      mutate(q_comp = p.adjust(p_comp, "BH"),
             q_H    = p.adjust(p_H,    "BH"),
             q_Bs   = p.adjust(p_Bs,   "BH")) %>%
      ungroup() %>%
      reclassify(pc = "q_comp", ph = "q_H", pb = "q_Bs")

    fdr_comparison <- ref %>% rename(pct_raw = pct_mean) %>% select(-sd_across_x) %>%
      left_join(pct_fam(runs_fdr, "family_alt") %>%
                  select(contrast, contrast_type, mode, family, pct_bh = pct_mean),
                by = c("contrast", "contrast_type", "mode", "family")) %>%
      mutate(pct_bh = coalesce(pct_bh, 0), shift_pp = pct_bh - pct_raw)
    write_csv(fdr_comparison, file.path(OUT, "M2_FDR_family_frequencies.csv"))

    set.seed(SEED)
    cell_list <- sort(unique(cls$cell)); K <- length(cell_list)

    base_counts <- cls %>%
      count(contrast, contrast_type, mode, x_threshold, !!rlang::sym(UNIT),
            cell, family, name = "n_cell")

    boot_one <- function(b) {
      w <- tibble(cell = cell_list, w = as.integer(rmultinom(1, K, rep(1 / K, K))))
      base_counts %>%
        inner_join(w, by = "cell") %>% filter(w > 0) %>%
        count(contrast, contrast_type, mode, x_threshold, !!rlang::sym(UNIT),
              family, wt = n_cell * w, name = "n") %>%
        group_by(contrast, contrast_type, mode, x_threshold, !!rlang::sym(UNIT)) %>%
        mutate(pct = 100 * n / sum(n)) %>%
        ungroup() %>%
        tidyr::complete(
          tidyr::nesting(contrast, contrast_type, mode, x_threshold,
                         !!rlang::sym(UNIT)),
          family = FAMILY_LEVELS, fill = list(n = 0L, pct = 0)) %>%
        group_by(contrast, contrast_type, mode, x_threshold, family) %>%
        summarise(pct = mean(pct), .groups = "drop") %>%
        group_by(contrast, contrast_type, mode, family) %>%
        summarise(pct = mean(pct), .groups = "drop")
    }

    boot <- map_dfr(seq_len(B_BOOT), boot_one)

    boot_ci <- boot %>%
      group_by(contrast, contrast_type, mode, family) %>%
      summarise(pct_mean = mean(pct), ci_lo = quantile(pct, .025),
                ci_hi = quantile(pct, .975), .groups = "drop")
    write_csv(boot_ci, file.path(OUT, "M3_bootstrap_CI_family_pct.csv"))

    overlap <- boot_ci %>%
      filter(family %in% c("Stability", "Substitution")) %>%
      group_by(mode, family, contrast_type) %>%
      summarise(lo = min(ci_lo), hi = max(ci_hi), .groups = "drop") %>%
      pivot_wider(names_from = contrast_type, values_from = c(lo, hi)) %>%
      mutate(disjoint = if_else(family == "Stability",
                                hi_between < lo_within, lo_between > hi_within),
             gap_pp = if_else(family == "Stability",
                              lo_within - hi_between, lo_between - hi_within))
    write_csv(overlap, file.path(OUT, "M3_CI_overlap_verdict.csv"))
    print(as.data.frame(overlap))

    alpha_sens <- expand_grid(a = c(.01, .025, .05, .075, .10),
                              tr = c(.05, .10, .15, .20)) %>%
      filter(tr > a) %>%
      pmap_dfr(function(a, tr)
        pct_fam(reclassify(runs, alpha = a, trend = tr), "family_alt") %>%
          mutate(alpha = a, trend = tr))
    write_csv(alpha_sens, file.path(OUT, "M4_alpha_trend_sensitivity.csv"))

    p_alpha <- alpha_sens %>% filter(trend == .10) %>%
      mutate(family = factor(family, levels = FAMILY_LEVELS)) %>%
      ggplot(aes(alpha, pct_mean, colour = contrast_type, group = contrast)) +
      geom_line(alpha = .8) + geom_point(size = 1) +
      facet_grid(mode ~ family) +
      labs(x = "alpha", y = "Classifiable cells (%)", colour = NULL,
           title = "Family frequencies across significance thresholds") +
      theme_bw(base_size = 10)
    ggsave(file.path(OUT, "M4_alpha_sensitivity.png"), p_alpha,
           width = 10, height = 5.5, dpi = 300)

    stability_cell <- cls %>%
      group_by(contrast, mode, cell) %>%
      summarise(n_thresholds = n(),
                modal_family = names(which.max(table(family))),
                support = max(table(family)) / n(),
                n_families = n_distinct(family), .groups = "drop")
    write_csv(stability_cell, file.path(OUT, "M5_cell_classification_stability.csv"))

    stability_summary <- stability_cell %>%
      group_by(contrast, mode) %>%
      summarise(pct_support_ge_80 = 100 * mean(support >= .80),
                pct_support_ge_60 = 100 * mean(support >= .60),
                median_support = median(support), .groups = "drop")
    write_csv(stability_summary, file.path(OUT, "M5_stability_summary.csv"))
    print(as.data.frame(stability_summary))

    set.seed(SEED + 77)

    agg_boot <- function(d, B = 2000) {
      cells <- sort(unique(d$cell)); K <- length(cells)
      if (K < 3) return(NULL)
      per_cell <- d %>%
        group_by(cell) %>%
        summarise(dH = mean(delta_H, na.rm = TRUE),
                  adH = mean(abs(delta_H), na.rm = TRUE),
                  dBs = mean(delta_Bs, na.rm = TRUE),
                  adBs = mean(abs(delta_Bs), na.rm = TRUE), .groups = "drop")
      draws <- replicate(B, {
        i <- sample.int(K, K, replace = TRUE)
        c(mean(per_cell$dH[i]), mean(per_cell$adH[i]),
          mean(per_cell$dBs[i]), mean(per_cell$adBs[i]))
      })
      tibble(
        n_cells   = K,
        mean_dH   = mean(per_cell$dH),
        dH_lo     = quantile(draws[1, ], .025), dH_hi = quantile(draws[1, ], .975),
        mean_absdH = mean(per_cell$adH),
        absdH_lo  = quantile(draws[2, ], .025), absdH_hi = quantile(draws[2, ], .975),
        mean_dBs  = mean(per_cell$dBs),
        dBs_lo    = quantile(draws[3, ], .025), dBs_hi = quantile(draws[3, ], .975),
        mean_absdBs = mean(per_cell$adBs),
        absdBs_lo = quantile(draws[4, ], .025), absdBs_hi = quantile(draws[4, ], .975))
    }

    agg_equiv <- cls %>%
      group_by(contrast, contrast_type, mode, family) %>%
      group_modify(~ { r <- agg_boot(.x, B_BOOT); if (is.null(r)) tibble() else r }) %>%
      ungroup() %>%
      mutate(equiv_H_at_margin  = dH_lo  > -MARGIN_H  & dH_hi  < MARGIN_H,
             equiv_Bs_at_margin = dBs_lo > -MARGIN_BS & dBs_hi < MARGIN_BS)
    write_csv(agg_equiv, file.path(OUT, "M6_aggregate_equivalence.csv"))
    cat("\n--- M6: aggregate equivalence for Substitution cells ---\n")
    print(as.data.frame(agg_equiv %>% filter(family == "Substitution") %>%
                          select(contrast, mode, n_cells, mean_dH, dH_lo, dH_hi, equiv_H_at_margin,
                                 mean_dBs, dBs_lo, dBs_hi, equiv_Bs_at_margin)))

    contrast_pairs <- cls %>%
      filter(family %in% c("Substitution", "Reorganisation")) %>%
      group_by(contrast, contrast_type, mode) %>%
      group_modify(function(d, k) {
        per_cell <- d %>% group_by(cell, family) %>%
          summarise(adH = mean(abs(delta_H), na.rm = TRUE), .groups = "drop")
        sub <- per_cell$adH[per_cell$family == "Substitution"]
        reo <- per_cell$adH[per_cell$family == "Reorganisation"]
        if (length(sub) < 3 || length(reo) < 3) return(tibble())
        dr <- replicate(B_BOOT,
                        mean(sample(reo, length(reo), TRUE)) - mean(sample(sub, length(sub), TRUE)))
        tibble(n_sub = length(sub), n_reo = length(reo),
               diff_absdH = mean(reo) - mean(sub),
               lo = quantile(dr, .025), hi = quantile(dr, .975),
               ordered = quantile(dr, .025) > 0)
      }) %>% ungroup()
    write_csv(contrast_pairs, file.path(OUT, "M6b_family_ordering.csv"))
    print(as.data.frame(contrast_pairs %>% filter(contrast_type == "between")))

    # Output files
    expected <- c(
      "M0_family_frequencies_reference.csv",
      "M1_effect_sizes_by_family.csv",
      "M1_TOST_substitution.csv",
      "M1b_TOST_margin_sensitivity.csv",
      "M1_effect_sizes.png",
      "M2_FDR_family_frequencies.csv",
      "M3_bootstrap_CI_family_pct.csv",
      "M3_CI_overlap_verdict.csv",
      "M4_alpha_trend_sensitivity.csv",
      "M4_alpha_sensitivity.png",
      "M5_cell_classification_stability.csv",
      "M5_stability_summary.csv",
      "M6_aggregate_equivalence.csv",
      "M6b_family_ordering.csv"
    )
    absent <- expected[!file.exists(file.path(OUT, expected))]
    if (length(absent)) {
      stop("Missing post-hoc outputs for ", SPATIAL_LEVEL, ": ",
           paste(absent, collapse = ", "))
    }
    message("Post-hoc summaries saved: ", OUT)
  })
}
