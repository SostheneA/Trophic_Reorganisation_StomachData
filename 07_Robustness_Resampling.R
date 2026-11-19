# Balanced sampling, null simulations, calibration and effort sensitivity
# Input: data/dat_classed.rda and saved main analyses
# Output: Sensitivity_*_1_BAL/, Robustness_*_1_BAL/, Effort_*_1/,
#         NullCalibration_all_gulf_1/ and Calibrated_all_gulf_1_T/

# Balanced-run input
find_balanced_file <- function(level) {
  file_name <- sprintf("all_runs_%s_1_BAL.rda", level)
  candidates <- file.path(c("data/Sensitivity", "Sensitivity"), file_name)
  found <- candidates[vapply(candidates, function(path) {
    file.exists(path) && !is.na(file.info(path)$size) && file.info(path)$size > 0
  }, logical(1))]
  if (length(found)) found[[1L]] else NA_character_
}

run_balanced <- function(level) {
  SPATIAL_LEVEL <- level
  PREY_FAMILY <- "_1"
  RUN_TAG <- "_T"
  N_CORES <- 2L
  # Balance effort across contrasts
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(purrr); library(readr); library(tibble)
  })

  if (!exists("PREY_FAMILY"))   PREY_FAMILY   <- "_1"
  if (!exists("SPATIAL_LEVEL")) SPATIAL_LEVEL <- "all_gulf"

  N_MIN <- 5; N_STRICT <- 25; MIN_SETS <- 3; R_PERM <- 999; ALPHA <- 0.05

  B_REP <- 10
  SEED  <- 20260911

  X_GRID <- seq(50, 1000, by = 100)

  BALANCE_ALL <- TRUE

  TARGET_RULE <- "min_within"

  if (!exists("N_CORES")) N_CORES <- max(1L, parallel::detectCores())

  source("03_Statistical_Pipeline.R", local = TRUE)

  OUT_DIR <- sprintf("Sensitivity_%s%s_BAL", SPATIAL_LEVEL, PREY_FAMILY)
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  dir.create("data/Sensitivity", showWarnings = FALSE, recursive = TRUE)

  WITHIN_SC <- c("P1", "P2", "P3")
  RUN_SC    <- if (BALANCE_ALL) names(SCENARIOS) else c("P4", "PT")

  load("data/dat_classed.rda")
  diet_clean <- as.data.frame(dat_classed)
  if (!is.na(SPATIAL_SOURCE)) {
    if (!SPATIAL_SOURCE %in% names(diet_clean))
      stop("Missing column '", SPATIAL_SOURCE, "' in dat_classed.")
    diet_clean[[SPATIAL_SOURCE]] <- as.factor(diet_clean[[SPATIAL_SOURCE]])
  }

  diet_clean$set_uid <- paste(diet_clean$year, diet_clean[["vessel.code"]],
                              diet_clean$set, sep = "_")

  unit_vec <- if (isTRUE(SPATIAL_POOL)) rep(SPATIAL_LABEL, nrow(diet_clean)) else
    as.character(diet_clean[[SPATIAL_SOURCE]])
  diet_clean$.unit <- unit_vec
  diet_clean <- diet_clean[!is.na(diet_clean$.unit), , drop = FALSE]

  inv_src <- diet_clean %>%
    filter(!is.na(somatic_length_cm), somatic_length_cm > 0, !is.na(size_class)) %>%
    distinct(species = predator_species_common_name,
             size = as.character(size_class), unit = .unit,
             set_uid, year)

  n_sets_for <- function(sc) {
    pm <- SCENARIOS[[sc]]
    bind_rows(
      inv_src %>% filter(year %in% pm[[1]]) %>%
        count(species, size, unit, name = "n_sets") %>% mutate(period = 1L),
      inv_src %>% filter(year %in% pm[[2]]) %>%
        count(species, size, unit, name = "n_sets") %>% mutate(period = 2L)
    ) %>% mutate(sc = sc)
  }

  availability <- map_dfr(names(SCENARIOS), n_sets_for)

  targets <- availability %>%
    filter(sc %in% if (TARGET_RULE == "min_within") WITHIN_SC else names(SCENARIOS)) %>%
    group_by(species, size, unit) %>%
    summarise(n_target = min(n_sets), .groups = "drop")

  cap <- availability %>%
    filter(sc %in% RUN_SC) %>%
    group_by(species, size, unit) %>%
    summarise(n_min_available = min(n_sets), n_sc = n_distinct(sc), .groups = "drop")

  design <- targets %>%
    left_join(cap, by = c("species", "size", "unit")) %>%
    mutate(n_min_available = coalesce(n_min_available, 0L),
           complete_sc     = coalesce(n_sc, 0L) == length(RUN_SC),
           n_target_final  = pmin(n_target, n_min_available),
           keep            = complete_sc & n_target_final >= MIN_SETS,
           drop_reason     = case_when(
             !complete_sc                   ~ "cell absent from at least one contrast",
             n_target_final < MIN_SETS      ~ sprintf("target %d < %d trawl sets", n_target_final, MIN_SETS),
             TRUE                           ~ NA_character_))

  write_csv(design, sprintf("resample_design_%s%s.csv", SPATIAL_LEVEL, PREY_FAMILY))

  cat(sprintf("\nCells: %d total | %d retained | %d excluded\n",
              nrow(design), sum(design$keep), sum(!design$keep)))
  print(design %>% filter(!keep) %>% count(drop_reason))
  n_keep <- sum(design$keep)
  if (n_keep == 0L)
    stop("No cells remain after effort balancing at this spatial level.\n",
         "All targets are below ", MIN_SETS, " trawl sets: the ",
         "within-period contrasts have insufficient sets per cell.\n",
         "See resample_design_", SPATIAL_LEVEL, PREY_FAMILY, ".csv.")

  cat(sprintf("Trawl-set target: median %d (min %d, max %d)\n",
              median(design$n_target_final[design$keep]),
              min(design$n_target_final[design$keep]),
              max(design$n_target_final[design$keep])))

  pct_drop <- 100 * (1 - n_keep / nrow(design))
  cat(sprintf("Excluded cells: %.1f %%\n", pct_drop))
  if (pct_drop > 40)
    warning(sprintf("%.0f %% of cells were excluded. Results reflect ",
                    pct_drop),
            "a better-sampled subset and should be reported accordingly.",
            call. = FALSE)

  design_keep <- design %>% filter(keep) %>% select(species, size, unit, n_target_final)

  prey_cols_all <- grep(paste0("^prey_category_\\d+", PREY_FAMILY, "$"),
                        names(diet_clean), value = TRUE)
  keep_prey <- paste0("prey_category_", X_GRID, PREY_FAMILY)
  keep_prey <- intersect(keep_prey, prey_cols_all)
  if (!length(keep_prey))
    stop("No X_GRID prey columns found in dat_classed. Available: ",
         paste(head(prey_cols_all, 20), collapse = ", "))
  drop_prey <- setdiff(prey_cols_all, keep_prey)
  cat(sprintf("Retained resolutions: %s\n", paste(X_GRID[paste0("prey_category_", X_GRID, PREY_FAMILY) %in% keep_prey], collapse = ", ")))

  diet_grid <- diet_clean[, setdiff(names(diet_clean), drop_prey), drop = FALSE]

  draw_sets <- function(sc) {
    pm <- SCENARIOS[[sc]]
    inv_src %>%
      mutate(period = case_when(year %in% pm[[1]] ~ 1L,
                                year %in% pm[[2]] ~ 2L,
                                TRUE ~ NA_integer_)) %>%
      filter(!is.na(period)) %>%
      inner_join(design_keep, by = c("species", "size", "unit")) %>%
      group_by(species, size, unit, period) %>%
      mutate(.rank = sample.int(dplyr::n())) %>%
      filter(.rank <= n_target_final) %>%
      ungroup() %>%
      select(species, size, unit, set_uid)
  }

  TASKS <- expand.grid(rep_id = seq_len(B_REP), sc = RUN_SC,
                       mode = c("biomass", "occurrence"),
                       stringsAsFactors = FALSE)
  TASKS <- TASKS[order(TASKS$rep_id, TASKS$sc, TASKS$mode), ]
  rownames(TASKS) <- NULL

  run_one_task <- function(k) {
    b  <- TASKS$rep_id[k]
    sc <- TASKS$sc[k]
    md <- TASKS$mode[k]
    pm <- SCENARIOS[[sc]]

    set.seed(SEED + 1000L * b + match(sc, RUN_SC))
    keep <- draw_sets(sc)

    d_b <- diet_grid %>%
      inner_join(keep,
                 by = c("predator_species_common_name" = "species",
                        "size_class" = "size", ".unit" = "unit",
                        "set_uid" = "set_uid"))
    if (!nrow(d_b)) return(NULL)

    dat_b <- make_dat_classed(d_b, period_map = pm)

    r <- tryCatch(
      run_pipeline(dat_b, mode = md,
                   period_1 = names(pm)[1], period_2 = names(pm)[2])$results,
      error = function(e) {
        warning(sprintf("rep %d | %s | %s : %s", b, sc, md, conditionMessage(e)),
                call. = FALSE); NULL })
    if (is.null(r) || !nrow(r)) return(NULL)

    r$rep_id <- b; r$sc_code <- sc
    saveRDS(r, file.path(OUT_DIR, sprintf("task_%04d_%s_%s.rds", b, sc, md)))
    r
  }

  cat(sprintf("\n%d replicates x %d contrasts x 2 currencies x %d resolutions\n",
              B_REP, length(RUN_SC), length(keep_prey)))
  cat(sprintf("%d tasks across %d cores\n", nrow(TASKS), N_CORES))

  t0 <- Sys.time()

  run_parallel <- function(idx, fun, n_cores) {
    if (n_cores <= 1L) return(lapply(idx, fun))

    wd <- getwd()
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterExport(cl, "wd", envir = environment())
    parallel::clusterExport(cl, c("SPATIAL_LEVEL", "PREY_FAMILY",
                                  "N_MIN", "N_STRICT", "MIN_SETS",
                                  "R_PERM", "ALPHA"),
                            envir = parent.env(environment()))
    parallel::clusterEvalQ(cl, {
      setwd(wd)
      suppressPackageStartupMessages({
        library(dplyr); library(tidyr); library(purrr); library(tibble)
      })
      sink(nullfile())
      source("03_Statistical_Pipeline.R", local = .GlobalEnv)
      sink()
      TRUE
    })

    parallel::clusterExport(cl, c("diet_grid", "inv_src", "design_keep",
                                  "RUN_SC", "SEED", "OUT_DIR", "draw_sets",
                                  "TASKS", "SCENARIOS"),
                            envir = parent.env(environment()))

    parallel::parLapplyLB(cl, idx, fun)
  }

  res_list <- run_parallel(seq_len(nrow(TASKS)), run_one_task, N_CORES)
  reps <- Filter(Negate(is.null), res_list)

  all_runs_balanced <- bind_rows(reps)
  cat(sprintf("\nCompleted in %.1f min | %d rows\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins")),
              nrow(all_runs_balanced)))

  attr(all_runs_balanced, "design") <- design
  attr(all_runs_balanced, "config") <- list(
    B_REP = B_REP, SEED = SEED, X_GRID = X_GRID, TARGET_RULE = TARGET_RULE,
    BALANCE_ALL = BALANCE_ALL, MIN_SETS = MIN_SETS, R_PERM = R_PERM,
    SPATIAL_LEVEL = SPATIAL_LEVEL, PREY_FAMILY = PREY_FAMILY)

  if (B_REP > 1L) {
    chk <- all_runs_balanced %>%
      filter(testable) %>%
      count(sc_code, mode, rep_id, diagnostic) %>%
      group_by(sc_code, mode, diagnostic) %>%
      summarise(sd_n = sd(n), .groups = "drop") %>%
      group_by(sc_code, mode) %>%
      summarise(varies = any(sd_n > 0, na.rm = TRUE), .groups = "drop")
    cat("\nVariation among replicates by contrast:\n")
    print(as.data.frame(chk))
    if (any(!chk$varies))
      warning("Some contrasts do not vary across replicates: sampling ",
              "does not introduce variation in these cells.", call. = FALSE)
  }

  save(all_runs_balanced,
       file = sprintf("data/Sensitivity/all_runs_%s%s_BAL.rda",
                      SPATIAL_LEVEL, PREY_FAMILY))
  cat("Saved: data/Sensitivity/all_runs_", SPATIAL_LEVEL, PREY_FAMILY, "_BAL.rda\n", sep = "")
}

run_summaries <- function(level) {
  SPATIAL_LEVEL <- level
  PREY_FAMILY <- "_1"
  RUN_TAG <- "_T"
  N_CORES <- 2L
  # Summarize balanced samples
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(purrr); library(readr)
    library(ggplot2); library(tibble)
  })

  if (!exists("PREY_FAMILY"))   PREY_FAMILY   <- "_1"
  if (!exists("SPATIAL_LEVEL")) SPATIAL_LEVEL <- "all_gulf"

  OUT <- sprintf("Robustness_%s%s_BAL", SPATIAL_LEVEL, PREY_FAMILY)
  dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

  source("R_helpers/Results.R", local = TRUE)

  tidy_runs <- function(r) {
    r$contrast <- contrast_from_periods(r$period_1, r$period_2)
    r$contrast_type <- ifelse(as.character(r$contrast) %in% INTRA_CONTRASTS,
                              "within", "between")
    r$family <- as.character(family_of(r$diagnostic))
    unit <- intersect(c("Area", "str"), names(r))[1]
    r$cell <- paste(r$species, r$size_class, r[[unit]], sep = "|")
    attr(r, "unit_col") <- unit
    r
  }

  saved_file <- find_balanced_file(SPATIAL_LEVEL)
  if (is.na(saved_file)) {
    stop("Missing balanced results for ", SPATIAL_LEVEL,
         ": data/Sensitivity/all_runs_", SPATIAL_LEVEL,
         "_1_BAL.rda (or Sensitivity/). Run balanced resampling first.")
  }
  load(saved_file)
  message("Balanced results loaded: ", saved_file)
  bal <- tidy_runs(all_runs_balanced)
  UNIT <- attr(bal, "unit_col")
  CFG  <- attr(all_runs_balanced, "config")
  bal_cls <- bal %>% filter(testable, !is.na(family))

  n_rep <- n_distinct(bal$rep_id)
  message(sprintf("%d replicates | %d rows | %d cells",
                  n_rep, nrow(bal), n_distinct(bal_cls$cell)))

  pct_by_rep <- bal_cls %>%
    count(rep_id, contrast, contrast_type, mode, x_threshold,
          !!rlang::sym(UNIT), family, name = "n") %>%
    group_by(rep_id, contrast, contrast_type, mode, x_threshold, !!rlang::sym(UNIT)) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    group_by(rep_id, contrast, contrast_type, mode, x_threshold, family) %>%
    summarise(pct = mean(pct), .groups = "drop") %>%
    group_by(rep_id, contrast, contrast_type, mode, family) %>%
    summarise(pct = mean(pct), .groups = "drop")

  fam_summary <- pct_by_rep %>%
    group_by(contrast, contrast_type, mode, family) %>%
    summarise(mean_pct = mean(pct), sd_pct = sd(pct),
              q025 = quantile(pct, .025), q500 = median(pct),
              q975 = quantile(pct, .975), .groups = "drop") %>%
    mutate(family = factor(family, levels = FAMILY_LEVELS)) %>%
    arrange(mode, contrast, family)

  write_csv(fam_summary, file.path(OUT, "balanced_family_frequencies.csv"))
  write_csv(pct_by_rep,  file.path(OUT, "balanced_family_per_replicate.csv"))

  inconclusive <- bal %>%
    filter(testable) %>%
    group_by(rep_id, contrast, mode) %>%
    summarise(pct_inconclusive = 100 * mean(is.na(family)), .groups = "drop")
  write_csv(inconclusive, file.path(OUT, "balanced_inconclusive_rate.csv"))

  verdict <- pct_by_rep %>%
    filter(family %in% c("Stability", "Substitution")) %>%
    group_by(rep_id, mode, family, contrast_type) %>%
    summarise(lo = min(pct), hi = max(pct), .groups = "drop") %>%
    pivot_wider(names_from = contrast_type, values_from = c(lo, hi)) %>%
    mutate(separated = if_else(family == "Stability",
                               hi_between < lo_within, lo_between > hi_within),
           gap_pp    = if_else(family == "Stability",
                               lo_within - hi_between, lo_between - hi_within)) %>%
    group_by(mode, family) %>%
    summarise(n_rep = n(),
              pct_replicates_separated = 100 * mean(separated, na.rm = TRUE),
              median_gap_pp = median(gap_pp, na.rm = TRUE),
              q025_gap_pp   = quantile(gap_pp, .025, na.rm = TRUE), .groups = "drop")

  write_csv(verdict, file.path(OUT, "balanced_separation_verdict.csv"))
  cat("\n=================== VERDICT ===================\n")
  print(as.data.frame(verdict))
  cat("===============================================\n\n")

  ref_ok <- tryCatch({
    ref <- load_all_runs(SPATIAL_LEVEL, PREY_FAMILY, tag = RUN_TAG, verbose = FALSE) %>%
      filter(testable, !is.na(family), x_threshold %in% CFG$X_GRID)
    fam_ref <- ref %>%
      count(contrast, contrast_type, mode, x_threshold, !!rlang::sym(UNIT), family, name = "n") %>%
      group_by(contrast, contrast_type, mode, x_threshold, !!rlang::sym(UNIT)) %>%
      mutate(pct = 100 * n / sum(n)) %>%
      group_by(contrast, contrast_type, mode, x_threshold, family) %>%
      summarise(pct = mean(pct), .groups = "drop") %>%
      group_by(contrast, contrast_type, mode, family) %>%
      summarise(pct_unbalanced = mean(pct), .groups = "drop")

    comparison <- fam_summary %>%
      left_join(fam_ref, by = c("contrast", "contrast_type", "mode", "family")) %>%
      mutate(shift_pp = mean_pct - pct_unbalanced)
    write_csv(comparison, file.path(OUT, "balanced_vs_unbalanced.csv"))
    fam_ref
  }, error = function(e) { message("Unbalanced reference results unavailable: ",
                                   conditionMessage(e)); NULL })

  p <- ggplot(pct_by_rep %>% mutate(family = factor(family, levels = FAMILY_LEVELS)),
              aes(contrast, pct, fill = contrast_type)) +
    geom_violin(scale = "width", alpha = .55, colour = NA) +
    geom_boxplot(width = .14, outlier.size = .3, alpha = .9) +
    facet_grid(mode ~ family, scales = "free_y") +
    labs(x = NULL, y = "Classifiable cells (%)", fill = "Contrast",
         title = sprintf("Balanced trawl-set effort - %s%s", SPATIAL_LEVEL, PREY_FAMILY),
         subtitle = sprintf("%d replicates%s", n_rep,
                            if (!is.null(ref_ok)) "; crosses = unbalanced estimates" else "")) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (!is.null(ref_ok)) {
    p <- p + geom_point(data = ref_ok %>% mutate(family = factor(family, levels = FAMILY_LEVELS)),
                        aes(contrast, pct_unbalanced), inherit.aes = FALSE,
                        shape = 4, size = 2.6, stroke = 1)
  }
  ggsave(file.path(OUT, "balanced_family_distributions.png"), p,
         width = 11, height = 6, dpi = 300)

  n_sc <- n_distinct(bal_cls$contrast)
  common <- bal_cls %>%
    distinct(cell, x_threshold, mode, contrast, rep_id) %>%
    count(cell, x_threshold, mode, name = "n_obs") %>%
    filter(n_obs == n_rep * n_sc)

  fam_common <- bal_cls %>%
    semi_join(common, by = c("cell", "x_threshold", "mode")) %>%
    count(rep_id, contrast, contrast_type, mode, x_threshold, family, name = "n") %>%
    group_by(rep_id, contrast, contrast_type, mode, x_threshold) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    group_by(contrast, contrast_type, mode, family) %>%
    summarise(mean_pct = mean(pct), q025 = quantile(pct, .025),
              q975 = quantile(pct, .975), .groups = "drop")
  write_csv(fam_common, file.path(OUT, "balanced_common_cells.csv"))
  message(sprintf("Common cells: %d", n_distinct(common$cell)))

  modal <- function(x) { t <- table(x); names(t)[which.max(t)] }
  UNITQ <- rlang::sym(UNIT)

  consensus <- bal %>%
    filter(testable) %>%
    group_by(species, size_class, !!UNITQ, x_threshold, mode, contrast) %>%
    summarise(
      prey_family = dplyr::first(prey_family),
      period_1 = dplyr::first(period_1), period_2 = dplyr::first(period_2),
      spatial_level = dplyr::first(spatial_level),
      n_rep_cell = dplyr::n(),
      testable = dplyr::n() >= 0.5 * n_rep,
      diagnostic = modal(diagnostic),
      rep_support = mean(diagnostic == modal(diagnostic)),
      across(any_of(c("p_comp", "p_H", "p_Bs", "p_disp", "BC", "R2_comp",
                      "delta_H", "delta_Bs", "t_H", "t_Bs", "se_Bs",
                      "n_sto_P1", "n_sto_P2", "n_set_P1", "n_set_P2")),
             ~ median(.x, na.rm = TRUE)),
      reliable = mean(reliable) >= 0.5,
      .groups = "drop")

  res <- list(results = as.data.frame(consensus),
              inventory = attr(all_runs_balanced, "design"),
              n_comp_failures = 0L)
  assign("res_balanced_consensus", res)
  save(res_balanced_consensus,
       file = sprintf("data/Sensitivity/all_runs_%s%s_BALCONS.rda",
                      SPATIAL_LEVEL, PREY_FAMILY))

  message("Balanced summaries saved: ", OUT,
          " | consensus : data/Sensitivity/all_runs_", SPATIAL_LEVEL,
          PREY_FAMILY, "_BALCONS.rda")
}

run_null <- function(level) {
  SPATIAL_LEVEL <- level
  PREY_FAMILY <- "_1"
  RUN_TAG <- "_T"
  N_CORES <- 2L
  # Permute period labels within cells
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(purrr); library(readr)
    library(ggplot2); library(tibble)
  })

  if (!exists("PREY_FAMILY"))   PREY_FAMILY   <- "_1"
  if (!exists("SPATIAL_LEVEL")) SPATIAL_LEVEL <- "all_gulf"

  N_MIN <- 5; N_STRICT <- 25; MIN_SETS <- 3; R_PERM <- 999; ALPHA <- 0.05

  B_PERM <- 200
  SEED   <- 20260911
  X_GRID <- seq(100, 1000, by = 100)
  SC     <- "PT"

  OUT <- sprintf("NullCalibration_%s%s", SPATIAL_LEVEL, PREY_FAMILY)
  dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

  source("R_helpers/Results.R", local = TRUE)

  null_file <- sprintf("data/null_runs_%s%s.rda", SPATIAL_LEVEL, PREY_FAMILY)
  if (file.exists(null_file)) {
    saved <- new.env(parent = baseenv())
    objects <- load(null_file, envir = saved)
    if (!"null_runs" %in% objects || !is.data.frame(saved$null_runs) ||
        nrow(saved$null_runs) == 0L) {
      stop("Invalid saved null results: ", null_file)
    }
    null_runs <- saved$null_runs
    message("Using saved null results: ", null_file)
  } else {
    source("03_Statistical_Pipeline.R", local = TRUE)

  load("data/dat_classed.rda")
  diet_clean <- as.data.frame(dat_classed)
  if (!is.na(SPATIAL_SOURCE)) {
    if (!SPATIAL_SOURCE %in% names(diet_clean))
      stop("Missing column '", SPATIAL_SOURCE, "' in dat_classed.")
    diet_clean[[SPATIAL_SOURCE]] <- as.factor(diet_clean[[SPATIAL_SOURCE]])
  }

  prey_all  <- grep(paste0("^prey_category_\\d+", PREY_FAMILY, "$"),
                    names(diet_clean), value = TRUE)
  prey_keep <- intersect(paste0("prey_category_", X_GRID, PREY_FAMILY), prey_all)
  if (!length(prey_keep))
    stop("No X_GRID prey columns found. Available: ",
         paste(head(prey_all, 15), collapse = ", "))
  diet_grid <- diet_clean[, setdiff(names(diet_clean), setdiff(prey_all, prey_keep)),
                          drop = FALSE]
  cat(sprintf("Resolutions: %d columns retained from %d\n",
              length(prey_keep), length(prey_all)))

  pm  <- SCENARIOS[[SC]]
  dat <- make_dat_classed(diet_grid, period_map = pm)
  dat <- prepare_spatial_strata(dat)
  dat$set_uid <- paste(dat$year, dat[["vessel.code"]], dat$set, sep = "_")

  PERIOD_LEVELS <- levels(dat$period)

  key <- dat %>%
    distinct(sp = predator_species_common_name,
             sz = as.character(size_class),
             un = as.character(.data[[SPATIAL_OUT]]),
             set_uid, period)

  permute_labels <- function() {
    key %>%
      group_by(sp, sz, un) %>%
      mutate(period_perm = sample(period)) %>%
      ungroup()
  }

  apply_labels <- function(kp) {
    dat %>%
      select(-period) %>%
      inner_join(
        kp %>% select(sp, sz, un, set_uid, period = period_perm),
        by = c("predator_species_common_name" = "sp",
               "size_class" = "sz", "set_uid" = "set_uid")) %>%
      filter(as.character(.data[[SPATIAL_OUT]]) == un) %>%
      select(-un) %>%
      mutate(period = factor(as.character(period), levels = PERIOD_LEVELS))
  }

  run_one_perm <- function(b) {
    set.seed(SEED + 1e5 + b)
    d_p <- apply_labels(permute_labels())
    out <- list()
    for (md in c("biomass", "occurrence")) {
      r <- tryCatch(
        run_pipeline(d_p, mode = md,
                     period_1 = names(pm)[1], period_2 = names(pm)[2])$results,
        error = function(e) { warning("perm ", b, " / ", md, " : ",
                                      conditionMessage(e), call. = FALSE); NULL })
      if (is.null(r) || !nrow(r)) next
      r$perm_id <- b
      out[[length(out) + 1L]] <- r
    }
    bind_rows(out)
  }

  cat(sprintf("\nNull model: %d permutations x 2 currencies x %d resolutions\n",
              B_PERM, length(prey_keep)))

  t0 <- Sys.time()


  run_parallel <- function(idx, fun, n_cores) {
    if (n_cores <= 1L) return(lapply(idx, fun))
    wd <- getwd()
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, "wd", envir = environment())
    parallel::clusterExport(cl, c("SPATIAL_LEVEL", "PREY_FAMILY",
                                  "N_MIN", "N_STRICT", "MIN_SETS", "R_PERM", "ALPHA"),
                            envir = parent.env(environment()))
    parallel::clusterEvalQ(cl, {
      setwd(wd)
      suppressPackageStartupMessages({
        library(dplyr); library(tidyr); library(purrr); library(tibble)
      })
      source("03_Statistical_Pipeline.R", local = .GlobalEnv)
      TRUE
    })
    parallel::clusterExport(cl, c("dat", "key", "pm", "SEED", "PERIOD_LEVELS",
                                  "permute_labels", "apply_labels"),
                            envir = parent.env(environment()))
    parallel::parLapplyLB(cl, idx, fun)
  }

  null_runs <- bind_rows(run_parallel(seq_len(B_PERM), run_one_perm, N_CORES))
  cat(sprintf("\nCompleted in %.1f min | %d rows\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(null_runs)))

  save(null_runs, file = null_file)
  }

  nr <- null_runs %>% filter(testable)

  fpr <- nr %>%
    group_by(mode) %>%
    summarise(n = n(),
              fpr_comp = 100 * mean(p_comp < ALPHA, na.rm = TRUE),
              fpr_H    = 100 * mean(p_H    < ALPHA, na.rm = TRUE),
              fpr_Bs   = 100 * mean(p_Bs   < ALPHA, na.rm = TRUE),
              fpr_disp = 100 * mean(p_disp < ALPHA, na.rm = TRUE), .groups = "drop")
  write_csv(fpr, file.path(OUT, "null_false_positive_rates.csv"))
  cat("\n--- Null rejection rates (nominal 5%) ---\n")
  print(as.data.frame(fpr))

  fpr_by_n <- nr %>%
    mutate(n_min = pmin(n_set_P1, n_set_P2),
           n_bin = cut(n_min, breaks = c(2, 3, 5, 8, 12, 20, Inf), right = FALSE)) %>%
    group_by(mode, n_bin) %>%
    summarise(n = n(),
              `Composition (PERMANOVA)` = 100 * mean(p_comp < ALPHA, na.rm = TRUE),
              `Diversite (H')`          = 100 * mean(p_H    < ALPHA, na.rm = TRUE),
              `Amplitude (Bs)`          = 100 * mean(p_Bs   < ALPHA, na.rm = TRUE),
              .groups = "drop")
  write_csv(fpr_by_n, file.path(OUT, "null_fpr_by_sample_size.csv"))

  p_cal <- fpr_by_n %>%
    pivot_longer(c(`Composition (PERMANOVA)`, `Diversite (H')`, `Amplitude (Bs)`),
                 names_to = "test", values_to = "fpr") %>%
    ggplot(aes(n_bin, fpr, colour = test, group = test)) +
    geom_hline(yintercept = 100 * ALPHA, linetype = 2, colour = "grey40") +
    geom_line() + geom_point() + facet_wrap(~ mode) +
    labs(x = "Trawl sets per period (minimum)",
         y = "Null rejection rate (%)", colour = NULL,
         title = "Calibration of three tests under permuted period labels",
         subtitle = "Dashed line: nominal 5% level") +
    theme_bw(base_size = 10)
  ggsave(file.path(OUT, "null_calibration.png"), p_cal, width = 9, height = 4.5, dpi = 300)

  nf <- nr %>%
    mutate(family = as.character(family_of(diagnostic))) %>%
    filter(!is.na(family)) %>%
    count(perm_id, mode, x_threshold, family) %>%
    group_by(perm_id, mode, x_threshold) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    group_by(perm_id, mode, family) %>%
    summarise(pct = mean(pct), .groups = "drop")

  null_fam <- nf %>%
    group_by(mode, family) %>%
    summarise(null_mean = mean(pct), null_q025 = quantile(pct, .025),
              null_q975 = quantile(pct, .975), .groups = "drop")
  write_csv(null_fam, file.path(OUT, "null_family_distribution.csv"))
  cat("\n--- Null distribution of diagnostic families ---\n")
  print(as.data.frame(null_fam))

  cat("\nCompare the observed Substitution frequency (M0) with its null\n",
      "upper 97.5% null quantile and the relative calibration of the tests.\n", sep = "")
}

run_effort <- function(level) {
  SPATIAL_LEVEL <- level
  PREY_FAMILY <- "_1"
  RUN_TAG <- "_T"
  N_CORES <- 2L
  # Summarize sampling effort
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(purrr); library(readr)
    library(ggplot2); library(tibble)
  })

  if (!exists("PREY_FAMILY"))   PREY_FAMILY   <- "_1"
  if (!exists("SPATIAL_LEVEL")) SPATIAL_LEVEL <- "all_gulf"
  if (!exists("RUN_TAG"))       RUN_TAG       <- "_T"

  OUT <- sprintf("Effort_%s%s", SPATIAL_LEVEL, PREY_FAMILY)
  dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

  source("R_helpers/Results.R", local = TRUE)

  runs <- load_all_runs(SPATIAL_LEVEL, PREY_FAMILY, tag = RUN_TAG)
  UNIT <- attr(runs, "unit_col")

  runs <- runs %>%
    filter(testable, !is.na(family)) %>%
    mutate(n_min = pmin(n_set_P1, n_set_P2),
           n_sto_min = pmin(n_sto_P1, n_sto_P2))

  pct_fam <- function(d) {
    d %>%
      count(contrast, contrast_type, mode, x_threshold, !!rlang::sym(UNIT),
            family, name = "n") %>%
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

  per_cell <- runs %>%
    group_by(contrast, contrast_type, mode, cell) %>%
    summarise(n_set_P1 = median(n_set_P1), n_set_P2 = median(n_set_P2),
              n_min = median(n_min), n_sto_min = median(n_sto_min),
              .groups = "drop")

  effort <- per_cell %>%
    group_by(contrast, contrast_type, mode) %>%
    summarise(
      n_cells    = n(),
      sets_P1_med = median(n_set_P1),
      sets_P2_med = median(n_set_P2),
      sets_min_med = median(n_min),
      sets_min_q25 = quantile(n_min, .25),
      sets_min_q75 = quantile(n_min, .75),
      sto_min_med  = median(n_sto_min),
      .groups = "drop")

  write_csv(effort, file.path(OUT, "A_effort_by_contrast.csv"))
  cat("\n=== A. Sampling effort by contrast (sets per period) ===\n")
  print(as.data.frame(effort %>%
    select(contrast, contrast_type, mode, n_cells,
           sets_P1_med, sets_P2_med, sets_min_med, sets_min_q25, sets_min_q75)))

  gradient <- runs %>%
    mutate(n_bin = cut(n_min, breaks = c(2, 4, 6, 9, 14, 25, Inf), right = FALSE)) %>%
    filter(!is.na(n_bin)) %>%
    count(contrast, contrast_type, mode, n_bin, family, name = "n") %>%
    group_by(contrast, contrast_type, mode, n_bin) %>%
    mutate(pct = 100 * n / sum(n), n_tot = sum(n)) %>%
    ungroup()

  write_csv(gradient, file.path(OUT, "B_family_by_effort_class.csv"))

  slopes <- runs %>%
    mutate(is_stab = as.integer(family == "Stability")) %>%
    group_by(contrast, contrast_type, mode) %>%
    group_modify(function(d, k) {
      if (n_distinct(d$is_stab) < 2 || n_distinct(d$n_min) < 3) return(tibble())
      fit <- try(glm(is_stab ~ log(n_min), family = binomial, data = d), silent = TRUE)
      if (inherits(fit, "try-error")) return(tibble())
      s <- summary(fit)$coefficients
      tibble(slope_log_nsets = s[2, 1], se = s[2, 2],
             p = s[2, 4], n_obs = nrow(d))
    }) %>% ungroup()

  write_csv(slopes, file.path(OUT, "B_stability_vs_effort_slope.csv"))
  cat("\n=== B. Stability against log sampling effort ===\n")
  cat("A significant negative slope would indicate lower Stability\n")
  cat("at higher sampling effort, consistent with a power effect.\n\n")
  print(as.data.frame(slopes))

  p_grad <- gradient %>%
    filter(family %in% c("Stability", "Substitution")) %>%
    ggplot(aes(n_bin, pct, colour = contrast, group = contrast,
               linetype = contrast_type)) +
    geom_line() + geom_point(size = 1.2) +
    facet_grid(mode ~ family) +
    labs(x = "Trawl sets per period (minimum)",
         y = "Classifiable cells (%)", colour = NULL, linetype = NULL,
         title = "Diagnostic families by sampling effort and temporal contrast",
         subtitle = "Flat trends suggest limited sensitivity to sampling effort") +
    theme_bw(base_size = 10)
  ggsave(file.path(OUT, "B_effort_gradient.png"), p_grad,
         width = 10, height = 6, dpi = 300)

  common_support <- function(q_lo, q_hi) {
    per_cell %>%
      group_by(contrast, mode) %>%
      summarise(lo = quantile(n_min, q_lo), hi = quantile(n_min, q_hi),
                .groups = "drop") %>%
      group_by(mode) %>%
      summarise(lo_common = max(lo), hi_common = min(hi), .groups = "drop")
  }

  support <- NULL
  for (qq in list(c(.10, .90), c(.05, .95), c(0, 1))) {
    cand <- common_support(qq[1], qq[2])
    if (all(cand$lo_common <= cand$hi_common)) {
      support <- cand %>% mutate(q_lo = qq[1], q_hi = qq[2]); break
    }
  }

  cat("\n=== C. Common sampling-effort support ===\n")
  if (is.null(support)) {
    cat("No common support: the contrasts have no overlapping effort range.\n")
    cat("Comparable-effort restriction is not possible.\n")
    cat("Use the balanced resampling analysis for this comparison.\n")
    runs_common <- runs[0, ]
  } else {
    print(as.data.frame(support))
    cat(sprintf("Percentile range: %.0f%% - %.0f%%\n",
                100 * support$q_lo[1], 100 * support$q_hi[1]))
    runs_common <- runs %>%
      left_join(support %>% select(mode, lo_common, hi_common), by = "mode") %>%
      filter(n_min >= lo_common, n_min <= hi_common)
  }

  if (!nrow(runs_common)) {
    message("Common-support analysis inconclusive: no eligible cells.")
  }

  n_kept <- runs_common %>%
    group_by(contrast, mode) %>%
    summarise(cells_kept = n_distinct(cell), .groups = "drop") %>%
    left_join(per_cell %>% count(contrast, mode, name = "cells_total"),
              by = c("contrast", "mode")) %>%
    mutate(pct_kept = round(100 * cells_kept / cells_total, 1))
  write_csv(n_kept, file.path(OUT, "C_cells_retained.csv"))
  cat("\nCells retained in the common-support range:\n")
  print(as.data.frame(n_kept))

  comparison <- pct_fam(runs) %>% rename(pct_all = pct) %>%
    left_join(pct_fam(runs_common) %>% rename(pct_common = pct),
              by = c("contrast", "contrast_type", "mode", "family")) %>%
    mutate(pct_common = coalesce(pct_common, 0),
           shift_pp = pct_common - pct_all,
           family = factor(family, levels = FAMILY_LEVELS)) %>%
    arrange(mode, contrast, family)
  write_csv(comparison, file.path(OUT, "C_common_support_vs_all.csv"))

  verdict <- comparison %>%
    filter(family %in% c("Stability", "Substitution")) %>%
    group_by(mode, family, contrast_type) %>%
    summarise(lo_all = min(pct_all), hi_all = max(pct_all),
              lo_com = min(pct_common), hi_com = max(pct_common), .groups = "drop") %>%
    pivot_wider(names_from = contrast_type,
                values_from = c(lo_all, hi_all, lo_com, hi_com)) %>%
    mutate(
      gap_all = if_else(family == "Stability",
                        lo_all_within - hi_all_between, lo_all_between - hi_all_within),
      gap_common = if_else(family == "Stability",
                           lo_com_within - hi_com_between, lo_com_between - hi_com_within))
  write_csv(verdict, file.path(OUT, "C_separation_verdict.csv"))

  cat("\n=== Gap between within- and between-period ranges ===\n")
  print(as.data.frame(verdict %>% select(mode, family, gap_all, gap_common)))
  cat("\nA positive common-support gap indicates that separation persists\n")
  cat("after restricting comparisons to similar sampling effort.\n")

  message("Effort summaries saved: ", OUT)
}

run_calibrated <- function(level) {
  if (level != "all_gulf") stop("Calibration is available for all_gulf only")
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(ggplot2)
  })

  source("R_helpers/Results.R", local = TRUE)

  null_path <- Sys.getenv("TROPHIC_NULL_RDA",
                          file.path("data", "null_runs_all_gulf_1.rda"))
  observed_path <- Sys.getenv("TROPHIC_OBS_RDA",
                              file.path("data", "Sensitivity", "all_runs_all_gulf_1_T.rda"))
  output_dir <- "Calibrated_all_gulf_1_T"

  read_saved_object <- function(path, preferred) {
    if (!file.exists(path)) {
      stop("Missing local input: ", path,
           ". Keep DFO-controlled files outside the public repository.")
    }
    e <- new.env(parent = emptyenv())
    loaded <- load(path, envir = e)
    if (preferred %in% loaded && is.data.frame(e[[preferred]]))
      return(e[[preferred]])
    if (length(loaded) == 1L && is.data.frame(e[[loaded[[1L]]]]))
      return(e[[loaded[[1L]]]])
    pieces <- lapply(loaded, function(nm) {
      item <- e[[nm]]
      if (is.list(item) && is.data.frame(item$results)) return(item$results)
      if (is.data.frame(item)) return(item)
      NULL
    })
    pieces <- Filter(Negate(is.null), pieces)
    if (!length(pieces)) stop("No analysis table in ", path)
    dplyr::bind_rows(pieces)
  }

  null_runs <- read_saved_object(null_path, "null_runs")
  observed <- read_saved_object(observed_path, "all_runs_all_gulf_1")

  required_null <- c("mode", "testable", "p_comp", "p_H", "p_Bs", "n_set_P1", "n_set_P2", "x_threshold")
  required_obs <- c("mode", "testable", "p_comp", "p_H", "p_Bs", "diagnostic",
                    "period_1", "period_2", "x_threshold")
  if (length(missing <- setdiff(required_null, names(null_runs)))) {
    stop("Missing null-result columns: ", paste(missing, collapse = ", "))
  }
  if (length(missing <- setdiff(required_obs, names(observed)))) {
    stop("Missing observed-result columns: ", paste(missing, collapse = ", "))
  }
  if (!setequal(unique(as.character(null_runs$mode)), c("biomass", "occurrence"))) {
    stop("Unexpected dietary currencies in null simulations.")
  }

  if ("perm_id" %in% names(null_runs) &&
      dplyr::n_distinct(null_runs$perm_id) != 200L) {
    stop("Expected 200 null permutations in the supplied saved runs.")
  }
  if (dplyr::n_distinct(null_runs$x_threshold) != 10L) {
    stop("Expected ten prey thresholds in the supplied null simulations.")
  }

  nr <- null_runs %>% filter(.data$testable)
  if (!nrow(nr)) stop("No testable null results.")

  alpha <- nr %>%
    group_by(.data$mode) %>%
    summarise(
      n_null = n(),
      alpha_C = unname(quantile(.data$p_comp, 0.05, na.rm = TRUE)),
      alpha_H = unname(quantile(.data$p_H, 0.05, na.rm = TRUE)),
      alpha_Bs = unname(quantile(.data$p_Bs, 0.05, na.rm = TRUE)),
      trend_C = unname(quantile(.data$p_comp, 0.10, na.rm = TRUE)),
      fpr_C_at_05 = 100 * mean(.data$p_comp < 0.05, na.rm = TRUE),
      fpr_H_at_05 = 100 * mean(.data$p_H < 0.05, na.rm = TRUE),
      fpr_Bs_at_05 = 100 * mean(.data$p_Bs < 0.05, na.rm = TRUE),
      .groups = "drop")

  alpha_by_size <- nr %>%
    mutate(n_min = pmin(.data$n_set_P1, .data$n_set_P2),
           n_bin = cut(.data$n_min, breaks = c(2, 3, 5, 8, 12, 20, Inf),
                       right = FALSE)) %>%
    group_by(.data$mode, .data$n_bin) %>%
    summarise(
      n = n(),
      alpha_C = unname(quantile(.data$p_comp, 0.05, na.rm = TRUE)),
      alpha_H = unname(quantile(.data$p_H, 0.05, na.rm = TRUE)),
      alpha_Bs = unname(quantile(.data$p_Bs, 0.05, na.rm = TRUE)),
      .groups = "drop")

  obs <- observed %>%
    filter(.data$testable) %>%
    mutate(contrast = contrast_from_periods(.data$period_1, .data$period_2),
           contrast_type = if_else(as.character(.data$contrast) %in% INTRA_CONTRASTS,
                                   "within", "between"),
           family_nominal = as.character(family_of(.data$diagnostic))) %>%
    left_join(alpha %>% select(.data$mode, .data$alpha_C, .data$alpha_H,
                               .data$alpha_Bs, .data$trend_C),
              by = "mode") %>%
    mutate(
      complete_tests = !is.na(.data$p_comp) & !is.na(.data$p_H) & !is.na(.data$p_Bs),
      C = .data$p_comp < .data$alpha_C,
      H = .data$p_H < .data$alpha_H,
      B = .data$p_Bs < .data$alpha_Bs,
      family_calibrated = case_when(
        !.data$complete_tests ~ NA_character_,
        .data$C & (.data$H | .data$B) ~ "Reorganisation",
        .data$C                        ~ "Substitution",
        .data$H | .data$B               ~ "Functional change",
        TRUE                            ~ "Stability"
      ))

  if (anyNA(obs$contrast)) stop("One or more period contrasts could not be identified.")
  if (!nrow(obs)) stop("No testable observed results.")

  frequency_table <- function(data, family_column) {
    data %>%
      filter(!is.na(.data[[family_column]])) %>%
      mutate(family = factor(.data[[family_column]], levels = FAMILY_LEVELS)) %>%
      count(.data$contrast, .data$contrast_type, .data$mode,
            .data$x_threshold, .data$family, name = "n") %>%
      group_by(.data$contrast, .data$contrast_type, .data$mode, .data$x_threshold) %>%
      tidyr::complete(family = factor(FAMILY_LEVELS, levels = FAMILY_LEVELS),
                      fill = list(n = 0L)) %>%
      mutate(pct = 100 * .data$n / sum(.data$n)) %>%
      ungroup() %>%
      group_by(.data$contrast, .data$contrast_type, .data$mode, .data$family) %>%
      summarise(pct = mean(.data$pct), .groups = "drop")
  }
  nom <- frequency_table(obs, "family_nominal") %>%
    rename(pct_nominal = pct)
  cal <- frequency_table(obs, "family_calibrated") %>%
    rename(pct_calibrated = pct)

  compare <- nom %>%
    full_join(cal, by = c("contrast", "contrast_type", "mode", "family")) %>%
    mutate(shift_pp = .data$pct_calibrated - .data$pct_nominal) %>%
    arrange(.data$mode, .data$contrast, .data$family)

  if (nrow(compare) != 40L || anyNA(compare$pct_calibrated)) {
    stop("Expected 40 complete combinations (5 contrasts x 2 currencies x 4 families).")
  }

  separation <- compare %>%
    filter(.data$family %in% c("Stability", "Substitution")) %>%
    group_by(.data$mode, .data$family, .data$contrast_type) %>%
    summarise(lo_nom = min(.data$pct_nominal),
              hi_nom = max(.data$pct_nominal),
              lo_cal = min(.data$pct_calibrated),
              hi_cal = max(.data$pct_calibrated),
              .groups = "drop") %>%
    pivot_wider(names_from = contrast_type,
                values_from = c(lo_nom, hi_nom, lo_cal, hi_cal),
                names_glue = "{.value}_{contrast_type}") %>%
    transmute(
      mode = .data$mode, family = .data$family,
      lo_nom_between = .data$lo_nom_between, lo_nom_within = .data$lo_nom_within,
      hi_nom_between = .data$hi_nom_between, hi_nom_within = .data$hi_nom_within,
      lo_cal_between = .data$lo_cal_between, lo_cal_within = .data$lo_cal_within,
      hi_cal_between = .data$hi_cal_between, hi_cal_within = .data$hi_cal_within,
      gap_nominal = if_else(.data$family == "Stability",
                            .data$lo_nom_within - .data$hi_nom_between,
                            .data$lo_nom_between - .data$hi_nom_within),
      gap_calibrated = if_else(.data$family == "Stability",
                               .data$lo_cal_within - .data$hi_cal_between,
                               .data$lo_cal_between - .data$hi_cal_within))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(alpha, file.path(output_dir, "calibrated_alphas.csv"))
  write_csv(alpha_by_size, file.path(output_dir, "calibrated_alphas_by_sample_size.csv"))
  write_csv(compare, file.path(output_dir, "calibrated_vs_nominal.csv"))
  write_csv(separation, file.path(output_dir, "calibrated_separation_verdict.csv"))

  plot_data <- compare %>%
    pivot_longer(cols = c(pct_nominal, pct_calibrated),
                 names_to = "calibration", values_to = "pct") %>%
    mutate(calibration = recode(.data$calibration,
                                pct_nominal = "Nominal",
                                pct_calibrated = "Calibrated"),
           family = factor(.data$family, levels = FAMILY_LEVELS))
  p <- ggplot(plot_data,
              aes(x = .data$contrast, y = .data$pct,
                  color = .data$calibration, group = .data$calibration)) +
    geom_line(linewidth = 0.5) + geom_point(size = 1.4) +
    facet_grid(mode ~ family) +
    labs(x = "Temporal contrast", y = "Family frequency (%)",
         color = NULL,
         title = "Diagnostic families before and after null-based calibration") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(output_dir, "calibrated_family_frequencies.png"),
         p, width = 12, height = 6, dpi = 300)

  message("Calibration reconstructed using ", nrow(nr),
          " null rows and ", nrow(obs), " observed rows.")
  print(alpha)
  message("Saved statistical summaries in: ", normalizePath(output_dir))
  message("Calibrated summaries are available.")
}


# Sequential robustness analyses
SPATIAL_LEVELS <- c("all_gulf", "ecoregion", "stratum")
failed <- character()
run_status <- data.frame(step = character(), status = character(),
                         stringsAsFactors = FALSE)

file_ready <- function(path) {
  file.exists(path) && !is.na(file.info(path)$size) && file.info(path)$size > 0
}

all_ready <- function(paths) {
  length(paths) > 0L && all(vapply(paths, file_ready, logical(1)))
}

effort_outputs <- function(level) {
  file.path(sprintf("Effort_%s_1", level), c(
    "A_effort_by_contrast.csv",
    "B_family_by_effort_class.csv",
    "B_stability_vs_effort_slope.csv",
    "B_effort_gradient.png",
    "C_cells_retained.csv",
    "C_common_support_vs_all.csv",
    "C_separation_verdict.csv"
  ))
}

balanced_summary_outputs <- function(level) {
  c(file.path(sprintf("Robustness_%s_1_BAL", level), c(
      "balanced_family_frequencies.csv",
      "balanced_family_per_replicate.csv",
      "balanced_inconclusive_rate.csv",
      "balanced_separation_verdict.csv",
      "balanced_vs_unbalanced.csv",
      "balanced_family_distributions.png",
      "balanced_common_cells.csv"
    )),
    sprintf("data/Sensitivity/all_runs_%s_1_BALCONS.rda", level))
}

null_outputs <- function() {
  c("data/null_runs_all_gulf_1.rda",
    file.path("NullCalibration_all_gulf_1", c(
      "null_false_positive_rates.csv",
      "null_fpr_by_sample_size.csv",
      "null_family_distribution.csv",
      "null_calibration.png")))
}

calibrated_outputs <- function() {
  file.path("Calibrated_all_gulf_1_T", c(
    "calibrated_alphas.csv",
    "calibrated_alphas_by_sample_size.csv",
    "calibrated_vs_nominal.csv",
    "calibrated_separation_verdict.csv",
    "calibrated_family_frequencies.png"
  ))
}

run_step <- function(label, task, completed, required = character()) {
  message("\n--- ", label, " ---")
  if (isTRUE(completed())) {
    message("SKIPPED: saved outputs already exist")
    run_status <<- rbind(run_status,
                         data.frame(step = label, status = "SKIPPED"))
    return(invisible(NULL))
  }

  missing <- required[!vapply(required, file_ready, logical(1))]
  if (length(missing)) {
    detail <- paste("Missing inputs:", paste(missing, collapse = ", "))
    message("BLOCKED: ", detail)
    failed <<- c(failed, paste0(label, ": ", detail))
    run_status <<- rbind(run_status,
                         data.frame(step = label, status = "BLOCKED"))
    return(invisible(NULL))
  }

  result <- tryCatch({
    task()
    if (!isTRUE(completed())) {
      stop("Some expected output files are missing or empty after execution")
    }
    message("COMPLETED: outputs saved")
    "COMPLETED"
  }, error = function(e) {
    failed <<- c(failed, paste0(label, ": ", conditionMessage(e)))
    message("FAILED: ", conditionMessage(e))
    "FAILED"
  })
  run_status <<- rbind(run_status, data.frame(step = label, status = result))
  invisible(NULL)
}

for (level in SPATIAL_LEVELS) {
  local({
    this_level <- level

    run_step(paste("Sampling effort", this_level),
             function() run_effort(this_level),
             function() all_ready(effort_outputs(this_level)))

    run_step(paste("Balanced resampling", this_level),
             function() run_balanced(this_level),
             function() !is.na(find_balanced_file(this_level)),
             required = if (is.na(find_balanced_file(this_level))) {
               c("data/dat_classed.rda", "03_Statistical_Pipeline.R")
             } else character())

    run_step(paste("Balanced summaries", this_level),
             function() run_summaries(this_level),
             function() all_ready(balanced_summary_outputs(this_level)),
             required = if (is.na(find_balanced_file(this_level))) {
               sprintf("data/Sensitivity/all_runs_%s_1_BAL.rda", this_level)
             } else character())
  })
}

null_saved <- "data/null_runs_all_gulf_1.rda"
run_step("Null simulations and summaries",
         function() run_null("all_gulf"),
         function() all_ready(null_outputs()),
         required = if (file_ready(null_saved)) character() else {
           c("data/dat_classed.rda", "03_Statistical_Pipeline.R")
         })

run_step("Calibrated results",
         function() run_calibrated("all_gulf"),
         function() all_ready(calibrated_outputs()),
         required = c(
           Sys.getenv("TROPHIC_NULL_RDA", "data/null_runs_all_gulf_1.rda"),
           Sys.getenv("TROPHIC_OBS_RDA",
                      "data/Sensitivity/all_runs_all_gulf_1_T.rda")
         ))

message("\n--- Script 07 status ---")
print(run_status, row.names = FALSE)
if (length(failed)) {
  warning("Incomplete robustness workflow:\n", paste("-", failed, collapse = "\n"),
          call. = FALSE)
} else {
  message("Robustness workflow complete. Existing results were preserved.")
}
