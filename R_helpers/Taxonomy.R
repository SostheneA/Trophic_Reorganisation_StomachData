# Taxonomic ranks and prey aggregation

ranknfile <- c("superdomain", "domain", "subdomain",
               "hyperkingdom", "superkingdom","kingdom","subkingdom","infrakingdom", "parvkingdom",
               "superphylum","phylum", "subphylum","infraphylum", "parvphylum",
               "gigaclass", "megaclass", "superclass","class","subclass","infraclass", "subterclass", "parvclass",
               "superdivision", "division", "subdivision", "infradivision",
               "superlegion","legion","sublegion","infralegion",
               "megacohort", "supercohort","cohort","subcohort","infracohort",
               "gigaorder", "magnorder","grandorder", "superorder","order","suborder","infraorder","parvorder",
               "section","subsection",
               "gigafamily", "superfamily","family","subfamily","infrafamily",
               "supertribe","tribe","subtribe","infratribe",
               "suprageneric_name","genus","subgenus","infragenus",

               "species_aggregate","species","infraspecific_name",
               "grex","subspecies","cultivar_group","convariety",
               "infrasubspecific_name","proles","race","natio","aberration","morph",
               "variety","subvariety","form","subform","pathovar","biovar","chemovar",
               "morphovar","phagovar","serovar","chemoform","forma_specialis",
               "cultivar","strain","other","unranked")

rank_counts <- function(dt, ranks, target_rows = NULL) {
  idx <- if (is.null(target_rows)) seq_len(nrow(dt)) else which(target_rows)
  out <- list()
  for (rank in ranks) {
    v <- dt[[rank]]
    rows <- idx[!is.na(v[idx]) & v[idx] != ""]
    if (!length(rows)) next
    cn <- dt[rows, .(n_stomach  = uniqueN(unlist(list_stomach_id)),
                     n_predator = uniqueN(unlist(list_predator_id))),
             by = c(rank)]
    data.table::setnames(cn, rank, "name")
    out[[rank]] <- cn
  }
  out
}

PreyCategory <- function(dt,
                         target_ranks = c("phylum", "subphylum", "class", "subclass", "order",
                                          "infraorder", "family", "genus", "species"),
                         full_taxonomy,
                         target_rows = NULL,
                         n_stomach_threshold = 100,
                         n_predator_threshold = 2,
                         counts = NULL) {

  n <- nrow(dt)
  level <- rep(NA_character_, n)
  name  <- rep(NA_character_, n)

  all_idx <- if (is.null(target_rows)) seq_len(n) else which(target_rows)
  rem <- all_idx
  if (!length(rem)) return(list(level = level, name = name))

  valid_cols  <- intersect(full_taxonomy[full_taxonomy %in% target_ranks], names(dt))
  spec_to_gen <- rev(valid_cols)
  if (is.null(counts)) counts <- rank_counts(dt, spec_to_gen, target_rows)

  for (rank in spec_to_gen) {
    if (!length(rem)) break
    cn <- counts[[rank]]
    if (is.null(cn)) next
    valid <- cn[n_stomach >= n_stomach_threshold & n_predator >= n_predator_threshold, name]
    if (!length(valid)) next
    v   <- as.character(dt[[rank]])
    hit <- rem[!is.na(v[rem]) & v[rem] != "" & v[rem] %chin% valid]
    if (length(hit)) {
      level[hit] <- rank
      name[hit]  <- v[hit]
      rem <- setdiff(rem, hit)
    }
  }

  if (length(rem)) {
    m <- lapply(valid_cols, function(cc) { x <- as.character(dt[[cc]][rem]); x[x == ""] <- NA_character_; x })
    level[rem] <- valid_cols[1]
    name[rem]  <- do.call(fcoalesce, m)
  }

  list(level = level, name = name)
}
