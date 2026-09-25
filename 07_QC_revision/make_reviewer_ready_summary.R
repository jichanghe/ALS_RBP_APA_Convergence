
library(data.table)

ROOT <- "."

OUT <- "results/QC_revision"

# ============================================================
# Existing master summary
# ============================================================

summary_file <- file.path(
  OUT,
  "01_revision_quantitative_summary.csv"
)

master <- fread(summary_file)

# ============================================================
# Original statistical files
# ============================================================

files <- list(

  TDP43_K263E =
    "data/QC_revision/TDP43_K263E_mechanistic_APA_events.tsv",

  FUS_P525L =
    "data/QC_revision/FUS_P525L_APA_site_mechanistic_results.tsv",

  MATR3_KO =
    "data/QC_revision/MATR3_KO_APA_site_mechanistic_results.tsv",

  TDP43_KO =
    "data/QC_revision/TDP43_KD_mechanistic_APA_events.tsv",

  FUS_KO =
    "data/QC_revision/FUS_KO_APA_site_mechanistic_results.tsv"
)

# ============================================================
# Function
# ============================================================

analyze_dataset <- function(dataset, f) {

  d <- fread(f)

  d[, P_value := as.numeric(P_value)]
  d[, FDR := as.numeric(FDR)]
  d[, DeltaPAU := as.numeric(DeltaPAU)]

  # ----------------------------------------------------------
  # Gene column
  # ----------------------------------------------------------

  gene_col <- if ("Gene_Name" %in% names(d)) {
    "Gene_Name"
  } else {
    "Gene"
  }

  # ----------------------------------------------------------
  # Direction column
  # ----------------------------------------------------------

  if ("APA_direction" %in% names(d)) {

    direction_col <- "APA_direction"

  } else if ("Mechanism" %in% names(d)) {

    direction_col <- "Mechanism"

  } else {

    stop(
      "No direction column in ",
      dataset
    )
  }

  direction <- as.character(
    d[[direction_col]]
  )

  nominal <- !is.na(d$P_value) &
             d$P_value < 0.05

  fdrsig <- !is.na(d$FDR) &
            d$FDR < 0.05

  # ----------------------------------------------------------
  # Gene counts
  # ----------------------------------------------------------

  genes_tested <- unique(
    d[[gene_col]][
      !is.na(d[[gene_col]]) &
      d[[gene_col]] != ""
    ]
  )

  nominal_genes <- unique(
    d[[gene_col]][nominal]
  )

  nominal_genes <- nominal_genes[
    !is.na(nominal_genes) &
    nominal_genes != ""
  ]

  fdr_genes <- unique(
    d[[gene_col]][fdrsig]
  )

  fdr_genes <- fdr_genes[
    !is.na(fdr_genes) &
    fdr_genes != ""
  ]

  # ----------------------------------------------------------
  # Counts
  # ----------------------------------------------------------

  nominal_n <- sum(nominal)
  fdr_n <- sum(fdrsig)

  nom_NS <- sum(
    nominal &
    direction == "NS",
    na.rm = TRUE
  )

  nom_short <- sum(
    nominal &
    direction == "Shortening",
    na.rm = TRUE
  )

  nom_long <- sum(
    nominal &
    direction == "Lengthening",
    na.rm = TRUE
  )

  fdr_NS <- sum(
    fdrsig &
    direction == "NS",
    na.rm = TRUE
  )

  fdr_short <- sum(
    fdrsig &
    direction == "Shortening",
    na.rm = TRUE
  )

  fdr_long <- sum(
    fdrsig &
    direction == "Lengthening",
    na.rm = TRUE
  )

  # ----------------------------------------------------------
  # Percent helper
  # ----------------------------------------------------------

  pct <- function(x, denominator) {

    if (
      is.na(denominator) ||
      denominator == 0
    ) {
      return(NA_real_)
    }

    round(
      100 * x / denominator,
      2
    )
  }

  # ----------------------------------------------------------
  # Effect sizes
  # ----------------------------------------------------------

  med_dpau_nom <- if (nominal_n > 0) {
    median(
      d$DeltaPAU[nominal],
      na.rm = TRUE
    )
  } else {
    NA_real_
  }

  med_abs_dpau_nom <- if (nominal_n > 0) {
    median(
      abs(d$DeltaPAU[nominal]),
      na.rm = TRUE
    )
  } else {
    NA_real_
  }

  med_dpau_fdr <- if (fdr_n > 0) {
    median(
      d$DeltaPAU[fdrsig],
      na.rm = TRUE
    )
  } else {
    NA_real_
  }

  med_abs_dpau_fdr <- if (fdr_n > 0) {
    median(
      abs(d$DeltaPAU[fdrsig]),
      na.rm = TRUE
    )
  } else {
    NA_real_
  }

  data.table(

    Dataset = dataset,

    Events_tested = nrow(d),

    Genes_tested =
      length(genes_tested),

    Nominal_events =
      nominal_n,

    Nominal_genes =
      length(nominal_genes),

    Nominal_NS =
      nom_NS,

    Nominal_shortening =
      nom_short,

    Nominal_lengthening =
      nom_long,

    Nominal_NS_percent =
      pct(nom_NS, nominal_n),

    Nominal_shortening_percent =
      pct(nom_short, nominal_n),

    Nominal_lengthening_percent =
      pct(nom_long, nominal_n),

    FDR_events =
      fdr_n,

    FDR_genes =
      length(fdr_genes),

    FDR_NS =
      fdr_NS,

    FDR_shortening =
      fdr_short,

    FDR_lengthening =
      fdr_long,

    FDR_NS_percent =
      pct(fdr_NS, fdr_n),

    FDR_shortening_percent =
      pct(fdr_short, fdr_n),

    FDR_lengthening_percent =
      pct(fdr_long, fdr_n),

    Median_DeltaPAU_nominal =
      med_dpau_nom,

    Median_abs_DeltaPAU_nominal =
      med_abs_dpau_nom,

    Median_DeltaPAU_FDR =
      med_dpau_fdr,

    Median_abs_DeltaPAU_FDR =
      med_abs_dpau_fdr
  )
}

# ============================================================
# Analyze all
# ============================================================

results <- rbindlist(
  lapply(
    names(files),
    function(x) {
      analyze_dataset(
        x,
        files[[x]]
      )
    }
  ),
  fill = TRUE
)

# ============================================================
# Add sample sizes and PCA
# ============================================================

keep_master <- master[
  ,
  .(
    Dataset,
    Control_n,
    Experimental_n,
    Total_samples,

    PCA_PC1_scaled_percent,
    PCA_PC2_scaled_percent,

    PCA_PC1_unscaled_percent,
    PCA_PC2_unscaled_percent
  )
]

final <- merge(
  keep_master,
  results,
  by = "Dataset",
  all.y = TRUE
)

# Desired ordering
order_levels <- c(
  "TDP43_KO",
  "TDP43_K263E",
  "FUS_KO",
  "FUS_P525L",
  "MATR3_KO"
)

final[
  ,
  sort_order := match(
    Dataset,
    order_levels
  )
]

setorder(
  final,
  sort_order
)

final[, sort_order := NULL]

# ============================================================
# Internal consistency checks
# ============================================================

final[
  ,
  Nominal_classified_total :=
    Nominal_NS +
    Nominal_shortening +
    Nominal_lengthening
]

final[
  ,
  FDR_classified_total :=
    FDR_NS +
    FDR_shortening +
    FDR_lengthening
]

final[
  ,
  Nominal_count_OK :=
    Nominal_classified_total ==
    Nominal_events
]

final[
  ,
  FDR_count_OK :=
    FDR_classified_total ==
    FDR_events
]

# ============================================================
# Save
# ============================================================

outfile <- file.path(
  OUT,
  "07_REVIEWER_READY_quantitative_summary.csv"
)

fwrite(
  final,
  outfile
)

cat("\n")
cat("============================================================\n")
cat("REVIEWER-READY TABLE\n")
cat("============================================================\n\n")

print(final)

cat("\nSaved:\n")
cat(outfile, "\n")

# ============================================================
# Fail loudly if classifications do not reconcile
# ============================================================

if (!all(final$Nominal_count_OK)) {

  warning(
    "Some nominal counts do not reconcile."
  )
}

if (!all(final$FDR_count_OK)) {

  warning(
    "Some FDR counts do not reconcile."
  )
}

