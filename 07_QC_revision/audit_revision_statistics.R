
library(data.table)

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

OUT <- "results/QC_revision"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

all_audit <- list()
top_events <- list()

for (nm in names(files)) {

  cat("\n")
  cat("============================================================\n")
  cat(nm, "\n")
  cat("============================================================\n")

  d <- fread(files[[nm]])

  cat("\nDimensions:\n")
  cat(nrow(d), "rows x", ncol(d), "columns\n")

  cat("\nColumn names:\n")
  print(names(d))

  # ----------------------------------------------------------
  # Basic significance
  # ----------------------------------------------------------

  d[, P_value := as.numeric(P_value)]
  d[, FDR := as.numeric(FDR)]
  d[, DeltaPAU := as.numeric(DeltaPAU)]

  d[, nominal := !is.na(P_value) & P_value < 0.05]
  d[, fdr_sig := !is.na(FDR) & FDR < 0.05]

  cat("\nP < 0.05 events:", sum(d$nominal), "\n")
  cat("FDR < 0.05 events:", sum(d$fdr_sig), "\n")

  # ----------------------------------------------------------
  # Show categorical columns
  # ----------------------------------------------------------

  possible_cat <- names(d)[
    grepl(
      "direction|mechan|class|position|site|type|utr",
      names(d),
      ignore.case = TRUE
    )
  ]

  cat("\nPotential classification columns:\n")
  print(possible_cat)

  for (cc in possible_cat) {

    cat("\n---------------- ", cc, " ----------------\n")

    x <- d[nominal == TRUE, .N, by = cc][order(-N)]

    print(x)
  }

  # ----------------------------------------------------------
  # Direction column
  # ----------------------------------------------------------

  direction_col <- names(d)[
    grepl(
      "^Direction$|APA.*Direction|Mechanism|Classification",
      names(d),
      ignore.case = TRUE
    )
  ]

  if (length(direction_col) > 0) {

    direction_col <- direction_col[1]

    cat("\n")
    cat("Direction column used:", direction_col, "\n")

    cat("\nAll events:\n")
    print(
      d[, .N, by = direction_col][order(-N)]
    )

    cat("\nNominal P < 0.05:\n")
    print(
      d[nominal == TRUE,
        .N,
        by = direction_col][order(-N)]
    )

    cat("\nFDR < 0.05:\n")
    print(
      d[fdr_sig == TRUE,
        .N,
        by = direction_col][order(-N)]
    )
  }

  # ----------------------------------------------------------
  # DeltaPAU distributions
  # ----------------------------------------------------------

  cat("\nDeltaPAU summary — all:\n")
  print(summary(d$DeltaPAU))

  cat("\nDeltaPAU summary — nominal P < 0.05:\n")
  print(summary(d[nominal == TRUE]$DeltaPAU))

  cat("\nDeltaPAU summary — FDR < 0.05:\n")
  print(summary(d[fdr_sig == TRUE]$DeltaPAU))

  # ----------------------------------------------------------
  # Top 10 nominal events
  # ----------------------------------------------------------

  gene_cols <- names(d)[
    grepl(
      "^Gene$|GeneName|Gene_Name|Symbol",
      names(d),
      ignore.case = TRUE
    )
  ]

  keep <- unique(c(
    gene_cols,
    possible_cat,
    "DeltaPAU",
    "P_value",
    "FDR"
  ))

  keep <- keep[keep %in% names(d)]

  top <- d[
    nominal == TRUE
  ][
    order(-abs(DeltaPAU))
  ][
    1:min(.N, 10),
    ..keep
  ]

  top[, Dataset := nm]

  top_events[[nm]] <- top

  # ----------------------------------------------------------
  # Audit row
  # ----------------------------------------------------------

  all_audit[[nm]] <- data.table(
    Dataset = nm,
    Events_total = nrow(d),
    Nominal_events = sum(d$nominal),
    FDR_events = sum(d$fdr_sig),
    Median_DeltaPAU_nominal =
      median(d[nominal == TRUE]$DeltaPAU, na.rm = TRUE),
    Median_abs_DeltaPAU_nominal =
      median(abs(d[nominal == TRUE]$DeltaPAU), na.rm = TRUE),
    Min_DeltaPAU_nominal =
      min(d[nominal == TRUE]$DeltaPAU, na.rm = TRUE),
    Max_DeltaPAU_nominal =
      max(d[nominal == TRUE]$DeltaPAU, na.rm = TRUE)
  )
}

audit <- rbindlist(all_audit, fill = TRUE)
tops <- rbindlist(top_events, fill = TRUE)

fwrite(
  audit,
  file.path(
    OUT,
    "05_statistics_audit_summary.csv"
  )
)

fwrite(
  tops,
  file.path(
    OUT,
    "06_top_nominal_APA_events.csv"
  )
)

cat("\n\n")
cat("============================================================\n")
cat("AUDIT SUMMARY\n")
cat("============================================================\n")

print(audit)

