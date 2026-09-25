# ============================================================
# TDP-43 KO QAPA GO enrichment analysis
#
# Converted from the full TDP-43 K263E GO-enrichment workflow.
#
# Input:
#   TDP43_KO_QAPA_PAU.tsv
#
# Comparison:
#   TDP-43 KO vs Control
#
# Samples:
#   Control: Cont_S1, Cont_S2, Cont_S3, Cont_S4
#   KO:      TDP43_S5, TDP43_S6, TDP43_S7
#
# Mechanistic APA classification from QAPA APA_ID:
#
#   Proximal PAS increased  -> Shortening
#   Distal PAS decreased    -> Shortening
#
#   Proximal PAS decreased  -> Lengthening
#   Distal PAS increased    -> Lengthening
#
# Event criteria:
#   nominal P < 0.05
#   |DeltaPAU| >= 10
#
# IMPORTANT:
#   _S (single-site) events are not used to call shortening/
#   lengthening because they do not define a proximal/distal switch.
#
# GO background:
#   All QAPA genes represented by P or D sites.
#
# Main output needed by the replot script:
#   GO_analysis/TDP43_KO_GO_enrichment.xlsx
#
# ============================================================


# ============================================================
# 1. Packages
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(gprofiler2)
  library(ggplot2)
  library(openxlsx)
})


# ============================================================
# 2. Paths
# ============================================================

BASE_DIR <-
  "results/01_TDP43_KD"

INPUT <- file.path(
  BASE_DIR,
  "TDP43_KO_QAPA_PAU.tsv"
)

OUT_DIR <- file.path(
  BASE_DIR,
  "GO_analysis"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


EVENT_TABLE_OUT <- file.path(
  OUT_DIR,
  "TDP43_KO_mechanistic_APA_events.tsv"
)

GENE_LIST_XLSX <- file.path(
  OUT_DIR,
  "TDP43_KO_GO_gene_lists.xlsx"
)

GO_EXCEL_OUT <- file.path(
  OUT_DIR,
  "TDP43_KO_GO_enrichment.xlsx"
)

SHORT_GENES_OUT <- file.path(
  OUT_DIR,
  "TDP43_KO_shortening_genes.txt"
)

LONG_GENES_OUT <- file.path(
  OUT_DIR,
  "TDP43_KO_lengthening_genes.txt"
)

DISCORDANT_OUT <- file.path(
  OUT_DIR,
  "TDP43_KO_discordant_genes.txt"
)

BACKGROUND_OUT <- file.path(
  OUT_DIR,
  "TDP43_KO_QAPA_background_genes.txt"
)


# ============================================================
# 3. Parameters
# ============================================================

P_CUTOFF <- 0.05
DPAU_CUTOFF <- 10
TOP_GO <- 15


# ============================================================
# 4. Read QAPA table
# ============================================================

cat("\n========================================\n")
cat("READING TDP-43 KO QAPA TABLE\n")
cat("========================================\n")

df <- read.delim(
  INPUT,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Total QAPA events:", nrow(df), "\n")
cat("Unique genes:", length(unique(df$Gene)), "\n")


# ============================================================
# 5. Sample columns
# ============================================================

CONTROL_PAU <- c(
  "Cont_S1.PAU",
  "Cont_S2.PAU",
  "Cont_S3.PAU",
  "Cont_S4.PAU"
)

KO_PAU <- c(
  "TDP43_S5.PAU",
  "TDP43_S6.PAU",
  "TDP43_S7.PAU"
)

CONTROL_TPM <- c(
  "Cont_S1.TPM",
  "Cont_S2.TPM",
  "Cont_S3.TPM",
  "Cont_S4.TPM"
)

KO_TPM <- c(
  "TDP43_S5.TPM",
  "TDP43_S6.TPM",
  "TDP43_S7.TPM"
)

ALL_PAU <- c(
  CONTROL_PAU,
  KO_PAU
)

ALL_TPM <- c(
  CONTROL_TPM,
  KO_TPM
)


# ============================================================
# 6. Check required columns
# ============================================================

required_cols <- c(
  "APA_ID",
  "Gene",
  "Gene_Name",
  ALL_PAU,
  ALL_TPM
)

missing_cols <- setdiff(
  required_cols,
  colnames(df)
)

if (length(missing_cols) > 0) {
  stop(
    paste(
      "Missing required columns:",
      paste(missing_cols, collapse = ", ")
    )
  )
}


# ============================================================
# 7. Convert numeric sample columns
# ============================================================

for (x in c(
  ALL_PAU,
  ALL_TPM,
  "Length",
  "Num_Events",
  "UTR3.Start",
  "UTR3.End"
)) {

  if (x %in% colnames(df)) {
    df[[x]] <- suppressWarnings(
      as.numeric(df[[x]])
    )
  }
}


# ============================================================
# 8. QAPA site class from APA_ID
#
# D = distal
# P = proximal
# S = single
# ============================================================

df$Site <- sub(
  ".*_",
  "",
  df$APA_ID
)

cat("\nQAPA PAS classes:\n")
print(
  table(df$Site)
)


# ============================================================
# 9. Mean PAU and DeltaPAU
#
# DeltaPAU = KO - Control
# ============================================================

df$Control_mean_PAU <- rowMeans(
  df[, CONTROL_PAU, drop = FALSE],
  na.rm = TRUE
)

df$KO_mean_PAU <- rowMeans(
  df[, KO_PAU, drop = FALSE],
  na.rm = TRUE
)

df$DeltaPAU <- (
  df$KO_mean_PAU -
  df$Control_mean_PAU
)


# ============================================================
# 10. Safe Welch t-test
# ============================================================

safe_ttest <- function(ko, control) {

  ko <- as.numeric(ko)
  control <- as.numeric(control)

  ko <- ko[
    is.finite(ko)
  ]

  control <- control[
    is.finite(control)
  ]

  if (
    length(ko) < 2 ||
    length(control) < 2
  ) {
    return(NA_real_)
  }

  if (
    length(unique(ko)) == 1 &&
    length(unique(control)) == 1
  ) {

    if (
      unique(ko) ==
      unique(control)
    ) {
      return(1)
    }
  }

  tryCatch(
    t.test(
      ko,
      control,
      paired = FALSE,
      var.equal = FALSE
    )$p.value,
    error = function(e) NA_real_
  )
}


# ============================================================
# 11. Event-level statistics
#
# This follows the original K263E GO-enrichment workflow:
# statistics are calculated across the QAPA PAU table directly.
# ============================================================

cat("\nCalculating event-level statistics...\n")

df$P_value <- vapply(
  seq_len(nrow(df)),
  function(i) {

    safe_ttest(
      df[i, KO_PAU],
      df[i, CONTROL_PAU]
    )
  },
  numeric(1)
)


# ============================================================
# 12. BH FDR
# ============================================================

df$FDR <- NA_real_

tested <- which(
  is.finite(df$P_value)
)

df$FDR[tested] <- p.adjust(
  df$P_value[tested],
  method = "BH"
)


# ============================================================
# 13. Mechanistic APA direction
#
# Shortening:
#   P increased
#   D decreased
#
# Lengthening:
#   P decreased
#   D increased
#
# Single-site S events are excluded from these calls.
# ============================================================

df$APA_direction <- "NS"

# Proximal increased -> shortening
df$APA_direction[
  df$Site == "P" &
  is.finite(df$P_value) &
  df$P_value < P_CUTOFF &
  df$DeltaPAU >= DPAU_CUTOFF
] <- "Shortening"

# Distal decreased -> shortening
df$APA_direction[
  df$Site == "D" &
  is.finite(df$P_value) &
  df$P_value < P_CUTOFF &
  df$DeltaPAU <= -DPAU_CUTOFF
] <- "Shortening"

# Proximal decreased -> lengthening
df$APA_direction[
  df$Site == "P" &
  is.finite(df$P_value) &
  df$P_value < P_CUTOFF &
  df$DeltaPAU <= -DPAU_CUTOFF
] <- "Lengthening"

# Distal increased -> lengthening
df$APA_direction[
  df$Site == "D" &
  is.finite(df$P_value) &
  df$P_value < P_CUTOFF &
  df$DeltaPAU >= DPAU_CUTOFF
] <- "Lengthening"


# ============================================================
# 14. Explicit four-category classification
# ============================================================

df$Mechanistic_category <- "NS"

df$Mechanistic_category[
  df$APA_direction == "Shortening" &
  df$Site == "P"
] <- "Proximal increased"

df$Mechanistic_category[
  df$APA_direction == "Shortening" &
  df$Site == "D"
] <- "Distal decreased"

df$Mechanistic_category[
  df$APA_direction == "Lengthening" &
  df$Site == "P"
] <- "Proximal decreased"

df$Mechanistic_category[
  df$APA_direction == "Lengthening" &
  df$Site == "D"
] <- "Distal increased"


cat("\nAPA events:\n")
print(
  table(df$APA_direction)
)

cat("\nFour mechanistic categories:\n")
print(
  table(df$Mechanistic_category)
)


# ============================================================
# 15. Save full event table
# ============================================================

write.table(
  df,
  EVENT_TABLE_OUT,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 16. Shortening and lengthening gene lists
# ============================================================

short_genes <- df %>%
  filter(
    APA_direction == "Shortening"
  ) %>%
  pull(
    Gene_Name
  ) %>%
  unique() %>%
  na.omit()

long_genes <- df %>%
  filter(
    APA_direction == "Lengthening"
  ) %>%
  pull(
    Gene_Name
  ) %>%
  unique() %>%
  na.omit()


# ============================================================
# 17. Remove discordant genes
#
# A gene with both shortening and lengthening events is removed
# from both clean GO query sets.
# ============================================================

discordant <- intersect(
  short_genes,
  long_genes
)

cat(
  "\nShortening genes before removing discordant:",
  length(short_genes),
  "\n"
)

cat(
  "Lengthening genes before removing discordant:",
  length(long_genes),
  "\n"
)

cat(
  "Discordant genes:",
  length(discordant),
  "\n"
)

short_genes_clean <- setdiff(
  short_genes,
  discordant
)

long_genes_clean <- setdiff(
  long_genes,
  discordant
)

cat(
  "\nFinal shortening genes:",
  length(short_genes_clean),
  "\n"
)

cat(
  "Final lengthening genes:",
  length(long_genes_clean),
  "\n"
)


# ============================================================
# 18. QAPA custom background
#
# Same logic as the original K263E GO analysis:
# all genes represented by multi-PAS P/D QAPA sites.
#
# Do NOT restrict this background by the later TPM filter,
# because that would change the enrichment universe and could
# change the GO pathways/ranking relative to the original
# K263E workflow.
# ============================================================

background <- df %>%
  filter(
    Site %in% c(
      "P",
      "D"
    )
  ) %>%
  pull(
    Gene_Name
  ) %>%
  unique() %>%
  na.omit()

cat(
  "QAPA P/D background genes:",
  length(background),
  "\n"
)


# ============================================================
# 19. Save gene lists
# ============================================================

writeLines(
  short_genes_clean,
  SHORT_GENES_OUT
)

writeLines(
  long_genes_clean,
  LONG_GENES_OUT
)

writeLines(
  discordant,
  DISCORDANT_OUT
)

writeLines(
  background,
  BACKGROUND_OUT
)


# ============================================================
# 20. Save gene-list workbook
# ============================================================

wb_gene <- createWorkbook()

addWorksheet(
  wb_gene,
  "Shortening"
)
writeData(
  wb_gene,
  "Shortening",
  data.frame(
    Gene = short_genes_clean
  )
)

addWorksheet(
  wb_gene,
  "Lengthening"
)
writeData(
  wb_gene,
  "Lengthening",
  data.frame(
    Gene = long_genes_clean
  )
)

addWorksheet(
  wb_gene,
  "Discordant"
)
writeData(
  wb_gene,
  "Discordant",
  data.frame(
    Gene = discordant
  )
)

addWorksheet(
  wb_gene,
  "QAPA_background"
)
writeData(
  wb_gene,
  "QAPA_background",
  data.frame(
    Gene = background
  )
)

saveWorkbook(
  wb_gene,
  GENE_LIST_XLSX,
  overwrite = TRUE
)


# ============================================================
# 21. GO enrichment
#
# significant = FALSE is intentional.
# It forces g:Profiler to return all tested terms, even if
# none survives multiple-testing correction. This ensures the
# workbook required by the publication replot script is made.
# ============================================================

run_GO <- function(genes) {

  if (
    length(genes) < 2
  ) {
    return(NULL)
  }

  gost(
    query = genes,
    organism = "hsapiens",
    ordered_query = FALSE,
    custom_bg = background,
    correction_method = "fdr",
    sources = c(
      "GO:BP",
      "GO:CC",
      "GO:MF"
    ),
    significant = FALSE
  )
}


# ============================================================
# 22. Run g:Profiler
# ============================================================

cat(
  "\n========================================\n"
)
cat(
  "RUNNING GO ENRICHMENT\n"
)
cat(
  "========================================\n"
)

GO_short <- run_GO(
  short_genes_clean
)

GO_long <- run_GO(
  long_genes_clean
)


# ============================================================
# 23. Extract GO results safely
# ============================================================

short_result <- data.frame()
long_result <- data.frame()

if (
  !is.null(GO_short) &&
  !is.null(GO_short$result) &&
  nrow(GO_short$result) > 0
) {
  short_result <- GO_short$result
}

if (
  !is.null(GO_long) &&
  !is.null(GO_long$result) &&
  nrow(GO_long$result) > 0
) {
  long_result <- GO_long$result
}

cat(
  "\nShortening GO terms returned:",
  nrow(short_result),
  "\n"
)

cat(
  "Lengthening GO terms returned:",
  nrow(long_result),
  "\n"
)

if (
  nrow(short_result) > 0 &&
  "significant" %in% colnames(short_result)
) {
  cat(
    "Shortening FDR-significant GO terms:",
    sum(
      short_result$significant,
      na.rm = TRUE
    ),
    "\n"
  )
}

if (
  nrow(long_result) > 0 &&
  "significant" %in% colnames(long_result)
) {
  cat(
    "Lengthening FDR-significant GO terms:",
    sum(
      long_result$significant,
      na.rm = TRUE
    ),
    "\n"
  )
}


# ============================================================
# 24. Write GO enrichment workbook
#
# The sheet names exactly match what
# TDP43_KO_GO_replot_from_original.R expects.
# ============================================================

wb_go <- createWorkbook()

sheet_names <- c(
  "Shortening_BP",
  "Shortening_CC",
  "Shortening_MF",
  "Lengthening_BP",
  "Lengthening_CC",
  "Lengthening_MF"
)

for (s in sheet_names) {
  addWorksheet(
    wb_go,
    s
  )
}


write_go_sheet <- function(
  wb,
  sheet,
  dat,
  source_name
) {

  if (
    is.null(dat) ||
    nrow(dat) == 0
  ) {
    return(invisible(NULL))
  }

  out <- dat %>%
    filter(
      source == source_name
    ) %>%
    arrange(
      p_value
    )

  if (
    nrow(out) > 0
  ) {
    writeData(
      wb,
      sheet,
      out
    )
  }

  invisible(NULL)
}


write_go_sheet(
  wb_go,
  "Shortening_BP",
  short_result,
  "GO:BP"
)

write_go_sheet(
  wb_go,
  "Shortening_CC",
  short_result,
  "GO:CC"
)

write_go_sheet(
  wb_go,
  "Shortening_MF",
  short_result,
  "GO:MF"
)

write_go_sheet(
  wb_go,
  "Lengthening_BP",
  long_result,
  "GO:BP"
)

write_go_sheet(
  wb_go,
  "Lengthening_CC",
  long_result,
  "GO:CC"
)

write_go_sheet(
  wb_go,
  "Lengthening_MF",
  long_result,
  "GO:MF"
)


# ------------------------------------------------------------
# Simple workbook formatting
# ------------------------------------------------------------

header_style <- createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  halign = "center",
  valign = "center",
  border = "Bottom"
)

for (s in sheet_names) {

  sheet_data <- readWorkbook(
    wb_go,
    sheet = s
  )

  if (
    ncol(sheet_data) > 0
  ) {
    addStyle(
      wb_go,
      sheet = s,
      style = header_style,
      rows = 1,
      cols = seq_len(
        ncol(sheet_data)
      ),
      gridExpand = TRUE
    )

    freezePane(
      wb_go,
      sheet = s,
      firstRow = TRUE
    )

    setColWidths(
      wb_go,
      sheet = s,
      cols = seq_len(
        min(
          ncol(sheet_data),
          25
        )
      ),
      widths = "auto"
    )
  }
}

saveWorkbook(
  wb_go,
  GO_EXCEL_OUT,
  overwrite = TRUE
)


# ============================================================
# 25. Optional quick GO lollipop plots
#
# These are NOT the final publication replot. They are only
# diagnostic plots so you can inspect the enrichment now.
#
# Final publication plots can then be generated by:
#   TDP43_KO_GO_replot_from_original.R
# ============================================================

plot_GO_quick <- function(
  res,
  direction,
  source_name,
  outfile,
  top_n = TOP_GO
) {

  if (
    is.null(res) ||
    nrow(res) == 0 ||
    !"p_value" %in% colnames(res)
  ) {
    return(NULL)
  }

  x <- res %>%
    filter(
      source == source_name,
      is.finite(p_value),
      p_value > 0
    ) %>%
    arrange(
      p_value
    ) %>%
    slice_head(
      n = top_n
    )

  if (
    nrow(x) == 0
  ) {
    return(NULL)
  }

  x <- x %>%
    mutate(
      GeneRatio =
        intersection_size /
        query_size,
      GeneCount =
        intersection_size
    )

  x$term_name <- factor(
    x$term_name,
    levels = rev(
      x$term_name
    )
  )

  ontology_text <- dplyr::case_when(
    source_name == "GO:BP" ~ "biological processes",
    source_name == "GO:CC" ~ "cellular components",
    source_name == "GO:MF" ~ "molecular functions",
    TRUE ~ source_name
  )

  p <- ggplot(
    x,
    aes(
      x = GeneRatio,
      y = term_name
    )
  ) +
    geom_segment(
      aes(
        x = 0,
        xend = GeneRatio,
        y = term_name,
        yend = term_name
      ),
      linewidth = 0.35,
      color = "grey50"
    ) +
    geom_point(
      aes(
        size = GeneCount,
        color = p_value
      )
    ) +
    scale_color_viridis_c(
      option = "viridis",
      direction = -1,
      name = "FDR adjusted\np-value"
    ) +
    scale_size_continuous(
      name = "Gene count",
      range = c(
        2.5,
        7
      )
    ) +
    labs(
      title = paste0(
        "GO ",
        ontology_text,
        " associated with TDP-43 KO 3'UTR ",
        tolower(direction),
        " genes"
      ),
      x = "GeneRatio",
      y = NULL
    ) +
    theme_bw(
      base_size = 8
    ) +
    theme(
      panel.background =
        element_rect(
          fill = "#EAF2F7",
          color = NA
        ),
      panel.grid.major.y =
        element_blank(),
      panel.grid.minor =
        element_blank(),
      plot.title =
        element_text(
          size = 8,
          hjust = 0.5
        ),
      axis.title =
        element_text(
          size = 8
        ),
      axis.text =
        element_text(
          size = 8,
          color = "black"
        ),
      legend.title =
        element_text(
          size = 8
        ),
      legend.text =
        element_text(
          size = 8
        )
    )

  ggsave(
    outfile,
    p,
    width = 8,
    height = 6,
    device = cairo_pdf
  )

  invisible(p)
}


plot_GO_quick(
  short_result,
  "Shortening",
  "GO:BP",
  file.path(
    OUT_DIR,
    "TDP43_KO_shortening_GO_BP_quick.pdf"
  )
)

plot_GO_quick(
  short_result,
  "Shortening",
  "GO:CC",
  file.path(
    OUT_DIR,
    "TDP43_KO_shortening_GO_CC_quick.pdf"
  )
)

plot_GO_quick(
  short_result,
  "Shortening",
  "GO:MF",
  file.path(
    OUT_DIR,
    "TDP43_KO_shortening_GO_MF_quick.pdf"
  )
)

plot_GO_quick(
  long_result,
  "Lengthening",
  "GO:BP",
  file.path(
    OUT_DIR,
    "TDP43_KO_lengthening_GO_BP_quick.pdf"
  )
)

plot_GO_quick(
  long_result,
  "Lengthening",
  "GO:CC",
  file.path(
    OUT_DIR,
    "TDP43_KO_lengthening_GO_CC_quick.pdf"
  )
)

plot_GO_quick(
  long_result,
  "Lengthening",
  "GO:MF",
  file.path(
    OUT_DIR,
    "TDP43_KO_lengthening_GO_MF_quick.pdf"
  )
)


# ============================================================
# 26. Final report
# ============================================================

cat(
  "\n========================================\n"
)
cat(
  "TDP-43 KO GO ANALYSIS COMPLETE\n"
)
cat(
  "========================================\n"
)

cat(
  "\nFinal shortening genes:",
  length(short_genes_clean),
  "\n"
)

cat(
  "Final lengthening genes:",
  length(long_genes_clean),
  "\n"
)

cat(
  "Discordant genes:",
  length(discordant),
  "\n"
)

cat(
  "QAPA P/D background genes:",
  length(background),
  "\n"
)

cat(
  "\nGO workbook created:\n",
  GO_EXCEL_OUT,
  "\n"
)

cat(
  "\nNext run:\n",
  "Rscript TDP43_KO_GO_replot_from_original.R\n"
)
