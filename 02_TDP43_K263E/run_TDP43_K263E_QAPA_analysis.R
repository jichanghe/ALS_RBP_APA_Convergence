# ============================================================
# TDP-43 K263E QAPA differential APA analysis
#
# Input:
# TDP43_K263E_QAPA_PAU.tsv
#
# Outputs:
# TDP43_K263E_differential_APA.xlsx
# TDP43_K263E_volcano.pdf
# TDP43_K263E_PCA.pdf
# TDP43_K263E_heatmap.pdf
# TDP43_K263E_shortening_lengthening.pdf
#
# Comparison:
# K263E vs WT
#
# Primary significance:
# FDR < 0.05 and |DeltaPAU| >= 10 percentage points
# ============================================================


# ============================================================
# 1. Packages
# ============================================================

packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2",
  "ggrepel",
  "openxlsx",
  "pheatmap"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing missing package: ", pkg)
    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
  }
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(openxlsx)
  library(pheatmap)
})


# ============================================================
# 2. Paths
# ============================================================

OUT_DIR <- "results/02_TDP43_K263E"

INPUT <- file.path(
  OUT_DIR,
  "TDP43_K263E_QAPA_PAU.tsv"
)

XLSX_OUT <- file.path(
  OUT_DIR,
  "TDP43_K263E_differential_APA.xlsx"
)

VOLCANO_OUT <- file.path(
  OUT_DIR,
  "TDP43_K263E_volcano.pdf"
)

PCA_OUT <- file.path(
  OUT_DIR,
  "TDP43_K263E_PCA.pdf"
)

HEATMAP_OUT <- file.path(
  OUT_DIR,
  "TDP43_K263E_heatmap.pdf"
)

SL_OUT <- file.path(
  OUT_DIR,
  "TDP43_K263E_shortening_lengthening.pdf"
)


# ============================================================
# 3. Parameters
# ============================================================

TPM_CUTOFF <- 1
MIN_TPM_SAMPLES <- 2

DPAU_CUTOFF <- 10
FDR_CUTOFF <- 0.05

TOP_LABELS <- 10
TOP_HEATMAP <- 50

# Gene-level weighted 3'UTR length threshold
# 100 nt difference required for shortening/lengthening call
LENGTH_CHANGE_CUTOFF <- 100


# ============================================================
# 4. Read QAPA PAU file
# ============================================================

message("Reading:")
message(INPUT)

df <- read.delim(
  INPUT,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

message("Rows loaded: ", nrow(df))
message("Unique genes: ", length(unique(df$Gene)))


# ============================================================
# 5. Define sample columns
# ============================================================

wt_pau <- c(
  "WT_Rep1.PAU",
  "WT_Rep2.PAU",
  "WT_Rep3.PAU"
)

mut_pau <- c(
  "K263E_Rep1.PAU",
  "K263E_Rep2.PAU",
  "K263E_Rep3.PAU"
)

wt_tpm <- c(
  "WT_Rep1.TPM",
  "WT_Rep2.TPM",
  "WT_Rep3.TPM"
)

mut_tpm <- c(
  "K263E_Rep1.TPM",
  "K263E_Rep2.TPM",
  "K263E_Rep3.TPM"
)

all_pau <- c(wt_pau, mut_pau)
all_tpm <- c(wt_tpm, mut_tpm)


# ============================================================
# 6. Convert numeric columns
# ============================================================

for (x in c(
  all_pau,
  all_tpm,
  "Length",
  "Num_Events",
  "UTR3.Start",
  "UTR3.End",
  "LastExon.Start",
  "LastExon.End"
)) {
  if (x %in% colnames(df)) {
    df[[x]] <- suppressWarnings(
      as.numeric(df[[x]])
    )
  }
}


# ============================================================
# 7. Expression filtering
#
# Keep APA events with TPM >= 1 in >=2 samples
# ============================================================

df$TPM_pass_n <- rowSums(
  as.data.frame(df[, all_tpm, drop = FALSE]) >= TPM_CUTOFF,
  na.rm = TRUE
)

df$Expression_Pass <- (
  df$TPM_pass_n >= MIN_TPM_SAMPLES
)

message(
  "Events passing TPM filter: ",
  sum(df$Expression_Pass),
  " / ",
  nrow(df)
)


# ============================================================
# 8. Mean PAU and DeltaPAU
# ============================================================

df$WT_mean_PAU <- rowMeans(
  df[, wt_pau, drop = FALSE],
  na.rm = TRUE
)

df$K263E_mean_PAU <- rowMeans(
  df[, mut_pau, drop = FALSE],
  na.rm = TRUE
)

df$DeltaPAU <- (
  df$K263E_mean_PAU -
  df$WT_mean_PAU
)


# ============================================================
# 9. Welch t-test for each APA event
# ============================================================

safe_ttest <- function(x, y) {

  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  if (length(x) < 2 || length(y) < 2) {
    return(NA_real_)
  }

  # If both groups have exactly no variation
  if (
    length(unique(x)) == 1 &&
    length(unique(y)) == 1
  ) {

    if (unique(x) == unique(y)) {
      return(1)
    }

    return(NA_real_)
  }

  tryCatch(
    t.test(
      x,
      y,
      paired = FALSE,
      var.equal = FALSE
    )$p.value,
    error = function(e) NA_real_
  )
}


message("Calculating event-level statistics...")

df$P_value <- vapply(
  seq_len(nrow(df)),
  function(i) {

    if (!df$Expression_Pass[i]) {
      return(NA_real_)
    }

    safe_ttest(
      as.numeric(df[i, mut_pau]),
      as.numeric(df[i, wt_pau])
    )
  },
  numeric(1)
)


# ============================================================
# 10. Multiple testing correction
# ============================================================

df$FDR <- NA_real_

tested <- which(
  df$Expression_Pass &
  is.finite(df$P_value)
)

df$FDR[tested] <- p.adjust(
  df$P_value[tested],
  method = "BH"
)


# ============================================================
# 11. Extract QAPA site class
#
# Typical APA_ID:
# *_P = proximal
# *_D = distal
# *_S = single
# ============================================================

df$APA_Class <- case_when(
  grepl("_P$", df$APA_ID) ~ "Proximal",
  grepl("_D$", df$APA_ID) ~ "Distal",
  grepl("_S$", df$APA_ID) ~ "Single",
  TRUE ~ "Other"
)


# ============================================================
# 12. Event significance
# ============================================================

df$Significance <- "NS"

df$Significance[
  df$Expression_Pass &
  !is.na(df$FDR) &
  df$FDR < FDR_CUTOFF &
  df$DeltaPAU >= DPAU_CUTOFF
] <- "Up"

df$Significance[
  df$Expression_Pass &
  !is.na(df$FDR) &
  df$FDR < FDR_CUTOFF &
  df$DeltaPAU <= -DPAU_CUTOFF
] <- "Down"

df$Significance <- factor(
  df$Significance,
  levels = c("Down", "NS", "Up")
)


# ============================================================
# 13. Summary statistics
# ============================================================

summary_table <- data.frame(
  Metric = c(
    "Total APA events",
    "Unique genes",
    "Events passing TPM filter",
    "Events statistically tested",
    "FDR < 0.05",
    "FDR < 0.05 and |DeltaPAU| >= 5",
    "FDR < 0.05 and |DeltaPAU| >= 10",
    "FDR < 0.05 and |DeltaPAU| >= 15",
    "FDR < 0.05 and |DeltaPAU| >= 20",
    "Significant increased PAU",
    "Significant decreased PAU"
  ),

  Value = c(
    nrow(df),
    length(unique(df$Gene)),
    sum(df$Expression_Pass),
    sum(is.finite(df$P_value)),
    sum(df$FDR < 0.05, na.rm = TRUE),

    sum(
      df$FDR < 0.05 &
      abs(df$DeltaPAU) >= 5,
      na.rm = TRUE
    ),

    sum(
      df$FDR < 0.05 &
      abs(df$DeltaPAU) >= 10,
      na.rm = TRUE
    ),

    sum(
      df$FDR < 0.05 &
      abs(df$DeltaPAU) >= 15,
      na.rm = TRUE
    ),

    sum(
      df$FDR < 0.05 &
      abs(df$DeltaPAU) >= 20,
      na.rm = TRUE
    ),

    sum(df$Significance == "Up"),

    sum(df$Significance == "Down")
  )
)


# ============================================================
# 14. Sort differential table
# ============================================================

df_out <- df %>%
  arrange(
    FDR,
    desc(abs(DeltaPAU))
  )


sig_df <- df_out %>%
  filter(
    Significance != "NS"
  )


# ============================================================
# 15. VOLCANO PLOT
# ============================================================

volcano <- df %>%
  filter(
    Expression_Pass,
    is.finite(P_value)
  )

volcano$minus_log10_P <- -log10(
  pmax(volcano$P_value, 1e-300)
)


# Label strongest genes separately in each direction

label_up <- volcano %>%
  filter(Significance == "Up") %>%
  arrange(FDR, desc(abs(DeltaPAU))) %>%
  distinct(Gene_Name, .keep_all = TRUE) %>%
  slice_head(n = TOP_LABELS)

label_down <- volcano %>%
  filter(Significance == "Down") %>%
  arrange(FDR, desc(abs(DeltaPAU))) %>%
  distinct(Gene_Name, .keep_all = TRUE) %>%
  slice_head(n = TOP_LABELS)

label_df <- bind_rows(
  label_up,
  label_down
)


p_volcano <- ggplot(
  volcano,
  aes(
    x = DeltaPAU,
    y = minus_log10_P
  )
) +

  geom_point(
    data = subset(volcano, Significance == "NS"),
    color = "grey70",
    size = 1.6,
    alpha = 0.6
  ) +

  geom_point(
    data = subset(volcano, Significance == "Down"),
    color = "#377EB8",
    size = 2.0,
    alpha = 0.85
  ) +

  geom_point(
    data = subset(volcano, Significance == "Up"),
    color = "#E41A1C",
    size = 2.0,
    alpha = 0.85
  ) +

  geom_vline(
    xintercept = c(
      -DPAU_CUTOFF,
      DPAU_CUTOFF
    ),
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey35"
  ) +

  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey35"
  ) +

  geom_text_repel(
    data = label_df,
    aes(label = Gene_Name),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0
  ) +

  labs(
    title = "TDP-43 K263E differential alternative polyadenylation",
    subtitle = paste0(
      "FDR < ",
      FDR_CUTOFF,
      " and |DeltaPAU| >= ",
      DPAU_CUTOFF
    ),
    x = "DeltaPAU (K263E - WT)",
    y = "-log10(P value)"
  ) +

  theme_classic(base_size = 13) +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    legend.position = "none"
  )


ggsave(
  VOLCANO_OUT,
  p_volcano,
  width = 7,
  height = 6.5,
  device = cairo_pdf
)


# ============================================================
# 16. PCA
#
# Restrict to expressed multi-APA events
# ============================================================

pca_df <- df %>%
  filter(
    Expression_Pass,
    Num_Events >= 2
  )

pca_mat <- as.matrix(
  pca_df[, all_pau, drop = FALSE]
)

rownames(pca_mat) <- make.unique(
  pca_df$APA_ID
)

# Remove rows containing NA
pca_mat <- pca_mat[
  apply(
    pca_mat,
    1,
    function(x) all(is.finite(x))
  ),
  ,
  drop = FALSE
]

# Remove invariant events
pca_var <- apply(
  pca_mat,
  1,
  var
)

pca_mat <- pca_mat[
  is.finite(pca_var) &
  pca_var > 0,
  ,
  drop = FALSE
]

# Keep most variable events if matrix is very large
if (nrow(pca_mat) > 5000) {

  vv <- apply(
    pca_mat,
    1,
    var
  )

  pca_mat <- pca_mat[
    order(
      vv,
      decreasing = TRUE
    )[1:5000],
    ,
    drop = FALSE
  ]
}

pca <- prcomp(
  t(pca_mat),
  center = TRUE,
  scale. = TRUE
)

variance <- (
  pca$sdev^2 /
  sum(pca$sdev^2)
) * 100

pca_plot_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
)

pca_plot_df$Group <- ifelse(
  grepl("^WT_", pca_plot_df$Sample),
  "WT",
  "TDP-43 K263E"
)

pca_plot_df$Group <- factor(
  pca_plot_df$Group,
  levels = c(
    "WT",
    "TDP-43 K263E"
  )
)


p_pca <- ggplot(
  pca_plot_df,
  aes(
    PC1,
    PC2,
    color = Group,
    label = Sample
  )
) +

  geom_point(
    size = 4
  ) +

  geom_text_repel(
    size = 3.5,
    show.legend = FALSE
  ) +

  stat_ellipse(
    aes(group = Group),
    type = "norm",
    linewidth = 0.8,
    linetype = 2,
    show.legend = FALSE
  ) +

  labs(
    title = "QAPA PAU principal component analysis",
    x = paste0(
      "PC1 (",
      round(variance[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(variance[2], 1),
      "%)"
    )
  ) +

  theme_classic(base_size = 13) +

  theme(
    plot.title = element_text(
      face = "bold"
    ),
    legend.title = element_blank()
  )


ggsave(
  PCA_OUT,
  p_pca,
  width = 7,
  height = 6,
  device = cairo_pdf
)


# ============================================================
# 17. Significant-event heatmap
# ============================================================

heat_df <- df %>%
  filter(
    Significance != "NS"
  ) %>%
  arrange(
    FDR,
    desc(abs(DeltaPAU))
  )

if (nrow(heat_df) == 0) {

  message(
    "No events meet primary FDR + DeltaPAU threshold.",
    " Using top events by P value for heatmap."
  )

  heat_df <- df %>%
    filter(
      Expression_Pass,
      is.finite(P_value)
    ) %>%
    arrange(
      P_value,
      desc(abs(DeltaPAU))
    )
}

heat_df <- heat_df %>%
  distinct(
    APA_ID,
    .keep_all = TRUE
  ) %>%
  slice_head(
    n = TOP_HEATMAP
  )


heat_mat <- as.matrix(
  heat_df[, all_pau, drop = FALSE]
)

rownames(heat_mat) <- make.unique(
  paste0(
    heat_df$Gene_Name,
    " | ",
    heat_df$APA_Class
  )
)


# row z-score
heat_z <- t(
  scale(
    t(heat_mat)
  )
)

heat_z[
  !is.finite(heat_z)
] <- 0


annotation_col <- data.frame(
  Group = c(
    rep("WT", 3),
    rep("TDP-43 K263E", 3)
  )
)

rownames(annotation_col) <- all_pau


pdf(
  HEATMAP_OUT,
  width = 8,
  height = 10,
  useDingbats = FALSE
)

pheatmap(
  heat_z,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 7,
  fontsize_col = 9,
  border_color = NA,
  main = paste0(
    "Top ",
    nrow(heat_z),
    " differential APA events"
  )
)

dev.off()


# ============================================================
# 18. Gene-level PAU-weighted mean 3'UTR length
#
# WeightedLength(sample) =
# sum(PAU_fraction * UTR Length)
#
# This directly estimates whether a gene shifts toward
# shorter or longer 3'UTR isoforms.
# ============================================================

multi_df <- df %>%
  filter(
    Expression_Pass,
    Num_Events >= 2,
    is.finite(Length),
    Length > 0
  )


calc_weighted_length <- function(dat, pau_col) {

  pau <- dat[[pau_col]]

  good <- (
    is.finite(pau) &
    is.finite(dat$Length)
  )

  if (sum(good) == 0) {
    return(NA_real_)
  }

  pau <- pau[good]
  lengths <- dat$Length[good]

  total <- sum(pau)

  if (!is.finite(total) || total <= 0) {
    return(NA_real_)
  }

  sum(
    (pau / total) *
    lengths
  )
}


gene_list <- split(
  multi_df,
  multi_df$Gene
)

message(
  "Calculating gene-level weighted 3'UTR lengths..."
)

length_list <- lapply(
  gene_list,
  function(g) {

    data.frame(
      Gene = unique(g$Gene)[1],
      Gene_Name = unique(g$Gene_Name)[1],
      Num_APA_Sites = nrow(g),

      WT_Rep1 = calc_weighted_length(g, "WT_Rep1.PAU"),
      WT_Rep2 = calc_weighted_length(g, "WT_Rep2.PAU"),
      WT_Rep3 = calc_weighted_length(g, "WT_Rep3.PAU"),

      K263E_Rep1 = calc_weighted_length(g, "K263E_Rep1.PAU"),
      K263E_Rep2 = calc_weighted_length(g, "K263E_Rep2.PAU"),
      K263E_Rep3 = calc_weighted_length(g, "K263E_Rep3.PAU"),

      stringsAsFactors = FALSE
    )
  }
)

length_df <- bind_rows(
  length_list
)


length_df$WT_mean_length <- rowMeans(
  length_df[, c(
    "WT_Rep1",
    "WT_Rep2",
    "WT_Rep3"
  )],
  na.rm = TRUE
)

length_df$K263E_mean_length <- rowMeans(
  length_df[, c(
    "K263E_Rep1",
    "K263E_Rep2",
    "K263E_Rep3"
  )],
  na.rm = TRUE
)

length_df$Delta_UTR_length <- (
  length_df$K263E_mean_length -
  length_df$WT_mean_length
)


length_df$P_value <- vapply(
  seq_len(nrow(length_df)),
  function(i) {

    safe_ttest(
      as.numeric(
        length_df[
          i,
          c(
            "K263E_Rep1",
            "K263E_Rep2",
            "K263E_Rep3"
          )
        ]
      ),
      as.numeric(
        length_df[
          i,
          c(
            "WT_Rep1",
            "WT_Rep2",
            "WT_Rep3"
          )
        ]
      )
    )
  },
  numeric(1)
)

length_df$FDR <- p.adjust(
  length_df$P_value,
  method = "BH"
)


length_df$Direction <- case_when(

  is.finite(length_df$Delta_UTR_length) &
  length_df$Delta_UTR_length <= -LENGTH_CHANGE_CUTOFF ~
    "Shortening",

  is.finite(length_df$Delta_UTR_length) &
  length_df$Delta_UTR_length >= LENGTH_CHANGE_CUTOFF ~
    "Lengthening",

  TRUE ~
    "Stable"
)


# ============================================================
# 19. Shortening / lengthening plot
# ============================================================

sl_counts <- length_df %>%
  filter(
    Direction != "Stable"
  ) %>%
  count(
    Direction,
    name = "Genes"
  )

# Ensure both directions exist
sl_counts <- data.frame(
  Direction = c(
    "Shortening",
    "Lengthening"
  )
) %>%
  left_join(
    sl_counts,
    by = "Direction"
  )

sl_counts$Genes[
  is.na(sl_counts$Genes)
] <- 0

sl_counts$Direction <- factor(
  sl_counts$Direction,
  levels = c(
    "Shortening",
    "Lengthening"
  )
)


p_sl <- ggplot(
  sl_counts,
  aes(
    x = Direction,
    y = Genes,
    fill = Direction
  )
) +

  geom_col(
    width = 0.65
  ) +

  geom_text(
    aes(label = Genes),
    vjust = -0.35,
    size = 4.5
  ) +

  scale_fill_manual(
    values = c(
      "Shortening" = "#377EB8",
      "Lengthening" = "#E41A1C"
    )
  ) +

  labs(
    title = "TDP-43 K263E alters 3'UTR length",
    subtitle = paste0(
      "Gene-level PAU-weighted mean 3'UTR length; |change| >= ",
      LENGTH_CHANGE_CUTOFF,
      " nt"
    ),
    x = NULL,
    y = "Number of genes"
  ) +

  theme_classic(base_size = 13) +

  theme(
    legend.position = "none",
    plot.title = element_text(
      face = "bold"
    )
  )


ggsave(
  SL_OUT,
  p_sl,
  width = 6,
  height = 5.5,
  device = cairo_pdf
)


# ============================================================
# 20. Excel workbook
# ============================================================

message("Writing Excel workbook...")

wb <- createWorkbook()


# -------------------------
# Summary
# -------------------------

addWorksheet(
  wb,
  "Summary"
)

writeData(
  wb,
  "Summary",
  summary_table
)


# -------------------------
# All differential events
# -------------------------

addWorksheet(
  wb,
  "All_APA_events"
)

writeData(
  wb,
  "All_APA_events",
  df_out
)


# -------------------------
# Significant events
# -------------------------

addWorksheet(
  wb,
  "Significant_APA"
)

writeData(
  wb,
  "Significant_APA",
  sig_df
)


# -------------------------
# Gene-level UTR length
# -------------------------

addWorksheet(
  wb,
  "UTR_length"
)

writeData(
  wb,
  "UTR_length",
  length_df %>%
    arrange(
      desc(abs(Delta_UTR_length))
    )
)


# -------------------------
# Shortening
# -------------------------

addWorksheet(
  wb,
  "Shortening"
)

writeData(
  wb,
  "Shortening",
  length_df %>%
    filter(
      Direction == "Shortening"
    ) %>%
    arrange(
      Delta_UTR_length
    )
)


# -------------------------
# Lengthening
# -------------------------

addWorksheet(
  wb,
  "Lengthening"
)

writeData(
  wb,
  "Lengthening",
  length_df %>%
    filter(
      Direction == "Lengthening"
    ) %>%
    arrange(
      desc(Delta_UTR_length)
    )
)


# ============================================================
# Excel formatting
# ============================================================

header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  fgFill = "#D9EAF7",
  border = "Bottom"
)

for (s in names(wb)) {

  addStyle(
    wb,
    sheet = s,
    style = header_style,
    rows = 1,
    cols = 1:ncol(
      readWorkbook(
        wb,
        sheet = s
      )
    ),
    gridExpand = TRUE
  )

  freezePane(
    wb,
    sheet = s,
    firstRow = TRUE
  )

  setColWidths(
    wb,
    sheet = s,
    cols = 1:50,
    widths = "auto"
  )
}


saveWorkbook(
  wb,
  XLSX_OUT,
  overwrite = TRUE
)


# ============================================================
# 21. Final report
# ============================================================

message("")
message("==============================================")
message("ANALYSIS COMPLETE")
message("==============================================")

message("")
message("Input:")
message(INPUT)

message("")
message("Outputs:")

message(XLSX_OUT)
message(VOLCANO_OUT)
message(PCA_OUT)
message(HEATMAP_OUT)
message(SL_OUT)

message("")
message("Summary:")

print(summary_table)

message("")
message("Gene-level 3'UTR direction:")

print(
  table(
    length_df$Direction
  )
)
