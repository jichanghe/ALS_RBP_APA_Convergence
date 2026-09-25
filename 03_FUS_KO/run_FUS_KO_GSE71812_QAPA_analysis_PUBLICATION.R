# ============================================================
# FUS KO GSE71812 QAPA differential APA analysis
#
# Input:
# FUS_KO_QAPA_PAU.tsv
#
# Outputs:
# FUS_KO_differential_APA.xlsx
# FUS_KO_volcano.pdf
# FUS_KO_PCA.pdf
# FUS_KO_heatmap.pdf
# FUS_KO_shortening_lengthening.pdf
#
# Comparison:
# FUS KO vs WT
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
  "pheatmap",
  "scales"
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
  library(scales)
})


# ============================================================
# 2. Paths
# ============================================================

OUT_DIR <- "results/03_FUS_KO"

INPUT <- file.path(
  OUT_DIR,
  "FUS_KO_QAPA_PAU.tsv"
)

XLSX_OUT <- file.path(
  OUT_DIR,
  "FUS_KO_differential_APA.xlsx"
)

VOLCANO_OUT <- file.path(
  OUT_DIR,
  "FUS_KO_volcano.pdf"
)

PCA_OUT <- file.path(
  OUT_DIR,
  "FUS_KO_PCA.pdf"
)

HEATMAP_OUT <- file.path(
  OUT_DIR,
  "FUS_KO_heatmap.pdf"
)

SL_OUT <- file.path(
  OUT_DIR,
  "FUS_KO_shortening_lengthening.pdf"
)


# ============================================================
# 3. Parameters
# ============================================================

TPM_CUTOFF <- 1
MIN_TPM_SAMPLES <- 2

DPAU_CUTOFF <- 10
FDR_CUTOFF <- 0.05

TOP_LABELS <- 10
TOP_HEATMAP_SHORT <- 20
TOP_HEATMAP_LONG <- 20

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

control_pau <- c(
  "WT_1.PAU",
  "WT_2.PAU",
  "WT_3.PAU",
  "WT_4.PAU"
)

ko_pau <- c(
  "A4_1.PAU",
  "A4_2.PAU",
  "A4_3.PAU",
  "A4_4.PAU",
  "A5_1.PAU",
  "A5_2.PAU",
  "A5_3.PAU",
  "A5_4.PAU"
)

control_tpm <- c(
  "WT_1.TPM",
  "WT_2.TPM",
  "WT_3.TPM",
  "WT_4.TPM"
)

ko_tpm <- c(
  "A4_1.TPM",
  "A4_2.TPM",
  "A4_3.TPM",
  "A4_4.TPM",
  "A5_1.TPM",
  "A5_2.TPM",
  "A5_3.TPM",
  "A5_4.TPM"
)


all_pau <- c(control_pau, ko_pau)
all_tpm <- c(control_tpm, ko_tpm)

# Validate expected sample columns before analysis
required_sample_cols <- c(all_pau, all_tpm)
missing_sample_cols <- setdiff(required_sample_cols, colnames(df))

if (length(missing_sample_cols) > 0) {
  message("")
  message("ERROR: The following expected PAU/TPM columns are missing:")
  message(paste(missing_sample_cols, collapse = ", "))
  message("")
  message("Available PAU/TPM columns in the input file:")
  message(
    paste(
      grep("\\.(PAU|TPM)$", colnames(df), value = TRUE),
      collapse = ", "
    )
  )
  stop("Sample column names do not match the FUS KO GSE71812 QAPA table.")
}


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
  df[, control_pau, drop = FALSE],
  na.rm = TRUE
)

df$FUS_KO_mean_PAU <- rowMeans(
  df[, ko_pau, drop = FALSE],
  na.rm = TRUE
)

df$DeltaPAU <- (
  df$FUS_KO_mean_PAU -
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
      as.numeric(df[i, ko_pau]),
      as.numeric(df[i, control_pau])
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
    title = "FUS KO differential alternative polyadenylation",
    subtitle = paste0(
      "FDR < ",
      FDR_CUTOFF,
      " and |DeltaPAU| >= ",
      DPAU_CUTOFF
    ),
    x = "DeltaPAU (FUS KO - WT)",
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
  "FUS KO"
)

pca_plot_df$Group <- factor(
  pca_plot_df$Group,
  levels = c(
    "WT",
    "FUS KO"
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
# 17. Heatmap will be generated after gene-level 3'UTR length
#     statistics are calculated.
# ============================================================


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

      WT_1 = calc_weighted_length(g, "WT_1.PAU"),
      WT_2 = calc_weighted_length(g, "WT_2.PAU"),
      WT_3 = calc_weighted_length(g, "WT_3.PAU"),
      WT_4 = calc_weighted_length(g, "WT_4.PAU"),

      A4_1 = calc_weighted_length(g, "A4_1.PAU"),
      A4_2 = calc_weighted_length(g, "A4_2.PAU"),
      A4_3 = calc_weighted_length(g, "A4_3.PAU"),
      A4_4 = calc_weighted_length(g, "A4_4.PAU"),
      A5_1 = calc_weighted_length(g, "A5_1.PAU"),
      A5_2 = calc_weighted_length(g, "A5_2.PAU"),
      A5_3 = calc_weighted_length(g, "A5_3.PAU"),
      A5_4 = calc_weighted_length(g, "A5_4.PAU"),

      stringsAsFactors = FALSE
    )
  }
)

length_df <- bind_rows(
  length_list
)


length_df$WT_mean_length <- rowMeans(
  length_df[, c(
    "WT_1", "WT_2", "WT_3", "WT_4"
  )],
  na.rm = TRUE
)

length_df$FUS_KO_mean_length <- rowMeans(
  length_df[, c(
    "A4_1", "A4_2", "A4_3", "A4_4",
    "A5_1", "A5_2", "A5_3", "A5_4"
  )],
  na.rm = TRUE
)

length_df$Delta_UTR_length <- (
  length_df$FUS_KO_mean_length -
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
            "A4_1", "A4_2", "A4_3", "A4_4",
            "A5_1", "A5_2", "A5_3", "A5_4"
          )
        ]
      ),
      as.numeric(
        length_df[
          i,
          c(
            "WT_1", "WT_2", "WT_3", "WT_4"
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
# 19. Gene-level 3'UTR shortening / lengthening scatter plot
#
# PANEL C FOR MAIN FIGURE
#
# Label:
#   Top 10 shortening genes
#   Top 10 lengthening genes
#
# Ranking:
#   lowest P value first, then largest absolute length change
#
# Colors:
#   Shortening = blue
#   Lengthening = red
#   Stable = gray
#
# Output:
#   Broader 10 x 6.5 inch publication PDF
# ============================================================

length_df$minus_log10_P <- -log10(
  pmax(
    length_df$P_value,
    1e-300
  )
)


# ------------------------------------------------------------
# Top 10 shortening genes
# ------------------------------------------------------------

top_short_genes <- length_df %>%

  filter(
    Direction == "Shortening",
    is.finite(P_value),
    is.finite(Delta_UTR_length),
    !is.na(Gene_Name),
    Gene_Name != ""
  ) %>%

  arrange(
    P_value,
    Delta_UTR_length
  ) %>%

  distinct(
    Gene_Name,
    .keep_all = TRUE
  ) %>%

  slice_head(
    n = 10
  )


# ------------------------------------------------------------
# Top 10 lengthening genes
# ------------------------------------------------------------

top_long_genes <- length_df %>%

  filter(
    Direction == "Lengthening",
    is.finite(P_value),
    is.finite(Delta_UTR_length),
    !is.na(Gene_Name),
    Gene_Name != ""
  ) %>%

  arrange(
    P_value,
    desc(Delta_UTR_length)
  ) %>%

  distinct(
    Gene_Name,
    .keep_all = TRUE
  ) %>%

  slice_head(
    n = 10
  )


message("")
message(
  "Panel C labels: ",
  nrow(top_short_genes),
  " shortening + ",
  nrow(top_long_genes),
  " lengthening"
)


message("")
message("Top 10 shortening genes:")

print(
  top_short_genes %>%
    select(
      Gene_Name,
      Delta_UTR_length,
      P_value,
      FDR
    )
)


message("")
message("Top 10 lengthening genes:")

print(
  top_long_genes %>%
    select(
      Gene_Name,
      Delta_UTR_length,
      P_value,
      FDR
    )
)


# ------------------------------------------------------------
# Broader publication-ready Panel C
# ------------------------------------------------------------

p_sl <- ggplot(
  length_df,
  aes(
    x = Delta_UTR_length / 1000,
    y = minus_log10_P
  )
) +

  # Stable genes
  geom_point(
    data = subset(
      length_df,
      Direction == "Stable"
    ),
    color = "grey75",
    size = 1.5,
    alpha = 0.50
  ) +

  # Shortening genes
  geom_point(
    data = subset(
      length_df,
      Direction == "Shortening"
    ),
    color = "#377EB8",
    size = 2.1,
    alpha = 0.88
  ) +

  # Lengthening genes
  geom_point(
    data = subset(
      length_df,
      Direction == "Lengthening"
    ),
    color = "#E41A1C",
    size = 2.1,
    alpha = 0.88
  ) +

  # +/-100 nt cutoffs
  geom_vline(
    xintercept = c(
      -LENGTH_CHANGE_CUTOFF / 1000,
       LENGTH_CHANGE_CUTOFF / 1000
    ),
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey35"
  ) +

  # Nominal P = 0.05 reference
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey35"
  ) +

  # Top 10 shortening labels
  geom_label_repel(
    data = top_short_genes,
    aes(
      label = Gene_Name
    ),
    fill = "#D5E9F3",
    color = "black",
    size = 3.0,
    label.size = 0.18,
    box.padding = 0.50,
    point.padding = 0.30,
    force = 5,
    force_pull = 0.8,
    segment.color = "grey40",
    segment.size = 0.30,
    min.segment.length = 0,
    max.overlaps = Inf,
    max.time = 10,
    max.iter = 50000,
    direction = "both",
    seed = 101
  ) +

  # Top 10 lengthening labels
  geom_label_repel(
    data = top_long_genes,
    aes(
      label = Gene_Name
    ),
    fill = "#F3D0CD",
    color = "black",
    size = 3.0,
    label.size = 0.18,
    box.padding = 0.50,
    point.padding = 0.30,
    force = 5,
    force_pull = 0.8,
    segment.color = "grey40",
    segment.size = 0.30,
    min.segment.length = 0,
    max.overlaps = Inf,
    max.time = 10,
    max.iter = 50000,
    direction = "both",
    seed = 202
  ) +

  scale_x_continuous(
    breaks = scales::pretty_breaks(
      n = 7
    ),
    expand = expansion(
      mult = c(
        0.10,
        0.10
      )
    )
  ) +

  scale_y_continuous(
    breaks = scales::pretty_breaks(
      n = 6
    ),
    expand = expansion(
      mult = c(
        0.02,
        0.08
      )
    )
  ) +

  labs(
    title = "FUS KO alters 3'UTR length",
    subtitle =
      "Gene-level PAU-weighted 3'UTR length; blue = shortening, red = lengthening",
    x = expression(
      Delta*"3'UTR length (FUS KO - WT, kb)"
    ),
    y = expression(
      -log[10]*"(P value)"
    )
  ) +

  theme_classic(
    base_size = 13
  ) +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 15,
      hjust = 0.5
    ),

    plot.subtitle = element_text(
      size = 10.5,
      hjust = 0.5,
      margin = margin(
        b = 7
      )
    ),

    axis.title = element_text(
      size = 12
    ),

    axis.text = element_text(
      size = 10.5,
      color = "black"
    ),

    plot.margin = margin(
      t = 10,
      r = 28,
      b = 10,
      l = 12
    ),

    legend.position = "none"
  )


# ------------------------------------------------------------
# Save broader Panel C PDF
# ------------------------------------------------------------

ggsave(
  SL_OUT,
  p_sl,
  width = 10,
  height = 6.5,
  units = "in",
  device = cairo_pdf
)


# Also save high-resolution PNG
ggsave(
  file.path(
    OUT_DIR,
    "FUS_KO_shortening_lengthening.png"
  ),
  p_sl,
  width = 10,
  height = 6.5,
  units = "in",
  dpi = 600
)


# ============================================================
# 19B. Heatmap:
#
#   Top 20 shortening genes
#   Top 20 lengthening genes
#
# Uses gene-level PAU-weighted 3'UTR length values.
# Ranking:
#   lowest P value first, then largest length change.
# ============================================================

heat_short <- length_df %>%
  filter(
    Direction == "Shortening",
    is.finite(P_value)
  ) %>%
  arrange(
    P_value,
    Delta_UTR_length
  ) %>%
  distinct(
    Gene_Name,
    .keep_all = TRUE
  ) %>%
  slice_head(
    n = TOP_HEATMAP_SHORT
  )

heat_long <- length_df %>%
  filter(
    Direction == "Lengthening",
    is.finite(P_value)
  ) %>%
  arrange(
    P_value,
    desc(Delta_UTR_length)
  ) %>%
  distinct(
    Gene_Name,
    .keep_all = TRUE
  ) %>%
  slice_head(
    n = TOP_HEATMAP_LONG
  )

heat_df <- bind_rows(
  heat_short,
  heat_long
)

message("")
message(
  "Heatmap genes: ",
  nrow(heat_short),
  " shortening + ",
  nrow(heat_long),
  " lengthening"
)

heat_mat <- as.matrix(
  heat_df[
    ,
    c(
      "WT_1", "WT_2", "WT_3", "WT_4",
      "A4_1", "A4_2", "A4_3", "A4_4",
      "A5_1", "A5_2", "A5_3", "A5_4"
    ),
    drop = FALSE
  ]
)

rownames(heat_mat) <- make.unique(
  paste0(
    heat_df$Gene_Name,
    " | ",
    heat_df$Direction
  )
)

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
    rep("WT", 4),
    rep("FUS KO", 8)
  )
)

rownames(annotation_col) <- colnames(
  heat_z
)


annotation_row <- data.frame(
  Direction = heat_df$Direction
)

rownames(annotation_row) <- rownames(
  heat_z
)


pdf(
  HEATMAP_OUT,
  width = 8.5,
  height = 10,
  useDingbats = FALSE
)

pheatmap(
  heat_z,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 7,
  fontsize_col = 9,
  border_color = NA,
  gaps_row = if (
    nrow(heat_short) > 0 &&
    nrow(heat_long) > 0
  ) {
    nrow(
      heat_short
    )
  } else {
    NULL
  },
  main = paste0(
    "Top ",
    nrow(
      heat_short
    ),
    " shortening and top ",
    nrow(
      heat_long
    ),
    " lengthening genes"
  )
)

dev.off()


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
