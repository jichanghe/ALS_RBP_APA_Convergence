# ============================================================
# FUS P525L QAPA
# Mechanistic APA-site volcano
#
# Four mechanistic categories:
#
# Proximal PAS increased  -> Shortening
# Proximal PAS decreased  -> Lengthening
# Distal PAS increased    -> Lengthening
# Distal PAS decreased    -> Shortening
#
# x = Delta PAU (FUS P525L - WT)
# y = -log10(P value)
#
# Shape:
#   circle   = Proximal
#   triangle = Distal
#
# Color:
#   blue = Shortening
#   red  = Lengthening
#   grey = NS
#
# Significance:
#   nominal P < 0.05
#   |Delta PAU| >= 10
#
# ============================================================


# ============================================================
# 1. Packages
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})


# ============================================================
# 2. Paths
# ============================================================

OUT_DIR <- "results/04_FUS_P525L"

INPUT <- file.path(
  OUT_DIR,
  "FUS_P525L_QAPA_PAU.tsv"
)

TABLE_OUT <- file.path(
  OUT_DIR,
  "FUS_P525L_APA_site_mechanistic_results.tsv"
)

PDF_OUT <- file.path(
  OUT_DIR,
  "FUS_P525L_APA_site_mechanistic_volcano.pdf"
)

PNG_OUT <- file.path(
  OUT_DIR,
  "FUS_P525L_APA_site_mechanistic_volcano.png"
)


# ============================================================
# 3. Parameters
# ============================================================

TPM_CUTOFF <- 1

MIN_TPM_SAMPLES <- 2

P_CUTOFF <- 0.05

DPAU_CUTOFF <- 10

TOP_N <- 10


# ============================================================
# 4. Read QAPA
# ============================================================

cat("\nReading QAPA PAU table...\n")

df <- read.delim(
  INPUT,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Total APA events:", nrow(df), "\n")
cat("Unique genes:", length(unique(df$Gene)), "\n")


# ============================================================
# 5. Sample columns
# ============================================================

WT_PAU <- c(
  "WT_Rep1.PAU",
  "WT_Rep2.PAU",
  "WT_Rep3.PAU"
)

MUT_PAU <- c(
  "P525L_Rep1.PAU",
  "P525L_Rep2.PAU",
  "P525L_Rep3.PAU"
)

WT_TPM <- c(
  "WT_Rep1.TPM",
  "WT_Rep2.TPM",
  "WT_Rep3.TPM"
)

MUT_TPM <- c(
  "P525L_Rep1.TPM",
  "P525L_Rep2.TPM",
  "P525L_Rep3.TPM"
)

ALL_PAU <- c(
  WT_PAU,
  MUT_PAU
)

ALL_TPM <- c(
  WT_TPM,
  MUT_TPM
)


# ============================================================
# 5B. Validate FUS P525L sample columns
# ============================================================

required_cols <- c(
  ALL_PAU,
  ALL_TPM
)

missing_cols <- setdiff(
  required_cols,
  colnames(df)
)

if (length(missing_cols) > 0) {

  cat(
    "\nERROR: Missing expected FUS P525L QAPA columns:\n"
  )

  print(
    missing_cols
  )

  cat(
    "\nAvailable PAU/TPM columns:\n"
  )

  print(
    grep(
      "\\.(PAU|TPM)$",
      colnames(df),
      value = TRUE
    )
  )

  stop(
    "Sample columns do not match FUS_P525L_Validation."
  )
}


# ============================================================
# 6. Convert numeric columns
# ============================================================

numeric_cols <- c(
  ALL_PAU,
  ALL_TPM,
  "Length",
  "Num_Events",
  "UTR3.Start",
  "UTR3.End"
)

for (x in numeric_cols) {

  if (x %in% colnames(df)) {

    df[[x]] <- suppressWarnings(
      as.numeric(df[[x]])
    )
  }
}


# ============================================================
# 7. Expression filter
#
# TPM >= 1 in at least 2/6 samples
# ============================================================

df$TPM_pass_n <- rowSums(
  df[, ALL_TPM, drop = FALSE] >= TPM_CUTOFF,
  na.rm = TRUE
)

df$Expression_Pass <- (
  df$TPM_pass_n >= MIN_TPM_SAMPLES
)


# ============================================================
# 8. Keep genes with >=2 APA events
# ============================================================

dat <- df %>%
  filter(
    Expression_Pass,
    Num_Events >= 2
  )

cat(
  "APA events after filtering:",
  nrow(dat),
  "\n"
)


# ============================================================
# 9. Mean PAU
# ============================================================

dat$WT_mean_PAU <- rowMeans(
  dat[, WT_PAU, drop = FALSE],
  na.rm = TRUE
)

dat$FUS_P525L_mean_PAU <- rowMeans(
  dat[, MUT_PAU, drop = FALSE],
  na.rm = TRUE
)


# ============================================================
# 10. Delta PAU
#
# positive:
# increased usage in FUS P525L
#
# negative:
# decreased usage in FUS P525L
# ============================================================

dat$DeltaPAU <- (
  dat$FUS_P525L_mean_PAU -
  dat$WT_mean_PAU
)


# ============================================================
# 11. Event-level statistical test
# ============================================================

safe_ttest <- function(mut, wt) {

  mut <- as.numeric(mut)
  wt  <- as.numeric(wt)

  mut <- mut[
    is.finite(mut)
  ]

  wt <- wt[
    is.finite(wt)
  ]

  if (
    length(mut) < 2 ||
    length(wt) < 2
  ) {

    return(
      NA_real_
    )
  }

  if (
    length(unique(mut)) == 1 &&
    length(unique(wt)) == 1
  ) {

    if (
      unique(mut) ==
      unique(wt)
    ) {

      return(1)

    }
  }

  tryCatch(

    t.test(
      mut,
      wt,
      paired = FALSE,
      var.equal = FALSE
    )$p.value,

    error = function(e) {

      NA_real_
    }
  )
}


dat$P_value <- vapply(

  seq_len(
    nrow(dat)
  ),

  function(i) {

    safe_ttest(

      dat[
        i,
        MUT_PAU
      ],

      dat[
        i,
        WT_PAU
      ]
    )
  },

  numeric(1)
)


# ============================================================
# 12. FDR
# ============================================================

dat$FDR <- p.adjust(
  dat$P_value,
  method = "BH"
)


# ============================================================
# 13. -log10(P)
# ============================================================

dat$minus_log10_P <- -log10(
  pmax(
    dat$P_value,
    1e-300
  )
)


# ============================================================
# 14. Identify PROXIMAL / DISTAL APA sites
#
# IMPORTANT:
#
# QAPA Length represents the 3'UTR length associated with
# each APA event.
#
# Within each gene:
#
# shortest 3'UTR = Proximal
# longest  3'UTR = Distal
#
# Intermediate APA sites are classified separately and are
# not used for the four-category mechanistic call.
# ============================================================

dat <- dat %>%
  group_by(Gene) %>%

  mutate(

    Min_Length = min(
      Length,
      na.rm = TRUE
    ),

    Max_Length = max(
      Length,
      na.rm = TRUE
    ),

    APA_Site = case_when(

      Length == Min_Length &
      Length < Max_Length ~
        "Proximal",

      Length == Max_Length &
      Length > Min_Length ~
        "Distal",

      TRUE ~
        "Intermediate"
    )
  ) %>%

  ungroup()


# ============================================================
# 15. Check proximal/distal numbers
# ============================================================

cat(
  "\nAPA site classification:\n"
)

print(
  table(
    dat$APA_Site
  )
)


# ============================================================
# 16. Four mechanistic categories
#
#                       DeltaPAU
#
#                   -             +
#
# Proximal     Lengthening    Shortening
#
# Distal       Shortening     Lengthening
#
# ============================================================

dat$Mechanism <- "NS"


# ------------------------------------------------------------
# Proximal increased
# -> more short isoform
# -> shortening
# ------------------------------------------------------------

dat$Mechanism[
  dat$APA_Site == "Proximal" &
  dat$DeltaPAU >= DPAU_CUTOFF &
  dat$P_value < P_CUTOFF
] <- "Shortening"


# ------------------------------------------------------------
# Proximal decreased
# -> less short isoform
# -> lengthening
# ------------------------------------------------------------

dat$Mechanism[
  dat$APA_Site == "Proximal" &
  dat$DeltaPAU <= -DPAU_CUTOFF &
  dat$P_value < P_CUTOFF
] <- "Lengthening"


# ------------------------------------------------------------
# Distal increased
# -> more long isoform
# -> lengthening
# ------------------------------------------------------------

dat$Mechanism[
  dat$APA_Site == "Distal" &
  dat$DeltaPAU >= DPAU_CUTOFF &
  dat$P_value < P_CUTOFF
] <- "Lengthening"


# ------------------------------------------------------------
# Distal decreased
# -> less long isoform
# -> shortening
# ------------------------------------------------------------

dat$Mechanism[
  dat$APA_Site == "Distal" &
  dat$DeltaPAU <= -DPAU_CUTOFF &
  dat$P_value < P_CUTOFF
] <- "Shortening"


# ============================================================
# 17. Explicit four-category labels
# ============================================================

dat$Four_Category <- "NS"


dat$Four_Category[
  dat$APA_Site == "Proximal" &
  dat$DeltaPAU >= DPAU_CUTOFF &
  dat$P_value < P_CUTOFF
] <- "Proximal increased"


dat$Four_Category[
  dat$APA_Site == "Proximal" &
  dat$DeltaPAU <= -DPAU_CUTOFF &
  dat$P_value < P_CUTOFF
] <- "Proximal decreased"


dat$Four_Category[
  dat$APA_Site == "Distal" &
  dat$DeltaPAU >= DPAU_CUTOFF &
  dat$P_value < P_CUTOFF
] <- "Distal increased"


dat$Four_Category[
  dat$APA_Site == "Distal" &
  dat$DeltaPAU <= -DPAU_CUTOFF &
  dat$P_value < P_CUTOFF
] <- "Distal decreased"


# ============================================================
# 18. Summary
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "FOUR MECHANISTIC APA CATEGORIES\n"
)

cat(
  "========================================\n"
)

print(
  table(
    dat$Four_Category
  )
)


cat(
  "\nBiological direction:\n"
)

print(
  table(
    dat$Mechanism
  )
)


# ============================================================
# 19. Plot only proximal/distal sites
#
# Intermediate sites are excluded from the volcano because
# they cannot be interpreted using the simple four-category
# proximal/distal model.
# ============================================================

plot_df <- dat %>%
  filter(
    APA_Site %in% c(
      "Proximal",
      "Distal"
    ),
    is.finite(DeltaPAU),
    is.finite(P_value)
  )


# ============================================================
# 20. Factors
# ============================================================

plot_df$APA_Site <- factor(
  plot_df$APA_Site,
  levels = c(
    "Proximal",
    "Distal"
  )
)

plot_df$Mechanism <- factor(
  plot_df$Mechanism,
  levels = c(
    "Shortening",
    "Lengthening",
    "NS"
  )
)


# ============================================================
# 21. Select TOP 10 SHORTENING genes
#
# One label per gene.
# Prioritize lowest P value.
# ============================================================

top_short <- plot_df %>%

  filter(
    Mechanism == "Shortening"
  ) %>%

  arrange(
    P_value,
    desc(abs(DeltaPAU))
  ) %>%

  distinct(
    Gene_Name,
    .keep_all = TRUE
  ) %>%

  slice_head(
    n = TOP_N
  )


# ============================================================
# 22. Select TOP 10 LENGTHENING genes
# ============================================================

top_long <- plot_df %>%

  filter(
    Mechanism == "Lengthening"
  ) %>%

  arrange(
    P_value,
    desc(abs(DeltaPAU))
  ) %>%

  distinct(
    Gene_Name,
    .keep_all = TRUE
  ) %>%

  slice_head(
    n = TOP_N
  )


# ============================================================
# 23. Print top genes
# ============================================================

cat(
  "\nTop 10 shortening:\n"
)

print(
  top_short %>%
    select(
      Gene_Name,
      APA_Site,
      DeltaPAU,
      P_value,
      FDR,
      Four_Category
    )
)


cat(
  "\nTop 10 lengthening:\n"
)

print(
  top_long %>%
    select(
      Gene_Name,
      APA_Site,
      DeltaPAU,
      P_value,
      FDR,
      Four_Category
    )
)


# ============================================================
# 24. Save complete results
# ============================================================

write.table(
  dat,
  TABLE_OUT,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 25. Volcano
# ============================================================

p <- ggplot(
  plot_df,
  aes(
    x = DeltaPAU,
    y = minus_log10_P
  )
) +

  # ----------------------------------------------------------
  # Neutral |DeltaPAU| < 10 region
  # ----------------------------------------------------------

  annotate(
    "rect",

    xmin = -DPAU_CUTOFF,
    xmax = DPAU_CUTOFF,

    ymin = -Inf,
    ymax = Inf,

    fill = "grey94"
  ) +

  # ----------------------------------------------------------
  # NS
  # ----------------------------------------------------------

  geom_point(

    data = subset(
      plot_df,
      Mechanism == "NS"
    ),

    aes(
      shape = APA_Site
    ),

    color = "grey78",

    size = 1.6,

    alpha = 0.55
  ) +

  # ----------------------------------------------------------
  # Shortening
  # ----------------------------------------------------------

  geom_point(

    data = subset(
      plot_df,
      Mechanism == "Shortening"
    ),

    aes(
      shape = APA_Site
    ),

    color = "#2878B5",

    size = 2.8,

    alpha = 0.95
  ) +

  # ----------------------------------------------------------
  # Lengthening
  # ----------------------------------------------------------

  geom_point(

    data = subset(
      plot_df,
      Mechanism == "Lengthening"
    ),

    aes(
      shape = APA_Site
    ),

    color = "#D63B32",

    size = 2.8,

    alpha = 0.95
  ) +

  # ----------------------------------------------------------
  # DeltaPAU = 0
  # ----------------------------------------------------------

  geom_vline(
    xintercept = 0,
    color = "grey45",
    linewidth = 0.45
  ) +

  # ----------------------------------------------------------
  # DeltaPAU +/-10
  # ----------------------------------------------------------

  geom_vline(
    xintercept = c(
      -DPAU_CUTOFF,
      DPAU_CUTOFF
    ),
    linetype = "dashed",
    linewidth = 0.5
  ) +

  # ----------------------------------------------------------
  # P = 0.05
  # ----------------------------------------------------------

  geom_hline(
    yintercept = -log10(P_CUTOFF),
    linetype = "dashed",
    linewidth = 0.5
  ) +

  # ==========================================================
  # Top shortening labels
  # ==========================================================

  geom_label_repel(

    data = top_short,

    aes(
      label = Gene_Name
    ),

    fill = "#D5E9F3",

    color = "black",

    size = 3,

    label.size = 0.18,

    box.padding = 0.5,

    point.padding = 0.25,

    force = 5,

    segment.color = "grey40",

    segment.size = 0.3,

    min.segment.length = 0,

    max.overlaps = Inf,

    seed = 100
  ) +

  # ==========================================================
  # Top lengthening labels
  # ==========================================================

  geom_label_repel(

    data = top_long,

    aes(
      label = Gene_Name
    ),

    fill = "#F3D0CD",

    color = "black",

    size = 3,

    label.size = 0.18,

    box.padding = 0.5,

    point.padding = 0.25,

    force = 5,

    segment.color = "grey40",

    segment.size = 0.3,

    min.segment.length = 0,

    max.overlaps = Inf,

    seed = 200
  ) +

  # ==========================================================
  # Shapes
  # ==========================================================

  scale_shape_manual(
    name = "APA site",
    values = c(
      "Proximal" = 16,
      "Distal" = 17
    )
  ) +

  # ==========================================================
  # X axis
  # ==========================================================

  scale_x_continuous(

    limits = c(
      -100,
      100
    ),

    oob = scales::squish,

    breaks = c(
      -100,
      -75,
      -50,
      -25,
      0,
      25,
      50,
      75,
      100
    ),

    expand = expansion(
      mult = c(
        0.02,
        0.02
      )
    )
  ) +

  # ==========================================================
  # Labels
  # ==========================================================

  labs(

    title =
      "FUS P525L alters alternative polyadenylation",

    subtitle =
      "Proximal and distal PAS usage reveals 3'UTR shortening and lengthening",

    x =
      expression(
        Delta*"PAU (FUS P525L - WT)"
      ),

    y =
      expression(
        -log[10]*"(P value)"
      )
  ) +

  # ==========================================================
  # Theme
  # ==========================================================

  theme_classic(
    base_size = 14
  ) +

  theme(

    plot.title = element_text(
      face = "bold",
      size = 17,
      hjust = 0.5
    ),

    plot.subtitle = element_text(
      size = 11,
      hjust = 0.5,
      margin = margin(
        b = 10
      )
    ),

    axis.title = element_text(
      size = 14
    ),

    axis.text = element_text(
      size = 11,
      color = "black"
    ),

    legend.title = element_text(
      face = "bold"
    ),

    legend.position = "right",

    plot.margin = margin(
      12,
      15,
      12,
      12
    )
  )


# ============================================================
# 26. Save
# ============================================================

ggsave(
  PDF_OUT,
  p,
  width = 10,
  height = 6.5,
  device = cairo_pdf
)

ggsave(
  PNG_OUT,
  p,
  width = 10,
  height = 6.5,
  dpi = 600
)


# ============================================================
# 27. Finished
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "DONE\n"
)

cat(
  "========================================\n"
)

cat(
  "\nResults:\n",
  TABLE_OUT,
  "\n"
)

cat(
  "\nPDF:\n",
  PDF_OUT,
  "\n"
)

cat(
  "\nPNG:\n",
  PNG_OUT,
  "\n"
)
