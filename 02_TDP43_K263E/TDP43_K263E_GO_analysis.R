library(dplyr)
library(tidyr)
library(gprofiler2)
library(ggplot2)
library(openxlsx)

# ============================================================
# Paths
# ============================================================

indir <- "results/02_TDP43_K263E"

infile <- file.path(
  indir,
  "TDP43_K263E_QAPA_PAU.tsv"
)

outdir <- file.path(indir, "GO_analysis")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Read QAPA
# ============================================================

df <- read.delim(
  infile,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Total QAPA events:", nrow(df), "\n")

# ============================================================
# Identify QAPA site type
# P = proximal
# D = distal
# S = single
# ============================================================

df$Site <- sub(".*_", "", df$APA_ID)

table(df$Site)

# ============================================================
# Calculate mean PAU
# ============================================================

df <- df %>%
  mutate(
    WT_mean = rowMeans(
      select(., WT_Rep1.PAU, WT_Rep2.PAU, WT_Rep3.PAU),
      na.rm = TRUE
    ),

    K263E_mean = rowMeans(
      select(., K263E_Rep1.PAU, K263E_Rep2.PAU, K263E_Rep3.PAU),
      na.rm = TRUE
    ),

    DeltaPAU = K263E_mean - WT_mean
  )

# ============================================================
# Event-level P values
# Welch t-test
# ============================================================

get_pvalue <- function(x) {

  wt <- as.numeric(x[c(
    "WT_Rep1.PAU",
    "WT_Rep2.PAU",
    "WT_Rep3.PAU"
  )])

  mut <- as.numeric(x[c(
    "K263E_Rep1.PAU",
    "K263E_Rep2.PAU",
    "K263E_Rep3.PAU"
  )])

  if (
    sum(is.finite(wt)) < 2 ||
    sum(is.finite(mut)) < 2
  ) {
    return(NA_real_)
  }

  tryCatch(
    t.test(mut, wt)$p.value,
    error = function(e) NA_real_
  )
}

df$P_value <- apply(df, 1, get_pvalue)
df$FDR <- p.adjust(df$P_value, method = "BH")

# ============================================================
# Significant mechanistic APA
#
# Shortening:
#   proximal ↑
#   distal ↓
#
# Lengthening:
#   proximal ↓
#   distal ↑
# ============================================================

df <- df %>%
  mutate(

    APA_direction = case_when(

      Site == "P" &
        DeltaPAU >= 10 &
        P_value < 0.05 ~ "Shortening",

      Site == "D" &
        DeltaPAU <= -10 &
        P_value < 0.05 ~ "Shortening",

      Site == "P" &
        DeltaPAU <= -10 &
        P_value < 0.05 ~ "Lengthening",

      Site == "D" &
        DeltaPAU >= 10 &
        P_value < 0.05 ~ "Lengthening",

      TRUE ~ "NS"
    )
  )

cat("\nAPA events:\n")
print(table(df$APA_direction))

# ============================================================
# Gene lists
# ============================================================

short_genes <- df %>%
  filter(APA_direction == "Shortening") %>%
  pull(Gene_Name) %>%
  unique() %>%
  na.omit()

long_genes <- df %>%
  filter(APA_direction == "Lengthening") %>%
  pull(Gene_Name) %>%
  unique() %>%
  na.omit()

# ============================================================
# Find genes appearing in BOTH categories
# ============================================================

discordant <- intersect(short_genes, long_genes)

cat("\nShortening genes before removing discordant:",
    length(short_genes), "\n")

cat("Lengthening genes before removing discordant:",
    length(long_genes), "\n")

cat("Discordant genes:",
    length(discordant), "\n")

# Remove conflicting genes
short_genes_clean <- setdiff(short_genes, discordant)
long_genes_clean  <- setdiff(long_genes, discordant)

cat("\nFinal shortening genes:",
    length(short_genes_clean), "\n")

cat("Final lengthening genes:",
    length(long_genes_clean), "\n")

# ============================================================
# QAPA background
#
# Only genes represented by multi-PAS QAPA events
# ============================================================

background <- df %>%
  filter(Site %in% c("P", "D")) %>%
  pull(Gene_Name) %>%
  unique() %>%
  na.omit()

cat("QAPA background genes:",
    length(background), "\n")

# ============================================================
# Save gene lists
# ============================================================

writeLines(
  short_genes_clean,
  file.path(outdir, "TDP43_K263E_shortening_genes.txt")
)

writeLines(
  long_genes_clean,
  file.path(outdir, "TDP43_K263E_lengthening_genes.txt")
)

writeLines(
  discordant,
  file.path(outdir, "TDP43_K263E_discordant_genes.txt")
)

writeLines(
  background,
  file.path(outdir, "TDP43_K263E_QAPA_background_genes.txt")
)

# ============================================================
# GO enrichment
# ============================================================

run_GO <- function(genes) {

  if (length(genes) < 2)
    return(NULL)

  gost(
    query = genes,
    organism = "hsapiens",
    ordered_query = FALSE,

    custom_bg = background,

    correction_method = "fdr",

    sources = c(
      "GO:BP",
      "GO:MF",
      "GO:CC"
    )
  )
}

GO_short <- run_GO(short_genes_clean)
GO_long  <- run_GO(long_genes_clean)

# ============================================================
# Extract results
# ============================================================

short_result <- if (!is.null(GO_short)) {
  GO_short$result
} else {
  data.frame()
}

long_result <- if (!is.null(GO_long)) {
  GO_long$result
} else {
  data.frame()
}

# ============================================================
# Excel
# ============================================================

wb <- createWorkbook()

addWorksheet(wb, "Shortening_GO")
writeData(wb, "Shortening_GO", short_result)

addWorksheet(wb, "Lengthening_GO")
writeData(wb, "Lengthening_GO", long_result)

addWorksheet(wb, "Shortening_genes")
writeData(
  wb,
  "Shortening_genes",
  data.frame(Gene = short_genes_clean)
)

addWorksheet(wb, "Lengthening_genes")
writeData(
  wb,
  "Lengthening_genes",
  data.frame(Gene = long_genes_clean)
)

addWorksheet(wb, "Discordant")
writeData(
  wb,
  "Discordant",
  data.frame(Gene = discordant)
)

saveWorkbook(
  wb,
  file.path(
    outdir,
    "TDP43_K263E_GO_enrichment.xlsx"
  ),
  overwrite = TRUE
)

# ============================================================
# GO plotting function
# ============================================================

plot_GO <- function(res, direction, source) {

  x <- res %>%
    filter(source == !!source) %>%
    arrange(p_value) %>%
    slice_head(n = 15)

  if (nrow(x) == 0)
    return(NULL)

  x <- x %>%
    mutate(
      term_name = factor(
        term_name,
        levels = rev(term_name)
      ),
      logP = -log10(p_value)
    )

  ggplot(
    x,
    aes(
      x = logP,
      y = term_name,
      fill = logP
    )
  ) +

    geom_col(width = 0.75) +

    scale_fill_viridis_c(
      option = "viridis",
      direction = -1
    ) +

    labs(
      title = paste(
        "TDP-43 K263E 3'UTR",
        direction
      ),
      subtitle = source,
      x = expression(-log[10](P)),
      y = NULL
    ) +

    theme_bw(base_size = 13) +

    theme(
      plot.title =
        element_text(
          face = "bold",
          size = 16
        ),

      plot.subtitle =
        element_text(
          size = 13
        ),

      legend.position = "none",

      panel.grid.major.y =
        element_blank()
    )
}

# ============================================================
# Generate PDFs
# ============================================================

for (src in c("GO:BP", "GO:MF", "GO:CC")) {

  p <- plot_GO(
    short_result,
    "Shortening",
    src
  )

  if (!is.null(p)) {

    ggsave(
      file.path(
        outdir,
        paste0(
          "TDP43_K263E_shortening_",
          gsub(":", "_", src),
          ".pdf"
        )
      ),
      p,
      width = 8,
      height = 6
    )
  }


  p <- plot_GO(
    long_result,
    "Lengthening",
    src
  )

  if (!is.null(p)) {

    ggsave(
      file.path(
        outdir,
        paste0(
          "TDP43_K263E_lengthening_",
          gsub(":", "_", src),
          ".pdf"
        )
      ),
      p,
      width = 8,
      height = 6
    )
  }
}

cat("\n====================================\n")
cat("GO ANALYSIS COMPLETE\n")
cat("====================================\n")

cat("Output directory:\n")
cat(outdir, "\n")
