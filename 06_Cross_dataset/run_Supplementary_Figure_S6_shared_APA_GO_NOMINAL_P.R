# ============================================================
# Supplementary Figure S6
# GO enrichment of APA-regulated genes shared across >=2 ALS-RBP models
#
# IMPORTANT:
# This version uses the SAME exploratory APA threshold used for the
# cross-model comparisons / mechanistic APA plots:
#
#       P < 0.05 AND |DeltaPAU| >= 10
#
# It does NOT require FDR < 0.05, because several models have zero
# FDR-significant APA events.
#
# Models:
#   TDP-43 KD
#   TDP-43 K263E
#   FUS KO
#   FUS P525L
#   MATR3 KO
#
# Output:
#   Shared_APA_genes_nominalP_ge2_models.tsv
#   Shared_APA_ge2models_GO_enrichment_nominalP.xlsx
#   Supplementary_Figure_S6_shared_APA_GO_nominalP.pdf
#   Supplementary_Figure_S6_shared_APA_GO_nominalP.png
# ============================================================


# ============================================================
# 1. Packages
# ============================================================

cran_packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "openxlsx",
  "patchwork"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

for (pkg in c("clusterProfiler", "org.Hs.eg.db")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
  library(patchwork)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})


# ============================================================
# 2. Parameters
# ============================================================

P_CUTOFF <- 0.05
DPAU_CUTOFF <- 10
MIN_MODELS <- 2
TOP_N <- 15


# ============================================================
# 3. Output directory
# ============================================================

OUT_DIR <- "results/CROSS_MODEL_SHARED_GO"

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

SHARED_TSV <- file.path(
  OUT_DIR,
  "Shared_APA_genes_nominalP_ge2_models.tsv"
)

GO_XLSX <- file.path(
  OUT_DIR,
  "Shared_APA_ge2models_GO_enrichment_nominalP.xlsx"
)

GO_PDF <- file.path(
  OUT_DIR,
  "Supplementary_Figure_S6_shared_APA_GO_nominalP.pdf"
)

GO_PNG <- file.path(
  OUT_DIR,
  "Supplementary_Figure_S6_shared_APA_GO_nominalP.png"
)


# ============================================================
# 4. Exact workbook paths
# ============================================================

workbooks <- c(

  "TDP-43 KD" =
    "data/TDP43_KD_differential_APA.xlsx",

  "TDP-43 K263E" =
    "data/TDP43_K263E_differential_APA.xlsx",

  "FUS KO" =
    "data/FUS_KO_differential_APA.xlsx",

  "FUS P525L" =
    "data/FUS_P525L_differential_APA.xlsx",

  "MATR3 KO" =
    "data/MATR3_KO_differential_APA.xlsx"
)


for (x in workbooks) {
  if (!file.exists(x)) {
    stop("Missing workbook: ", x)
  }
}


# ============================================================
# 5. Helper: gene-name column
# ============================================================

find_gene_col <- function(dat) {

  candidates <- c(
    "Gene_Name",
    "gene_name",
    "GeneName",
    "SYMBOL",
    "Symbol",
    "symbol"
  )

  hit <- candidates[
    candidates %in% colnames(dat)
  ]

  if (length(hit) == 0) {
    stop(
      "No recognizable gene-name column. Available columns: ",
      paste(colnames(dat), collapse = ", ")
    )
  }

  hit[1]
}


# ============================================================
# 6. Read All_APA_events and apply SAME threshold to every model
# ============================================================

read_model_APA <- function(xlsx, model_name) {

  sheets <- getSheetNames(xlsx)

  if (!"All_APA_events" %in% sheets) {
    stop(
      "Workbook does not contain All_APA_events: ",
      xlsx
    )
  }

  dat <- read.xlsx(
    xlsx,
    sheet = "All_APA_events"
  )

  gene_col <- find_gene_col(dat)

  required <- c(
    gene_col,
    "DeltaPAU",
    "P_value"
  )

  missing_cols <- setdiff(
    required,
    colnames(dat)
  )

  if (length(missing_cols) > 0) {
    stop(
      model_name,
      " is missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  dat[[gene_col]] <- trimws(
    as.character(
      dat[[gene_col]]
    )
  )

  dat$DeltaPAU <- suppressWarnings(
    as.numeric(
      dat$DeltaPAU
    )
  )

  dat$P_value <- suppressWarnings(
    as.numeric(
      dat$P_value
    )
  )

  # Expression/tested background
  bg <- dat

  if ("Expression_Pass" %in% colnames(bg)) {

    pass <- tolower(
      as.character(
        bg$Expression_Pass
      )
    )

    bg <- bg[
      pass %in% c(
        "true",
        "t",
        "1"
      ),
      ,
      drop = FALSE
    ]
  }

  background_genes <- unique(
    bg[[gene_col]]
  )

  background_genes <- background_genes[
    !is.na(background_genes) &
    background_genes != "" &
    background_genes != "NA"
  ]

  # Unified exploratory APA significance
  sig <- dat %>%

    filter(
      is.finite(P_value),
      is.finite(DeltaPAU),
      P_value < P_CUTOFF,
      abs(DeltaPAU) >= DPAU_CUTOFF
    )

  # If Expression_Pass exists, require it
  if ("Expression_Pass" %in% colnames(sig)) {

    pass <- tolower(
      as.character(
        sig$Expression_Pass
      )
    )

    sig <- sig[
      pass %in% c(
        "true",
        "t",
        "1"
      ),
      ,
      drop = FALSE
    ]
  }

  sig_genes <- unique(
    sig[[gene_col]]
  )

  sig_genes <- sig_genes[
    !is.na(sig_genes) &
    sig_genes != "" &
    sig_genes != "NA"
  ]

  list(
    significant_genes = sig_genes,
    background_genes = background_genes,
    significant_events = sig
  )
}


# ============================================================
# 7. Read all models
# ============================================================

model_results <- lapply(
  names(workbooks),
  function(model_name) {

    message("")
    message("Reading: ", model_name)
    message(workbooks[[model_name]])

    read_model_APA(
      workbooks[[model_name]],
      model_name
    )
  }
)

names(model_results) <- names(workbooks)


# ============================================================
# 8. Print model counts
# ============================================================

message("")
message("==============================================")
message("Unified APA criterion:")
message("P < ", P_CUTOFF, " and |DeltaPAU| >= ", DPAU_CUTOFF)
message("==============================================")

for (model in names(model_results)) {

  message(
    model,
    ": ",
    length(
      model_results[[model]]$significant_genes
    ),
    " genes; ",
    nrow(
      model_results[[model]]$significant_events
    ),
    " APA events"
  )
}


# ============================================================
# 9. Build gene x model membership table
# ============================================================

membership <- bind_rows(
  lapply(
    names(model_results),
    function(model) {

      genes <- model_results[[model]]$significant_genes

      if (length(genes) == 0) {
        return(
          data.frame(
            Gene = character(0),
            Model = character(0),
            stringsAsFactors = FALSE
          )
        )
      }

      data.frame(
        Gene = genes,
        Model = rep(
          model,
          length(genes)
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)


if (nrow(membership) == 0) {
  stop(
    "No significant APA-regulated genes found using unified threshold."
  )
}


shared_table <- membership %>%

  distinct(
    Gene,
    Model
  ) %>%

  group_by(
    Gene
  ) %>%

  summarise(
    N_models = n(),
    Models = paste(
      sort(
        unique(Model)
      ),
      collapse = "; "
    ),
    .groups = "drop"
  ) %>%

  arrange(
    desc(N_models),
    Gene
  )


shared_ge2 <- shared_table %>%
  filter(
    N_models >= MIN_MODELS
  )


message("")
message(
  "Genes significant in >= ",
  MIN_MODELS,
  " models: ",
  nrow(shared_ge2)
)

message("")
message("Distribution of number of models per gene:")
print(
  table(
    shared_table$N_models
  )
)


write.table(
  shared_ge2,
  SHARED_TSV,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 10. Common tested background across all 5 models
# ============================================================

background_lists <- lapply(
  model_results,
  function(x) {
    x$background_genes
  }
)

common_background <- Reduce(
  intersect,
  background_lists
)


message("")
message(
  "Common tested-gene background: ",
  length(common_background)
)


if (length(common_background) < 1000) {

  warning(
    "Common background is relatively small: ",
    length(common_background),
    " genes."
  )
}


# ============================================================
# 11. Shared genes used for GO
# ============================================================

shared_genes <- intersect(
  shared_ge2$Gene,
  common_background
)

message(
  "Shared >=2 genes present in common background: ",
  length(shared_genes)
)


if (length(shared_genes) < 10) {

  stop(
    "Only ",
    length(shared_genes),
    " shared genes remain in common background. ",
    "Too few for reliable GO analysis."
  )
}


# ============================================================
# 12. Map SYMBOL -> ENTREZID
# ============================================================

shared_map <- bitr(
  shared_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)


background_map <- bitr(
  common_background,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)


shared_entrez <- unique(
  shared_map$ENTREZID
)

background_entrez <- unique(
  background_map$ENTREZID
)


message("")
message(
  "Mapped shared genes: ",
  length(shared_entrez)
)

message(
  "Mapped background genes: ",
  length(background_entrez)
)


# ============================================================
# 13. GO enrichment
# ============================================================

run_go <- function(ontology) {

  enrichGO(
    gene = shared_entrez,
    universe = background_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )
}


bp_df <- as.data.frame(
  run_go("BP")
)

cc_df <- as.data.frame(
  run_go("CC")
)

mf_df <- as.data.frame(
  run_go("MF")
)


# ============================================================
# 14. Add numeric GeneRatio
# ============================================================

ratio_to_numeric <- function(x) {

  vapply(
    strsplit(
      as.character(x),
      "/",
      fixed = TRUE
    ),
    function(z) {

      if (length(z) != 2) {
        return(NA_real_)
      }

      as.numeric(z[1]) /
        as.numeric(z[2])
    },
    numeric(1)
  )
}


prepare_go <- function(dat) {

  if (nrow(dat) == 0) {
    return(dat)
  }

  dat$GeneRatio_numeric <- ratio_to_numeric(
    dat$GeneRatio
  )

  dat
}


bp_df <- prepare_go(bp_df)
cc_df <- prepare_go(cc_df)
mf_df <- prepare_go(mf_df)


# ============================================================
# 15. Select top terms
# ============================================================

select_top <- function(dat) {

  if (nrow(dat) == 0) {
    return(dat)
  }

  dat %>%
    arrange(
      p.adjust,
      pvalue,
      desc(Count)
    ) %>%
    slice_head(
      n = min(
        TOP_N,
        nrow(dat)
      )
    )
}


bp_plot <- select_top(bp_df)
cc_plot <- select_top(cc_df)
mf_plot <- select_top(mf_df)


# ============================================================
# 16. Publication-style GO plot
# ============================================================

make_go_plot <- function(dat, title_text) {

  if (nrow(dat) == 0) {

    return(
      ggplot() +
        annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = "No GO terms returned",
          size = 5
        ) +
        xlim(0, 1) +
        ylim(0, 1) +
        labs(
          title = title_text
        ) +
        theme_void()
    )
  }

  dat <- dat %>%

    mutate(
      Description = factor(
        Description,
        levels = rev(
          Description
        )
      )
    )


  n_sig <- sum(
    dat$p.adjust < 0.05,
    na.rm = TRUE
  )


  subtitle_text <- if (n_sig > 0) {

    paste0(
      n_sig,
      " displayed GO terms significant at FDR < 0.05"
    )

  } else {

    "Top exploratory terms; none significant after FDR correction"
  }


  ggplot(
    dat,
    aes(
      x = GeneRatio_numeric,
      y = Description
    )
  ) +

    geom_segment(
      aes(
        x = 0,
        xend = GeneRatio_numeric,
        y = Description,
        yend = Description
      ),
      linewidth = 0.4,
      color = "grey65"
    ) +

    geom_point(
      aes(
        size = Count,
        color = p.adjust
      ),
      alpha = 0.95
    ) +

    scale_color_viridis_c(
      option = "D",
      direction = -1,
      name = "FDR adjusted\nP value"
    ) +

    scale_size_continuous(
      name = "Gene count",
      range = c(2.5, 7)
    ) +

    scale_x_continuous(
      expand = expansion(
        mult = c(
          0,
          0.08
        )
      )
    ) +

    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "GeneRatio",
      y = NULL
    ) +

    theme_bw(
      base_size = 9
    ) +

    theme(
      plot.title = element_text(
        face = "bold",
        size = 11,
        hjust = 0.5
      ),

      plot.subtitle = element_text(
        size = 8,
        hjust = 0.5
      ),

      axis.text.y = element_text(
        size = 7.5,
        color = "black"
      ),

      axis.text.x = element_text(
        size = 8,
        color = "black"
      ),

      panel.grid.minor = element_blank()
    )
}


p_bp <- make_go_plot(
  bp_plot,
  "GO Biological Process"
)

p_cc <- make_go_plot(
  cc_plot,
  "GO Cellular Component"
)

p_mf <- make_go_plot(
  mf_plot,
  "GO Molecular Function"
)


# ============================================================
# 17. Combine S6
# ============================================================

combined <- (
  p_bp /
  p_cc /
  p_mf
) +

  plot_annotation(

    title =
      "GO enrichment of APA-regulated genes shared across >=2 ALS-RBP models",

    subtitle = paste0(
      "APA criterion: P < ",
      P_CUTOFF,
      " and |DeltaPAU| >= ",
      DPAU_CUTOFF,
      "; shared genes = ",
      length(shared_genes),
      "; common background = ",
      length(common_background)
    ),

    theme = theme(

      plot.title = element_text(
        face = "bold",
        size = 14,
        hjust = 0.5
      ),

      plot.subtitle = element_text(
        size = 10,
        hjust = 0.5
      )
    )
  )


# ============================================================
# 18. Save PDF + PNG
# ============================================================

ggsave(
  GO_PDF,
  combined,
  width = 11,
  height = 15,
  units = "in",
  device = cairo_pdf
)


ggsave(
  GO_PNG,
  combined,
  width = 11,
  height = 15,
  units = "in",
  dpi = 600
)


# ============================================================
# 19. Excel workbook
# ============================================================

wb <- createWorkbook()


addWorksheet(
  wb,
  "Shared_ge2_models"
)

writeData(
  wb,
  "Shared_ge2_models",
  shared_ge2
)


addWorksheet(
  wb,
  "Common_background"
)

writeData(
  wb,
  "Common_background",
  data.frame(
    Gene = common_background
  )
)


addWorksheet(
  wb,
  "GO_BP"
)

writeData(
  wb,
  "GO_BP",
  bp_df
)


addWorksheet(
  wb,
  "GO_CC"
)

writeData(
  wb,
  "GO_CC",
  cc_df
)


addWorksheet(
  wb,
  "GO_MF"
)

writeData(
  wb,
  "GO_MF",
  mf_df
)


model_summary <- bind_rows(
  lapply(
    names(model_results),
    function(model) {

      data.frame(
        Model = model,

        Significant_APA_events =
          nrow(
            model_results[[model]]$significant_events
          ),

        Significant_APA_genes =
          length(
            model_results[[model]]$significant_genes
          ),

        stringsAsFactors = FALSE
      )
    }
  )
)


addWorksheet(
  wb,
  "Model_counts"
)

writeData(
  wb,
  "Model_counts",
  model_summary
)


saveWorkbook(
  wb,
  GO_XLSX,
  overwrite = TRUE
)


# ============================================================
# 20. Final report
# ============================================================

message("")
message("==================================================")
message("SUPPLEMENTARY FIGURE S6 COMPLETE")
message("==================================================")

message("")
print(model_summary)

message("")
message(
  "Genes shared across >=2 models: ",
  nrow(shared_ge2)
)

message(
  "Shared genes used for GO: ",
  length(shared_genes)
)

message(
  "Common tested background: ",
  length(common_background)
)

message("")
message(
  "FDR-significant GO BP terms: ",
  sum(
    bp_df$p.adjust < 0.05,
    na.rm = TRUE
  )
)

message(
  "FDR-significant GO CC terms: ",
  sum(
    cc_df$p.adjust < 0.05,
    na.rm = TRUE
  )
)

message(
  "FDR-significant GO MF terms: ",
  sum(
    mf_df$p.adjust < 0.05,
    na.rm = TRUE
  )
)

message("")
message("Outputs:")
message(SHARED_TSV)
message(GO_XLSX)
message(GO_PDF)
message(GO_PNG)
