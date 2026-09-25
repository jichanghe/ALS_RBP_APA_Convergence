# ============================================================
# 5-MODEL APA INTEGRATION
#
# Models:
#   TDP-43 KO
#   TDP-43 K263E
#   FUS KO
#   FUS P525L
#   MATR3 KO
#
# Input:
#   *_differential_APA.xlsx
#   sheet: Mechanistic_APA (preferred) or All_APA_events
#
# Main mechanistic definition:
#   Proximal increased  -> Shortening
#   Distal decreased    -> Shortening
#   Proximal decreased  -> Lengthening
#   Distal increased    -> Lengthening
#
# Cross-model significance:
#   nominal P < 0.05 and |DeltaPAU| >= 10
#
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(UpSetR)
  library(scales)
  library(gprofiler2)
})

OUT_DIR <- "results/APA_5model_integration"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

MODELS <- c(
  "TDP_KO",
  "TDP_K263E",
  "FUS_KO",
  "FUS_P525L",
  "MATR3_KO"
)

MODEL_LABELS <- c(
  TDP_KO = "TDP-43 KO",
  TDP_K263E = "TDP-43 K263E",
  FUS_KO = "FUS KO",
  FUS_P525L = "FUS P525L",
  MATR3_KO = "MATR3 KO"
)

APA_FILES <- c(
  TDP_KO =
    "data/TDP43_KD_differential_APA.xlsx",
  TDP_K263E =
    "data/TDP43_K263E_differential_APA.xlsx",
  FUS_KO =
    "data/FUS_KO_differential_APA.xlsx",
  FUS_P525L =
    "data/FUS_P525L_differential_APA.xlsx",
  MATR3_KO =
    "data/MATR3_KO_differential_APA.xlsx"
)

P_CUTOFF <- 0.05
FDR_CUTOFF <- 0.05
DPAU_CUTOFF <- 10
HEATMAP_MIN_MODELS <- 2
RUN_GO <- TRUE
GO_SOURCES <- c("GO:BP","GO:CC","GO:MF")
TOP_GO_TERMS <- 40
COMMON_GO_MIN_MODELS <- 2

find_column <- function(df, candidates) {
  original_names <- colnames(df)
  clean_names <- tolower(gsub("[^a-z0-9]","",original_names))
  candidate_clean <- tolower(gsub("[^a-z0-9]","",candidates))
  idx <- match(candidate_clean,clean_names)
  idx <- idx[!is.na(idx)]
  if (length(idx)==0) return(NA_character_)
  original_names[idx[1]]
}

read_apa_file <- function(model) {

  infile <- APA_FILES[[model]]

  if (!file.exists(infile)) {
    stop("Missing APA workbook: ", infile)
  }

  message("Reading: ", infile)

  sheets <- excel_sheets(infile)

  sheet_to_read <- if ("Mechanistic_APA" %in% sheets) {
    "Mechanistic_APA"
  } else if ("All_APA_events" %in% sheets) {
    "All_APA_events"
  } else {
    stop("No Mechanistic_APA or All_APA_events sheet in: ", infile)
  }

  df <- as.data.frame(
    read_excel(infile,sheet=sheet_to_read),
    stringsAsFactors=FALSE
  )

  apa_id_col <- find_column(df,c("APA_ID","APAID","event_id"))
  gene_col <- find_column(df,c("Gene_Name","gene_name","gene","symbol"))
  dpau_col <- find_column(df,c("DeltaPAU","delta_PAU","dPAU"))
  p_col <- find_column(df,c("P_value","Pvalue","p_value","pvalue","P"))
  fdr_col <- find_column(df,c("FDR","padj","adj_p","adjusted_p_value"))
  site_col <- find_column(df,c("Site","APA_Class","APAClass"))

  if (is.na(apa_id_col)) stop("APA_ID missing: ", infile)
  if (is.na(gene_col)) stop("Gene_Name missing: ", infile)
  if (is.na(dpau_col)) stop("DeltaPAU missing: ", infile)
  if (is.na(p_col)) stop("P_value missing: ", infile)

  out <- data.frame(
    Model=model,
    ModelLabel=unname(MODEL_LABELS[model]),
    APA_ID=as.character(df[[apa_id_col]]),
    Gene=as.character(df[[gene_col]]),
    DeltaPAU=suppressWarnings(as.numeric(df[[dpau_col]])),
    P_value=suppressWarnings(as.numeric(df[[p_col]])),
    stringsAsFactors=FALSE
  )

  out$FDR <- if (!is.na(fdr_col)) {
    suppressWarnings(as.numeric(df[[fdr_col]]))
  } else {
    NA_real_
  }

  site_raw <- if (!is.na(site_col)) {
    as.character(df[[site_col]])
  } else {
    rep(NA_character_,nrow(df))
  }

  out$Site <- case_when(
    grepl("^P$|Proximal",site_raw,ignore.case=TRUE) ~ "P",
    grepl("^D$|Distal",site_raw,ignore.case=TRUE) ~ "D",
    grepl("^S$|Single",site_raw,ignore.case=TRUE) ~ "S",
    grepl("_P$",out$APA_ID) ~ "P",
    grepl("_D$",out$APA_ID) ~ "D",
    grepl("_S$",out$APA_ID) ~ "S",
    TRUE ~ "Other"
  )

  out$Mechanistic_category <- "NS"

  out$Mechanistic_category[
    out$Site=="P" &
    is.finite(out$P_value) &
    out$P_value<P_CUTOFF &
    out$DeltaPAU>=DPAU_CUTOFF
  ] <- "Proximal increased"

  out$Mechanistic_category[
    out$Site=="D" &
    is.finite(out$P_value) &
    out$P_value<P_CUTOFF &
    out$DeltaPAU<=-DPAU_CUTOFF
  ] <- "Distal decreased"

  out$Mechanistic_category[
    out$Site=="P" &
    is.finite(out$P_value) &
    out$P_value<P_CUTOFF &
    out$DeltaPAU<=-DPAU_CUTOFF
  ] <- "Proximal decreased"

  out$Mechanistic_category[
    out$Site=="D" &
    is.finite(out$P_value) &
    out$P_value<P_CUTOFF &
    out$DeltaPAU>=DPAU_CUTOFF
  ] <- "Distal increased"

  out$Direction <- case_when(
    out$Mechanistic_category %in%
      c("Proximal increased","Distal decreased") ~ "Shortening",
    out$Mechanistic_category %in%
      c("Proximal decreased","Distal increased") ~ "Lengthening",
    TRUE ~ "NS"
  )

  out$MechanisticSignificant <- out$Direction!="NS"

  out$FDRSignificant <- (
    is.finite(out$FDR) &
    out$FDR<FDR_CUTOFF &
    abs(out$DeltaPAU)>=DPAU_CUTOFF
  )

  out
}

all_data <- bind_rows(
  lapply(MODELS,read_apa_file)
)

message("Total APA rows: ", format(nrow(all_data),big.mark=","))

# ============================================================
# COUNTS
# ============================================================

direction_counts <- all_data %>%
  count(Model,Direction,name="N_events")

mechanistic_counts <- all_data %>%
  filter(MechanisticSignificant) %>%
  count(Model,Mechanistic_category,Direction,name="N_events")

write.csv(direction_counts,
          file.path(OUT_DIR,"APA_direction_counts.csv"),
          row.names=FALSE)

write.csv(mechanistic_counts,
          file.path(OUT_DIR,"APA_mechanistic_counts.csv"),
          row.names=FALSE)

# ============================================================
# EVENT x MODEL DeltaPAU MATRIX
# ============================================================

event_model_best <- all_data %>%
  filter(!is.na(APA_ID),APA_ID!="") %>%
  arrange(P_value) %>%
  group_by(Model,APA_ID) %>%
  slice(1) %>%
  ungroup()

dpau_matrix_df <- event_model_best %>%
  select(APA_ID,Model,DeltaPAU) %>%
  pivot_wider(names_from=Model,values_from=DeltaPAU)

for (m in MODELS) {
  if (!m %in% colnames(dpau_matrix_df)) dpau_matrix_df[[m]] <- NA_real_
}

dpau_matrix_df <- dpau_matrix_df %>%
  select(APA_ID,all_of(MODELS))

write.csv(
  dpau_matrix_df,
  file.path(OUT_DIR,"All_APA_events_5model_DeltaPAU_matrix.csv"),
  row.names=FALSE
)

# ============================================================
# EVENT MEMBERSHIP + UPSET
# ============================================================

sig_event_membership <- all_data %>%
  filter(MechanisticSignificant) %>%
  distinct(APA_ID,Model) %>%
  mutate(Present=1L) %>%
  pivot_wider(names_from=Model,values_from=Present,values_fill=0L)

for (m in MODELS) {
  if (!m %in% colnames(sig_event_membership)) sig_event_membership[[m]] <- 0L
}

sig_event_membership <- sig_event_membership %>%
  select(APA_ID,all_of(MODELS)) %>%
  mutate(N_models=rowSums(across(all_of(MODELS))))

write.csv(
  sig_event_membership,
  file.path(OUT_DIR,"Significant_APA_event_membership.csv"),
  row.names=FALSE
)

pdf(file.path(OUT_DIR,"Figure_APA_Event_Level_UpSet.pdf"),
    width=11,height=7)

UpSetR::upset(
  as.data.frame(sig_event_membership[,MODELS]),
  sets=MODELS,
  keep.order=TRUE,
  order.by="freq",
  nsets=length(MODELS),
  nintersects=40,
  mainbar.y.label="Shared differential APA events",
  sets.x.label="Differential APA events per model",
  text.scale=1.3
)

dev.off()

# ============================================================
# GENE MEMBERSHIP + UPSET
# ============================================================

sig_gene_membership <- all_data %>%
  filter(MechanisticSignificant,!is.na(Gene),Gene!="") %>%
  distinct(Gene,Model) %>%
  mutate(Present=1L) %>%
  pivot_wider(names_from=Model,values_from=Present,values_fill=0L)

for (m in MODELS) {
  if (!m %in% colnames(sig_gene_membership)) sig_gene_membership[[m]] <- 0L
}

sig_gene_membership <- sig_gene_membership %>%
  select(Gene,all_of(MODELS)) %>%
  mutate(N_models=rowSums(across(all_of(MODELS)))) %>%
  arrange(desc(N_models),Gene)

write.csv(
  sig_gene_membership,
  file.path(OUT_DIR,"Significant_APA_gene_membership.csv"),
  row.names=FALSE
)

pdf(file.path(OUT_DIR,"Figure_APA_Gene_Level_UpSet.pdf"),
    width=11,height=7)

UpSetR::upset(
  as.data.frame(sig_gene_membership[,MODELS]),
  sets=MODELS,
  keep.order=TRUE,
  order.by="freq",
  nsets=length(MODELS),
  nintersects=40,
  mainbar.y.label="Shared APA-regulated genes",
  sets.x.label="APA-regulated genes per model",
  text.scale=1.3
)

dev.off()

# ============================================================
# SHORTENING / LENGTHENING GENE MEMBERSHIP
# ============================================================

build_direction_membership <- function(direction) {

  x <- all_data %>%
    filter(Direction==direction,!is.na(Gene),Gene!="") %>%
    distinct(Gene,Model) %>%
    mutate(Present=1L) %>%
    pivot_wider(names_from=Model,values_from=Present,values_fill=0L)

  for (m in MODELS) {
    if (!m %in% colnames(x)) x[[m]] <- 0L
  }

  x %>%
    select(Gene,all_of(MODELS)) %>%
    mutate(N_models=rowSums(across(all_of(MODELS)))) %>%
    arrange(desc(N_models),Gene)
}

short_membership <- build_direction_membership("Shortening")
long_membership  <- build_direction_membership("Lengthening")

write.csv(short_membership,
          file.path(OUT_DIR,"Shortening_gene_membership.csv"),
          row.names=FALSE)

write.csv(long_membership,
          file.path(OUT_DIR,"Lengthening_gene_membership.csv"),
          row.names=FALSE)

plot_direction_upset <- function(membership,direction,outfile) {

  pdf(outfile,width=11,height=7)

  UpSetR::upset(
    as.data.frame(membership[,MODELS]),
    sets=MODELS,
    keep.order=TRUE,
    order.by="freq",
    nsets=length(MODELS),
    nintersects=40,
    mainbar.y.label=paste0(
      "Shared 3'UTR ",
      tolower(direction),
      " genes"
    ),
    sets.x.label=paste0(direction," genes per model"),
    text.scale=1.3
  )

  dev.off()
}

plot_direction_upset(
  short_membership,
  "Shortening",
  file.path(OUT_DIR,"Figure_Shortening_Gene_UpSet.pdf")
)

plot_direction_upset(
  long_membership,
  "Lengthening",
  file.path(OUT_DIR,"Figure_Lengthening_Gene_UpSet.pdf")
)

# ============================================================
# 5-MODEL DeltaPAU HEATMAP
# ============================================================

heatmap_events <- sig_event_membership %>%
  filter(N_models>=HEATMAP_MIN_MODELS) %>%
  pull(APA_ID)

heatmap_df <- dpau_matrix_df %>%
  filter(APA_ID %in% heatmap_events)

heatmap_mat <- as.matrix(heatmap_df[,MODELS,drop=FALSE])
rownames(heatmap_mat) <- heatmap_df$APA_ID

keep_rows <- apply(
  heatmap_mat,
  1,
  function(x) any(is.finite(x))
)

heatmap_mat <- heatmap_mat[keep_rows,,drop=FALSE]

heatmap_plot_mat <- heatmap_mat
heatmap_plot_mat[is.na(heatmap_plot_mat)] <- 0
heatmap_plot_mat[heatmap_plot_mat>50] <- 50
heatmap_plot_mat[heatmap_plot_mat< -50] <- -50

if (nrow(heatmap_plot_mat)>0) {

  pdf(
    file.path(OUT_DIR,"Figure_5Model_DeltaPAU_Heatmap.pdf"),
    width=7,
    height=max(7,min(16,nrow(heatmap_plot_mat)*0.05))
  )

  pheatmap(
    heatmap_plot_mat,
    cluster_rows=TRUE,
    cluster_cols=FALSE,
    color=colorRampPalette(c("#2166AC","white","#B2182B"))(100),
    breaks=seq(-50,50,length.out=101),
    border_color=NA,
    show_rownames=FALSE,
    labels_col=unname(MODEL_LABELS[MODELS]),
    main=paste0(
      "Differential APA events significant in ≥",
      HEATMAP_MIN_MODELS,
      " models"
    )
  )

  dev.off()
}

# ============================================================
# CORRELATION MATRIX
# ============================================================

cor_matrix <- cor(
  dpau_matrix_df[,MODELS,drop=FALSE],
  use="pairwise.complete.obs",
  method="pearson"
)

write.csv(
  cor_matrix,
  file.path(OUT_DIR,"DeltaPAU_Pearson_correlation_matrix.csv")
)

pdf(
  file.path(OUT_DIR,"Figure_DeltaPAU_Correlation_Matrix.pdf"),
  width=7,
  height=6.5
)

pheatmap(
  cor_matrix,
  cluster_rows=TRUE,
  cluster_cols=TRUE,
  display_numbers=TRUE,
  number_format="%.2f",
  color=colorRampPalette(c("#2166AC","white","#B2182B"))(100),
  breaks=seq(-1,1,length.out=101),
  border_color="white",
  main="Cross-model DeltaPAU correlation"
)

dev.off()

# ============================================================
# KO vs MUTANT CONCORDANCE
#
# MAIN-FIGURE DESIGN:
#   Only APA events significant in BOTH conditions are plotted.
#
# Classification:
#   Concordant shortening = both models show shortening
#   Concordant lengthening = both models show lengthening
#   Discordant = opposite directional APA effects
#
# Colors:
#   blue = concordant shortening
#   red  = concordant lengthening
#   gray = discordant
#
# Global all-event correlation is still exported separately,
# but it is NOT emphasized in the publication scatterplot.
# ============================================================


# ============================================================
# Helper: retrieve mechanistic direction for each APA event/model
# ============================================================

direction_lookup <- all_data %>%

  filter(
    MechanisticSignificant,
    !is.na(APA_ID),
    APA_ID != ""
  ) %>%

  arrange(
    P_value
  ) %>%

  group_by(
    Model,
    APA_ID
  ) %>%

  slice(1) %>%

  ungroup() %>%

  select(
    Model,
    APA_ID,
    Gene,
    Site,
    DeltaPAU,
    P_value,
    FDR,
    Mechanistic_category,
    Direction
  )


# ============================================================
# Publication-ready KO-vs-mutant function
# ============================================================

make_shared_significant_concordance_plot <- function(
  model_x,
  model_y,
  outfile_prefix
) {

  # ----------------------------------------------------------
  # Significant APA events in both models
  # ----------------------------------------------------------

  x_df <- direction_lookup %>%

    filter(
      Model == model_x
    ) %>%

    transmute(
      APA_ID,
      Gene_x = Gene,
      Site_x = Site,
      DeltaPAU_x = DeltaPAU,
      P_x = P_value,
      FDR_x = FDR,
      Mechanistic_x = Mechanistic_category,
      Direction_x = Direction
    )


  y_df <- direction_lookup %>%

    filter(
      Model == model_y
    ) %>%

    transmute(
      APA_ID,
      Gene_y = Gene,
      Site_y = Site,
      DeltaPAU_y = DeltaPAU,
      P_y = P_value,
      FDR_y = FDR,
      Mechanistic_y = Mechanistic_category,
      Direction_y = Direction
    )


  shared_df <- inner_join(
    x_df,
    y_df,
    by = "APA_ID"
  )


  # ----------------------------------------------------------
  # Unified gene / site annotation
  # ----------------------------------------------------------

  shared_df$Gene <- ifelse(
    !is.na(shared_df$Gene_x) &
    shared_df$Gene_x != "",
    shared_df$Gene_x,
    shared_df$Gene_y
  )


  shared_df$Site <- ifelse(
    !is.na(shared_df$Site_x) &
    shared_df$Site_x != "",
    shared_df$Site_x,
    shared_df$Site_y
  )


  # ----------------------------------------------------------
  # Concordance class
  # ----------------------------------------------------------

  shared_df$Concordance <- case_when(

    shared_df$Direction_x == "Shortening" &
    shared_df$Direction_y == "Shortening" ~
      "Concordant shortening",

    shared_df$Direction_x == "Lengthening" &
    shared_df$Direction_y == "Lengthening" ~
      "Concordant lengthening",

    TRUE ~
      "Discordant"
  )


  shared_df$Concordance <- factor(
    shared_df$Concordance,
    levels = c(
      "Concordant shortening",
      "Concordant lengthening",
      "Discordant"
    )
  )


  # ----------------------------------------------------------
  # Directional concordance percentage
  # ----------------------------------------------------------

  n_shared <- nrow(
    shared_df
  )


  n_concordant <- sum(
    shared_df$Concordance %in% c(
      "Concordant shortening",
      "Concordant lengthening"
    ),
    na.rm = TRUE
  )


  concordance_percent <- if (
    n_shared > 0
  ) {

    100 *
      n_concordant /
      n_shared

  } else {

    NA_real_
  }


  # ----------------------------------------------------------
  # Shared-event Pearson correlation
  # ----------------------------------------------------------

  shared_r <- if (
    n_shared >= 3
  ) {

    cor(
      shared_df$DeltaPAU_x,
      shared_df$DeltaPAU_y,
      use = "complete.obs",
      method = "pearson"
    )

  } else {

    NA_real_
  }


  # ----------------------------------------------------------
  # Overall correlation across the union of significant events
  # retained only for supplemental statistics
  # ----------------------------------------------------------

  sig_union <- sig_event_membership %>%

    filter(
      .data[[model_x]] == 1 |
      .data[[model_y]] == 1
    ) %>%

    pull(
      APA_ID
    )


  union_df <- dpau_matrix_df %>%

    filter(
      APA_ID %in% sig_union
    ) %>%

    select(
      APA_ID,
      all_of(
        c(
          model_x,
          model_y
        )
      )
    ) %>%

    filter(
      !is.na(.data[[model_x]]),
      !is.na(.data[[model_y]])
    )


  union_r <- if (
    nrow(union_df) >= 3
  ) {

    cor(
      union_df[[model_x]],
      union_df[[model_y]],
      use = "complete.obs",
      method = "pearson"
    )

  } else {

    NA_real_
  }


  # ----------------------------------------------------------
  # Symmetric plot limits
  # ----------------------------------------------------------

  max_abs <- max(
    abs(
      c(
        shared_df$DeltaPAU_x,
        shared_df$DeltaPAU_y
      )
    ),
    na.rm = TRUE
  )


  if (
    !is.finite(max_abs) ||
    max_abs <= 0
  ) {
    max_abs <- 10
  }


  plot_limit <- ceiling(
    max_abs / 10
  ) * 10


  # ----------------------------------------------------------
  # Publication scatter
  # ----------------------------------------------------------

  p <- ggplot(
    shared_df,
    aes(
      x = DeltaPAU_x,
      y = DeltaPAU_y
    )
  ) +

    # Light quadrant shading:
    # lower-left = shared shortening
    annotate(
      "rect",
      xmin = -Inf,
      xmax = 0,
      ymin = -Inf,
      ymax = 0,
      fill = "#EAF4FA",
      alpha = 0.55
    ) +

    # upper-right = shared lengthening
    annotate(
      "rect",
      xmin = 0,
      xmax = Inf,
      ymin = 0,
      ymax = Inf,
      fill = "#FBEDEC",
      alpha = 0.55
    ) +

    # Reference axes
    geom_hline(
      yintercept = 0,
      linewidth = 0.40,
      color = "grey45"
    ) +

    geom_vline(
      xintercept = 0,
      linewidth = 0.40,
      color = "grey45"
    ) +

    # y = x reference
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      linewidth = 0.45,
      color = "grey55"
    ) +

    geom_point(
      aes(
        color = Concordance
      ),
      size = 2.35,
      alpha = 0.90
    ) +

    scale_color_manual(
      values = c(
        "Concordant shortening" = "#2C7FB8",
        "Concordant lengthening" = "#D73027",
        "Discordant" = "grey55"
      ),
      name = NULL
    ) +

    scale_x_continuous(
      limits = c(
        -plot_limit,
        plot_limit
      ),
      breaks = scales::pretty_breaks(
        n = 5
      ),
      expand = expansion(
        mult = c(
          0.02,
          0.02
        )
      )
    ) +

    scale_y_continuous(
      limits = c(
        -plot_limit,
        plot_limit
      ),
      breaks = scales::pretty_breaks(
        n = 5
      ),
      expand = expansion(
        mult = c(
          0.02,
          0.02
        )
      )
    ) +

    coord_equal() +

    labs(
      title = paste0(
        MODEL_LABELS[
          model_x
        ],
        " vs ",
        MODEL_LABELS[
          model_y
        ]
      ),

      subtitle = paste0(
        format(
          concordance_percent,
          digits = 3
        ),
        "% directionally concordant among ",
        n_shared,
        " shared significant APA events"
      ),

      x = paste0(
        MODEL_LABELS[
          model_x
        ],
        " DeltaPAU"
      ),

      y = paste0(
        MODEL_LABELS[
          model_y
        ],
        " DeltaPAU"
      )
    ) +

    theme_classic(
      base_size = 10
    ) +

    theme(
      plot.title = element_text(
        face = "bold",
        size = 12,
        hjust = 0.5,
        margin = margin(
          b = 3
        )
      ),

      plot.subtitle = element_text(
        size = 9,
        hjust = 0.5,
        color = "grey25",
        margin = margin(
          b = 7
        )
      ),

      axis.title = element_text(
        size = 10
      ),

      axis.text = element_text(
        size = 9,
        color = "black"
      ),

      legend.position = "top",
      legend.direction = "horizontal",
      legend.text = element_text(
        size = 8.5
      ),

      axis.line = element_line(
        linewidth = 0.45
      ),

      axis.ticks = element_line(
        linewidth = 0.40
      ),

      plot.margin = margin(
        8,
        10,
        8,
        8
      )
    ) +

    guides(
      color = guide_legend(
        override.aes = list(
          size = 3,
          alpha = 1
        )
      )
    )


  # ----------------------------------------------------------
  # Save main-figure PDF / PNG
  # ----------------------------------------------------------

  ggsave(
    file.path(
      OUT_DIR,
      paste0(
        outfile_prefix,
        "_shared_significant.pdf"
      )
    ),
    p,
    width = 6.2,
    height = 5.6,
    device = cairo_pdf
  )


  ggsave(
    file.path(
      OUT_DIR,
      paste0(
        outfile_prefix,
        "_shared_significant.png"
      )
    ),
    p,
    width = 6.2,
    height = 5.6,
    dpi = 600
  )


  # ----------------------------------------------------------
  # Save the exact events
  # ----------------------------------------------------------

  write.csv(
    shared_df,
    file.path(
      OUT_DIR,
      paste0(
        outfile_prefix,
        "_shared_significant_events.csv"
      )
    ),
    row.names = FALSE
  )


  write.csv(
    union_df,
    file.path(
      OUT_DIR,
      paste0(
        outfile_prefix,
        "_union_events_for_supplement.csv"
      )
    ),
    row.names = FALSE
  )


  # ----------------------------------------------------------
  # Counts by concordance class
  # ----------------------------------------------------------

  concordance_counts <- shared_df %>%

    count(
      Concordance,
      name = "N_events"
    ) %>%

    mutate(
      Percent =
        100 *
        N_events /
        sum(
          N_events
        )
    )


  write.csv(
    concordance_counts,
    file.path(
      OUT_DIR,
      paste0(
        outfile_prefix,
        "_concordance_counts.csv"
      )
    ),
    row.names = FALSE
  )


  # ----------------------------------------------------------
  # Return summary stats
  # ----------------------------------------------------------

  data.frame(

    Model_1 =
      model_x,

    Model_2 =
      model_y,

    N_shared_significant =
      n_shared,

    N_concordant =
      n_concordant,

    Direction_concordance_percent =
      concordance_percent,

    Shared_event_Pearson_r =
      shared_r,

    Union_event_Pearson_r =
      union_r,

    stringsAsFactors = FALSE
  )
}


# ============================================================
# TDP-43 KO vs K263E
# ============================================================

tdp_stats <- make_shared_significant_concordance_plot(

  model_x = "TDP_KO",

  model_y = "TDP_K263E",

  outfile_prefix =
    "Figure_TDP_KO_vs_K263E_DeltaPAU"
)


# ============================================================
# FUS KO vs P525L
# ============================================================

fus_stats <- make_shared_significant_concordance_plot(

  model_x = "FUS_KO",

  model_y = "FUS_P525L",

  outfile_prefix =
    "Figure_FUS_KO_vs_P525L_DeltaPAU"
)


concordance_stats <- bind_rows(
  tdp_stats,
  fus_stats
)


write.csv(
  concordance_stats,
  file.path(
    OUT_DIR,
    "KO_vs_mutant_APA_concordance_statistics.csv"
  ),
  row.names = FALSE
)


# ============================================================
# Main-figure directional concordance summary
# ============================================================

concordance_summary_plot_df <- concordance_stats %>%

  mutate(
    Comparison = case_when(

      Model_1 == "TDP_KO" ~
        "TDP-43 KO vs K263E",

      Model_1 == "FUS_KO" ~
        "FUS KO vs P525L",

      TRUE ~
        paste(
          Model_1,
          "vs",
          Model_2
        )
    )
  )


concordance_summary_plot_df$Comparison <- factor(
  concordance_summary_plot_df$Comparison,
  levels = c(
    "TDP-43 KO vs K263E",
    "FUS KO vs P525L"
  )
)


p_concordance_summary <- ggplot(
  concordance_summary_plot_df,
  aes(
    x = Comparison,
    y = Direction_concordance_percent
  )
) +

  geom_col(
    width = 0.62,
    fill = "grey35"
  ) +

  geom_text(
    aes(
      label = paste0(
        round(
          Direction_concordance_percent,
          1
        ),
        "%"
      )
    ),
    vjust = -0.45,
    size = 3.2
  ) +

  scale_y_continuous(
    limits = c(
      0,
      105
    ),
    breaks = seq(
      0,
      100,
      by = 20
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  labs(
    title = "Directional concordance of shared APA events",
    x = NULL,
    y = "Concordant shared events (%)"
  ) +

  theme_classic(
    base_size = 10
  ) +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 11,
      hjust = 0.5
    ),

    axis.text.x = element_text(
      size = 9
    ),

    axis.title.y = element_text(
      size = 10
    )
  )


ggsave(
  file.path(
    OUT_DIR,
    "Figure_KO_Mutant_Directional_Concordance_Summary.pdf"
  ),
  p_concordance_summary,
  width = 5.2,
  height = 4.5,
  device = cairo_pdf
)


ggsave(
  file.path(
    OUT_DIR,
    "Figure_KO_Mutant_Directional_Concordance_Summary.png"
  ),
  p_concordance_summary,
  width = 5.2,
  height = 4.5,
  dpi = 600
)


# ============================================================
# SHARED EVENTS / GENES
# ============================================================

event_annotation <- all_data %>%
  select(APA_ID,Gene,Site) %>%
  filter(!is.na(APA_ID)) %>%
  distinct(APA_ID,.keep_all=TRUE)

core_events <- sig_event_membership %>%
  left_join(event_annotation,by="APA_ID") %>%
  arrange(desc(N_models),Gene)

write.csv(
  core_events,
  file.path(OUT_DIR,"Core_shared_APA_events_all.csv"),
  row.names=FALSE
)

write.csv(
  core_events %>% filter(N_models>=3),
  file.path(OUT_DIR,"Core_shared_APA_events_3plus_models.csv"),
  row.names=FALSE
)

write.csv(
  sig_gene_membership,
  file.path(OUT_DIR,"Core_shared_APA_genes.csv"),
  row.names=FALSE
)

write.csv(
  short_membership %>% filter(N_models>=2),
  file.path(OUT_DIR,"Core_shared_shortening_genes_2plus_models.csv"),
  row.names=FALSE
)

write.csv(
  long_membership %>% filter(N_models>=2),
  file.path(OUT_DIR,"Core_shared_lengthening_genes_2plus_models.csv"),
  row.names=FALSE
)

# ============================================================
# SHORTENING / LENGTHENING COUNTS
# ============================================================

direction_plot_df <- all_data %>%
  filter(Direction %in% c("Shortening","Lengthening")) %>%
  count(Model,Direction,name="N")

direction_plot_df$Model <- factor(
  direction_plot_df$Model,
  levels=MODELS
)

p_counts <- ggplot(
  direction_plot_df,
  aes(
    x=Model,
    y=N,
    fill=Direction
  )
) +
  geom_col(
    position=position_dodge(width=0.8),
    width=0.7
  ) +
  scale_fill_manual(
    values=c(
      "Shortening"="#2C7FB8",
      "Lengthening"="#D73027"
    )
  ) +
  scale_x_discrete(labels=MODEL_LABELS) +
  labs(
    title="Alternative polyadenylation remodeling across ALS-RBP models",
    x=NULL,
    y="Number of differential APA events",
    fill=NULL
  ) +
  theme_classic(base_size=12) +
  theme(
    axis.text.x=element_text(angle=30,hjust=1),
    plot.title=element_text(face="bold",hjust=0.5),
    legend.position="top"
  )

ggsave(
  file.path(OUT_DIR,"Figure_APA_Shortening_Lengthening_Counts.pdf"),
  p_counts,
  width=8,
  height=6
)

# ============================================================
# GO GENE LISTS
# ============================================================

short_gene_lists <- list()
long_gene_lists <- list()

for (model in MODELS) {

  short_genes <- all_data %>%
    filter(Model==model,Direction=="Shortening",!is.na(Gene),Gene!="") %>%
    pull(Gene) %>%
    unique()

  long_genes <- all_data %>%
    filter(Model==model,Direction=="Lengthening",!is.na(Gene),Gene!="") %>%
    pull(Gene) %>%
    unique()

  discordant <- intersect(short_genes,long_genes)

  short_gene_lists[[model]] <- setdiff(short_genes,discordant)
  long_gene_lists[[model]] <- setdiff(long_genes,discordant)
}

go_background <- all_data %>%
  filter(Site %in% c("P","D"),!is.na(Gene),Gene!="") %>%
  pull(Gene) %>%
  unique()

# ============================================================
# GO ENRICHMENT
# ============================================================

go_all <- NULL

if (RUN_GO) {

  go_results <- list()

  for (model in MODELS) {

    for (direction in c("Shortening","Lengthening")) {

      genes <- if (direction=="Shortening") {
        short_gene_lists[[model]]
      } else {
        long_gene_lists[[model]]
      }

      message(
        "GO: ",
        model,
        " / ",
        direction,
        " (",
        length(genes),
        " genes)"
      )

      if (length(genes)<5) next

      gost_res <- tryCatch(
        gost(
          query=genes,
          organism="hsapiens",
          ordered_query=FALSE,
          multi_query=FALSE,
          significant=FALSE,
          correction_method="fdr",
          sources=GO_SOURCES,
          custom_bg=go_background
        ),
        error=function(e) NULL
      )

      if (is.null(gost_res) || is.null(gost_res$result)) next

      tmp <- gost_res$result
      tmp$Model <- model
      tmp$Direction <- direction
      tmp$minusLog10FDR <- -log10(pmax(tmp$p_value,1e-300))

      go_results[[paste(model,direction,sep="_")]] <- tmp
    }
  }

  if (length(go_results)>0) {

    go_all <- bind_rows(go_results)

    go_all_export <- go_all %>%
      mutate(
        across(
          where(is.list),
          ~vapply(
            .x,
            function(z) {
              if (length(z)==0 || all(is.na(z))) return(NA_character_)
              paste(as.character(z),collapse=";")
            },
            character(1)
          )
        )
      )

    write.csv(
      go_all_export,
      file.path(
        OUT_DIR,
        "GO_enrichment_all_models_shortening_lengthening.csv"
      ),
      row.names=FALSE
    )

    common_go <- go_all %>%
      filter(p_value<0.05) %>%
      distinct(
        Model,Direction,source,term_id,term_name,
        .keep_all=TRUE
      ) %>%
      group_by(Direction,source,term_id,term_name) %>%
      summarise(
        N_models=n_distinct(Model),
        MeanScore=mean(minusLog10FDR,na.rm=TRUE),
        .groups="drop"
      ) %>%
      filter(N_models>=COMMON_GO_MIN_MODELS) %>%
      arrange(Direction,desc(N_models),desc(MeanScore))

    write.csv(
      common_go,
      file.path(
        OUT_DIR,
        "Common_GO_terms_shortening_lengthening.csv"
      ),
      row.names=FALSE
    )

    for (direction in c("Shortening","Lengthening")) {

      top_go <- common_go %>%
        filter(Direction==direction) %>%
        slice_head(n=TOP_GO_TERMS)

      if (nrow(top_go)==0) next

      selected_terms <- top_go$term_id

      go_heat <- go_all %>%
        filter(
          Direction==direction,
          term_id %in% selected_terms
        ) %>%
        select(
          term_id,term_name,source,Model,minusLog10FDR
        ) %>%
        group_by(
          term_id,term_name,source,Model
        ) %>%
        summarise(
          score=max(minusLog10FDR,na.rm=TRUE),
          .groups="drop"
        ) %>%
        pivot_wider(
          names_from=Model,
          values_from=score,
          values_fill=0
        )

      for (m in MODELS) {
        if (!m %in% colnames(go_heat)) go_heat[[m]] <- 0
      }

      go_matrix <- as.matrix(
        go_heat[,MODELS,drop=FALSE]
      )

      rownames(go_matrix) <- make.unique(
        paste0(go_heat$source," | ",go_heat$term_name)
      )

      pdf(
        file.path(
          OUT_DIR,
          paste0(
            "Figure_Common_GO_",
            direction,
            "_Heatmap.pdf"
          )
        ),
        width=9,
        height=max(8,nrow(go_matrix)*0.25)
      )

      pheatmap(
        go_matrix,
        cluster_rows=TRUE,
        cluster_cols=FALSE,
        labels_col=unname(MODEL_LABELS[MODELS]),
        color=colorRampPalette(
          c("white","yellow","orange","red")
        )(100),
        border_color=NA,
        main=paste0(
          "Common GO enrichment: 3'UTR ",
          tolower(direction)
        ),
        fontsize_row=8
      )

      dev.off()
    }
  }
}


# ============================================================
# PANEL F: CONVERGENT APA GO PATHWAYS
#
# Purpose:
#   Generate ONE publication-ready functional summary panel
#   from the five-model integrated APA analysis.
#
# File:
#   Figure_Convergent_APA_GO_Pathways.pdf
#
# Definition of convergent pathway:
#   GO term with nominally enriched result (g:Profiler adjusted
#   p-value < 0.05) in >= COMMON_GO_MIN_MODELS models.
#
# Plot:
#   x-axis  = mean -log10(FDR-adjusted p-value) across models
#   y-axis  = GO term
#   size    = number of models showing enrichment
#   color   = APA direction
#
#   blue = Shortening
#   red  = Lengthening
#
# If no GO term is shared by >=2 models, the script still
# generates the PDF using the strongest recurrent/exploratory
# terms and clearly states this in the subtitle.
# ============================================================

PANEL_F_PDF <- file.path(
  OUT_DIR,
  "Figure_Convergent_APA_GO_Pathways.pdf"
)

PANEL_F_PNG <- file.path(
  OUT_DIR,
  "Figure_Convergent_APA_GO_Pathways.png"
)

panel_f_data <- data.frame()

if (
  RUN_GO &&
  !is.null(go_all) &&
  nrow(go_all) > 0
) {

  # ----------------------------------------------------------
  # Build a model-level GO table
  # Keep the best result for each model/direction/GO term.
  # ----------------------------------------------------------

  go_model_level <- go_all %>%

    filter(
      is.finite(p_value),
      p_value > 0
    ) %>%

    distinct(
      Model,
      Direction,
      source,
      term_id,
      term_name,
      .keep_all = TRUE
    ) %>%

    mutate(
      minusLog10FDR =
        -log10(
          pmax(
            p_value,
            1e-300
          )
        ),

      Enriched =
        p_value < 0.05
    )


  # ----------------------------------------------------------
  # Summarize recurrence across models
  # ----------------------------------------------------------

  panel_f_all <- go_model_level %>%

    group_by(
      Direction,
      source,
      term_id,
      term_name
    ) %>%

    summarise(
      N_models =
        sum(
          Enriched,
          na.rm = TRUE
        ),

      Mean_minusLog10FDR =
        mean(
          minusLog10FDR[
            Enriched
          ],
          na.rm = TRUE
        ),

      Max_minusLog10FDR =
        max(
          minusLog10FDR[
            Enriched
          ],
          na.rm = TRUE
        ),

      Models =
        paste(
          MODEL_LABELS[
            unique(
              Model[
                Enriched
              ]
            )
          ],
          collapse = "; "
        ),

      .groups = "drop"
    )


  # Replace NaN / -Inf created when no model passes p<0.05
  panel_f_all$Mean_minusLog10FDR[
    !is.finite(
      panel_f_all$Mean_minusLog10FDR
    )
  ] <- 0

  panel_f_all$Max_minusLog10FDR[
    !is.finite(
      panel_f_all$Max_minusLog10FDR
    )
  ] <- 0


  # ----------------------------------------------------------
  # Primary selection:
  # GO terms enriched in >=2 models.
  # ----------------------------------------------------------

  panel_f_convergent <- panel_f_all %>%

    filter(
      N_models >= COMMON_GO_MIN_MODELS
    )


  panel_f_has_convergent <- (
    nrow(
      panel_f_convergent
    ) > 0
  )


  # ----------------------------------------------------------
  # Select top terms separately for shortening / lengthening.
  #
  # Primary ranking:
  #   1. number of models
  #   2. mean enrichment strength
  #
  # Up to 10 shortening + 10 lengthening terms.
  # ----------------------------------------------------------

  if (
    panel_f_has_convergent
  ) {

    panel_f_data <- panel_f_convergent %>%

      group_by(
        Direction
      ) %>%

      arrange(
        desc(N_models),
        desc(Mean_minusLog10FDR),
        .by_group = TRUE
      ) %>%

      slice_head(
        n = 10
      ) %>%

      ungroup()

    panel_f_subtitle <- paste0(
      "GO terms enriched in \u2265",
      COMMON_GO_MIN_MODELS,
      " ALS-RBP models"
    )

  } else {

    # --------------------------------------------------------
    # Graceful fallback:
    # If no term reaches recurrence in >=2 models, choose the
    # strongest terms observed in at least one model.
    # --------------------------------------------------------

    panel_f_data <- panel_f_all %>%

      filter(
        N_models >= 1
      ) %>%

      group_by(
        Direction
      ) %>%

      arrange(
        desc(N_models),
        desc(Mean_minusLog10FDR),
        .by_group = TRUE
      ) %>%

      slice_head(
        n = 10
      ) %>%

      ungroup()

    panel_f_subtitle <- paste0(
      "Exploratory GO summary; no term was enriched in \u2265",
      COMMON_GO_MIN_MODELS,
      " models"
    )
  }


  # ----------------------------------------------------------
  # Add readable GO labels
  # ----------------------------------------------------------

  if (
    nrow(
      panel_f_data
    ) > 0
  ) {

    panel_f_data <- panel_f_data %>%

      mutate(
        GO_label =
          paste0(
            source,
            ": ",
            term_name
          ),

        GO_label =
          stringr::str_wrap(
            GO_label,
            width = 48
          ),

        Direction =
          factor(
            Direction,
            levels = c(
              "Shortening",
              "Lengthening"
            )
          )
      )


    # --------------------------------------------------------
    # Order y-axis:
    # shortening and lengthening kept visually grouped, with
    # strongest terms at the top of each group.
    # --------------------------------------------------------

    panel_f_data <- panel_f_data %>%

      arrange(
        Direction,
        N_models,
        Mean_minusLog10FDR
      )


    panel_f_data$GO_label <- factor(
      panel_f_data$GO_label,
      levels = unique(
        panel_f_data$GO_label
      )
    )


    # --------------------------------------------------------
    # Publication-ready dot/lollipop plot
    # --------------------------------------------------------

    p_panel_f <- ggplot(
      panel_f_data,
      aes(
        x = Mean_minusLog10FDR,
        y = GO_label
      )
    ) +

      geom_segment(
        aes(
          x = 0,
          xend = Mean_minusLog10FDR,
          y = GO_label,
          yend = GO_label,
          color = Direction
        ),
        linewidth = 0.45,
        alpha = 0.55
      ) +

      geom_point(
        aes(
          size = N_models,
          color = Direction
        ),
        alpha = 0.95
      ) +

      scale_color_manual(
        values = c(
          "Shortening" = "#2C7FB8",
          "Lengthening" = "#D73027"
        ),
        name = "3'UTR change"
      ) +

      scale_size_continuous(
        name = "Models",
        range = c(
          2.8,
          7.0
        ),
        breaks = seq(
          1,
          length(MODELS),
          by = 1
        )
      ) +

      scale_x_continuous(
        expand = expansion(
          mult = c(
            0,
            0.08
          )
        ),
        breaks = scales::pretty_breaks(
          n = 5
        )
      ) +

      labs(
        title = "Functional pathways associated with convergent APA remodeling",
        subtitle = panel_f_subtitle,
        x = expression(
          "Mean " * -log[10] * "(FDR-adjusted p-value)"
        ),
        y = NULL
      ) +

      theme_classic(
        base_size = 8
      ) +

      theme(
        plot.title = element_text(
          size = 8,
          face = "bold",
          hjust = 0.5,
          margin = margin(
            b = 3
          )
        ),

        plot.subtitle = element_text(
          size = 8,
          hjust = 0.5,
          color = "grey30",
          margin = margin(
            b = 7
          )
        ),

        axis.title.x = element_text(
          size = 8,
          margin = margin(
            t = 6
          )
        ),

        axis.text.x = element_text(
          size = 8,
          color = "black"
        ),

        axis.text.y = element_text(
          size = 8,
          color = "black",
          lineheight = 0.95
        ),

        axis.line.y = element_blank(),

        axis.ticks.y = element_blank(),

        legend.position = "right",

        legend.title = element_text(
          size = 8
        ),

        legend.text = element_text(
          size = 8
        ),

        plot.margin = margin(
          t = 8,
          r = 10,
          b = 8,
          l = 24
        )
      )


    # --------------------------------------------------------
    # Save Panel F
    # --------------------------------------------------------

    ggsave(
      PANEL_F_PDF,
      p_panel_f,
      width = 8.8,
      height = 6.8,
      device = cairo_pdf
    )


    ggsave(
      PANEL_F_PNG,
      p_panel_f,
      width = 8.8,
      height = 6.8,
      dpi = 600
    )


    # --------------------------------------------------------
    # Save underlying Panel F data
    # --------------------------------------------------------

    write.csv(
      panel_f_data,
      file.path(
        OUT_DIR,
        "Figure_Convergent_APA_GO_Pathways_data.csv"
      ),
      row.names = FALSE
    )

  } else {

    # --------------------------------------------------------
    # Absolute fallback: create a PDF message rather than
    # silently failing to generate Panel F.
    # --------------------------------------------------------

    pdf(
      PANEL_F_PDF,
      width = 8.8,
      height = 6.8,
      useDingbats = FALSE
    )

    plot.new()

    text(
      0.5,
      0.58,
      "Functional pathways associated with convergent APA remodeling",
      cex = 1.2,
      font = 2
    )

    text(
      0.5,
      0.48,
      "No GO terms were available for the current shortening/lengthening gene sets.",
      cex = 1
    )

    dev.off()
  }

} else {

  # ----------------------------------------------------------
  # GO analysis unavailable / empty:
  # still generate Panel F PDF explaining why.
  # ----------------------------------------------------------

  pdf(
    PANEL_F_PDF,
    width = 8.8,
    height = 6.8,
    useDingbats = FALSE
  )

  plot.new()

  text(
    0.5,
    0.58,
    "Functional pathways associated with convergent APA remodeling",
    cex = 1.2,
    font = 2
  )

  text(
    0.5,
    0.48,
    "GO analysis returned no usable results.",
    cex = 1
  )

  dev.off()
}


# ============================================================
# EXCEL WORKBOOK
# ============================================================

wb <- createWorkbook()

addWorksheet(wb,"Summary")
summary_df <- data.frame(
  Parameter=c(
    "Nominal P cutoff",
    "FDR cutoff",
    "Absolute DeltaPAU cutoff",
    "Heatmap minimum models",
    "Models",
    "Direction definition"
  ),
  Value=c(
    P_CUTOFF,
    FDR_CUTOFF,
    DPAU_CUTOFF,
    HEATMAP_MIN_MODELS,
    paste(MODELS,collapse="; "),
    "P increased or D decreased = shortening; P decreased or D increased = lengthening"
  )
)
writeData(wb,"Summary",summary_df)

addWorksheet(wb,"APA_Counts")
writeData(wb,"APA_Counts",mechanistic_counts)

addWorksheet(wb,"DeltaPAU_Matrix")
writeData(wb,"DeltaPAU_Matrix",dpau_matrix_df)

addWorksheet(wb,"Event_Membership")
writeData(wb,"Event_Membership",core_events)

addWorksheet(wb,"Gene_Membership")
writeData(wb,"Gene_Membership",sig_gene_membership)

addWorksheet(wb,"Shortening_Membership")
writeData(wb,"Shortening_Membership",short_membership)

addWorksheet(wb,"Lengthening_Membership")
writeData(wb,"Lengthening_Membership",long_membership)

addWorksheet(wb,"KO_Mutant_Concordance")
writeData(wb,"KO_Mutant_Concordance",concordance_stats)

if (RUN_GO && !is.null(go_all)) {

  addWorksheet(wb,"GO_All")
  writeData(wb,"GO_All",go_all_export)

  if (exists("common_go")) {
    addWorksheet(wb,"GO_Common")
    writeData(wb,"GO_Common",common_go)
  }

  if (
    exists("panel_f_data") &&
    nrow(panel_f_data) > 0
  ) {
    addWorksheet(wb,"PanelF_GO")
    writeData(wb,"PanelF_GO",panel_f_data)
  }
}

saveWorkbook(
  wb,
  file.path(
    OUT_DIR,
    "APA_5model_integrated_results.xlsx"
  ),
  overwrite=TRUE
)

cat("\n============================================================\n")
cat("5-MODEL APA ANALYSIS FINISHED\n")
cat("============================================================\n")
cat("\nOutput directory:\n",OUT_DIR,"\n")
cat("\nMain outputs:\n")
cat("  Figure_APA_Event_Level_UpSet.pdf\n")
cat("  Figure_APA_Gene_Level_UpSet.pdf\n")
cat("  Figure_Shortening_Gene_UpSet.pdf\n")
cat("  Figure_Lengthening_Gene_UpSet.pdf\n")
cat("  Figure_5Model_DeltaPAU_Heatmap.pdf\n")
cat("  Figure_DeltaPAU_Correlation_Matrix.pdf\n")
cat("  Figure_TDP_KO_vs_K263E_DeltaPAU_shared_significant.pdf\n")
cat("  Figure_FUS_KO_vs_P525L_DeltaPAU_shared_significant.pdf\n")
cat("  Figure_APA_Shortening_Lengthening_Counts.pdf\n")
cat("  Figure_KO_Mutant_Directional_Concordance_Summary.pdf\n")
cat("  Figure_Convergent_APA_GO_Pathways.pdf\n")
cat("  APA_5model_integrated_results.xlsx\n")
cat("============================================================\n")
