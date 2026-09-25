# ============================================================
# MATR3 KO
# QAPA differential APA + shortening/lengthening + GO
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(openxlsx)
  library(pheatmap)
  library(gprofiler2)
  library(stringr)
  library(scales)
})

BASE_DIR <- "results/05_MATR3_KO"
INPUT <- file.path(BASE_DIR, "MATR3_KO_QAPA_PAU.tsv")

XLSX_OUT <- file.path(BASE_DIR, "MATR3_KO_differential_APA.xlsx")
VOLCANO_OUT <- file.path(BASE_DIR, "MATR3_KO_volcano.pdf")
PCA_OUT <- file.path(BASE_DIR, "MATR3_KO_PCA.pdf")
HEATMAP_OUT <- file.path(BASE_DIR, "MATR3_KO_heatmap.pdf")
SL_OUT <- file.path(BASE_DIR, "MATR3_KO_shortening_lengthening.pdf")
MECH_VOLCANO_OUT <- file.path(BASE_DIR, "MATR3_KO_APA_site_mechanistic_volcano.pdf")

GO_DIR <- file.path(BASE_DIR, "GO_analysis")
dir.create(GO_DIR, recursive = TRUE, showWarnings = FALSE)
GO_XLSX <- file.path(GO_DIR, "MATR3_KO_GO_enrichment.xlsx")
GO_GENES_XLSX <- file.path(GO_DIR, "MATR3_KO_GO_gene_lists.xlsx")

TPM_CUTOFF <- 1
MIN_TPM_SAMPLES <- 2
DPAU_CUTOFF <- 10
FDR_CUTOFF <- 0.05
NOMINAL_P_CUTOFF <- 0.05
LENGTH_CHANGE_CUTOFF <- 100
TOP_LABELS <- 10
TOP_HEATMAP <- 50
TOP_GO <- 15

control_pau <- c("WT1_Rep1.PAU", "WT1_Rep2.PAU", "WT1_Rep3.PAU", "WT7_Rep1.PAU", "WT7_Rep2.PAU")
case_pau <- c("MATR3_KO1_Rep1.PAU", "MATR3_KO1_Rep2.PAU", "MATR3_KO5_Rep1.PAU")
control_tpm <- c("WT1_Rep1.TPM", "WT1_Rep2.TPM", "WT1_Rep3.TPM", "WT7_Rep1.TPM", "WT7_Rep2.TPM")
case_tpm <- c("MATR3_KO1_Rep1.TPM", "MATR3_KO1_Rep2.TPM", "MATR3_KO5_Rep1.TPM")
all_pau <- c(control_pau, case_pau)
all_tpm <- c(control_tpm, case_tpm)

df <- read.delim(INPUT, header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)

required_cols <- c("APA_ID","Gene","Gene_Name","Length","Num_Events",all_pau,all_tpm)
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) stop("Missing columns: ", paste(missing_cols, collapse=", "))

for (x in c(all_pau,all_tpm,"Length","Num_Events","UTR3.Start","UTR3.End","LastExon.Start","LastExon.End")) {
  if (x %in% colnames(df)) df[[x]] <- suppressWarnings(as.numeric(df[[x]]))
}

df$TPM_pass_n <- rowSums(as.data.frame(df[,all_tpm,drop=FALSE]) >= TPM_CUTOFF, na.rm=TRUE)
df$Expression_Pass <- df$TPM_pass_n >= MIN_TPM_SAMPLES

df$Control_mean_PAU <- rowMeans(df[,control_pau,drop=FALSE], na.rm=TRUE)
df$Case_mean_PAU <- rowMeans(df[,case_pau,drop=FALSE], na.rm=TRUE)
df$DeltaPAU <- df$Case_mean_PAU - df$Control_mean_PAU

safe_ttest <- function(x,y) {
  x <- as.numeric(x); y <- as.numeric(y)
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x)<2 || length(y)<2) return(NA_real_)
  if (length(unique(x))==1 && length(unique(y))==1 && unique(x)==unique(y)) return(1)
  tryCatch(t.test(x,y,paired=FALSE,var.equal=FALSE)$p.value, error=function(e) NA_real_)
}

message("Calculating event-level statistics...")
df$P_value <- vapply(seq_len(nrow(df)), function(i) {
  if (!df$Expression_Pass[i]) return(NA_real_)
  safe_ttest(df[i,case_pau], df[i,control_pau])
}, numeric(1))

df$FDR <- NA_real_
tested <- which(df$Expression_Pass & is.finite(df$P_value))
df$FDR[tested] <- p.adjust(df$P_value[tested], method="BH")

df$APA_Class <- case_when(
  grepl("_P$",df$APA_ID) ~ "Proximal",
  grepl("_D$",df$APA_ID) ~ "Distal",
  grepl("_S$",df$APA_ID) ~ "Single",
  TRUE ~ "Other"
)
df$Site <- case_when(
  df$APA_Class=="Proximal" ~ "P",
  df$APA_Class=="Distal" ~ "D",
  df$APA_Class=="Single" ~ "S",
  TRUE ~ "Other"
)

df$Significance <- "NS"
df$Significance[df$Expression_Pass & is.finite(df$FDR) & df$FDR<FDR_CUTOFF & df$DeltaPAU>=DPAU_CUTOFF] <- "Up"
df$Significance[df$Expression_Pass & is.finite(df$FDR) & df$FDR<FDR_CUTOFF & df$DeltaPAU<=-DPAU_CUTOFF] <- "Down"
df$Significance <- factor(df$Significance, levels=c("Down","NS","Up"))

summary_table <- data.frame(
  Metric=c("Total APA events","Unique genes","Events passing TPM filter","Events statistically tested",
           "FDR < 0.05","FDR < 0.05 and |DeltaPAU| >= 10","Significant increased PAU","Significant decreased PAU"),
  Value=c(nrow(df),length(unique(df$Gene)),sum(df$Expression_Pass),sum(is.finite(df$P_value)),
          sum(df$FDR<0.05,na.rm=TRUE),sum(df$FDR<0.05 & abs(df$DeltaPAU)>=10,na.rm=TRUE),
          sum(df$Significance=="Up"),sum(df$Significance=="Down"))
)

df_out <- df %>% arrange(FDR, desc(abs(DeltaPAU)))
sig_df <- df_out %>% filter(Significance!="NS")

# ---------------- Mechanistic shortening/lengthening volcano ----------------
#
# Blue = Shortening
# Red  = Lengthening
# Gray = NS
#
# Mechanistic definitions:
#   Proximal increased -> Shortening
#   Distal decreased   -> Shortening
#   Proximal decreased -> Lengthening
#   Distal increased   -> Lengthening
#
# Thresholds:
#   nominal P < 0.05
#   |DeltaPAU| >= 10
#
# Labels:
#   top 10 unique shortening genes
#   top 10 unique lengthening genes
# ============================================================

volcano <- df %>%
  filter(
    Expression_Pass,
    is.finite(P_value)
  )

volcano$Mechanistic_category <- "NS"

# Proximal increased -> shortening
volcano$Mechanistic_category[
  volcano$Site == "P" &
  volcano$P_value < NOMINAL_P_CUTOFF &
  volcano$DeltaPAU >= DPAU_CUTOFF
] <- "Proximal increased"

# Distal decreased -> shortening
volcano$Mechanistic_category[
  volcano$Site == "D" &
  volcano$P_value < NOMINAL_P_CUTOFF &
  volcano$DeltaPAU <= -DPAU_CUTOFF
] <- "Distal decreased"

# Proximal decreased -> lengthening
volcano$Mechanistic_category[
  volcano$Site == "P" &
  volcano$P_value < NOMINAL_P_CUTOFF &
  volcano$DeltaPAU <= -DPAU_CUTOFF
] <- "Proximal decreased"

# Distal increased -> lengthening
volcano$Mechanistic_category[
  volcano$Site == "D" &
  volcano$P_value < NOMINAL_P_CUTOFF &
  volcano$DeltaPAU >= DPAU_CUTOFF
] <- "Distal increased"

volcano$Direction <- case_when(
  volcano$Mechanistic_category %in% c(
    "Proximal increased",
    "Distal decreased"
  ) ~ "Shortening",

  volcano$Mechanistic_category %in% c(
    "Proximal decreased",
    "Distal increased"
  ) ~ "Lengthening",

  TRUE ~ "NS"
)

volcano$minus_log10_P <- -log10(
  pmax(
    volcano$P_value,
    1e-300
  )
)

# ============================================================
# Top 10 shortening genes
# ============================================================

label_short <- volcano %>%
  filter(
    Direction == "Shortening",
    !is.na(Gene_Name),
    Gene_Name != ""
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
    n = TOP_LABELS
  )

# ============================================================
# Top 10 lengthening genes
# ============================================================

label_long <- volcano %>%
  filter(
    Direction == "Lengthening",
    !is.na(Gene_Name),
    Gene_Name != ""
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
    n = TOP_LABELS
  )

message("")
message("Top 10 shortening genes:")
print(
  label_short %>%
    select(
      Gene_Name,
      APA_Class,
      DeltaPAU,
      P_value,
      FDR,
      Mechanistic_category
    )
)

message("")
message("Top 10 lengthening genes:")
print(
  label_long %>%
    select(
      Gene_Name,
      APA_Class,
      DeltaPAU,
      P_value,
      FDR,
      Mechanistic_category
    )
)

# ============================================================
# Plot
# ============================================================

p_volcano <- ggplot(
  volcano,
  aes(
    x = DeltaPAU,
    y = minus_log10_P
  )
) +

  # NS = gray
  geom_point(
    data = subset(
      volcano,
      Direction == "NS"
    ),
    color = "grey72",
    size = 1.5,
    alpha = 0.55
  ) +

  # Shortening = blue
  geom_point(
    data = subset(
      volcano,
      Direction == "Shortening"
    ),
    color = "#377EB8",
    size = 2.1,
    alpha = 0.90
  ) +

  # Lengthening = red
  geom_point(
    data = subset(
      volcano,
      Direction == "Lengthening"
    ),
    color = "#E41A1C",
    size = 2.1,
    alpha = 0.90
  ) +

  # DeltaPAU thresholds
  geom_vline(
    xintercept = c(
      -DPAU_CUTOFF,
      DPAU_CUTOFF
    ),
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey35"
  ) +

  # P = 0.05 threshold
  geom_hline(
    yintercept = -log10(
      NOMINAL_P_CUTOFF
    ),
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey35"
  ) +

  # Top 10 shortening labels
  geom_label_repel(
    data = label_short,
    aes(
      label = Gene_Name
    ),
    fill = "#D5E9F3",
    color = "black",
    size = 3,
    label.size = 0.18,
    box.padding = 0.55,
    point.padding = 0.30,
    force = 5,
    segment.color = "grey40",
    segment.size = 0.30,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 101
  ) +

  # Top 10 lengthening labels
  geom_label_repel(
    data = label_long,
    aes(
      label = Gene_Name
    ),
    fill = "#F3D0CD",
    color = "black",
    size = 3,
    label.size = 0.18,
    box.padding = 0.55,
    point.padding = 0.30,
    force = 5,
    segment.color = "grey40",
    segment.size = 0.30,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 202
  ) +

  labs(
    title = "MATR3 KO differential alternative polyadenylation",
    subtitle = paste0(
      "Nominal P < ",
      NOMINAL_P_CUTOFF,
      " and |DeltaPAU| >= ",
      DPAU_CUTOFF,
      "; blue = shortening, red = lengthening"
    ),
    x = "DeltaPAU (MATR3 KO - Control)",
    y = "-log10(P value)"
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
      hjust = 0.5
    ),
    legend.position = "none"
  )

ggsave(
  VOLCANO_OUT,
  p_volcano,
  width = 8.5,
  height = 6.5,
  device = cairo_pdf
)


# ---------------- PCA ----------------
pca_df <- df %>% filter(Expression_Pass,Num_Events>=2)
pca_mat <- as.matrix(pca_df[,all_pau,drop=FALSE])
rownames(pca_mat) <- make.unique(pca_df$APA_ID)
pca_mat <- pca_mat[apply(pca_mat,1,function(x) all(is.finite(x))),,drop=FALSE]
vv <- apply(pca_mat,1,var)
pca_mat <- pca_mat[is.finite(vv)&vv>0,,drop=FALSE]
if (nrow(pca_mat)>5000) {
  vv <- apply(pca_mat,1,var)
  pca_mat <- pca_mat[order(vv,decreasing=TRUE)[1:5000],,drop=FALSE]
}
pca <- prcomp(t(pca_mat),center=TRUE,scale.=TRUE)
variance <- (pca$sdev^2/sum(pca$sdev^2))*100
pca_plot_df <- data.frame(Sample=rownames(pca$x),PC1=pca$x[,1],PC2=pca$x[,2],stringsAsFactors=FALSE)
pca_plot_df$Group <- ifelse(pca_plot_df$Sample %in% control_pau,"Control","MATR3 KO")
p_pca <- ggplot(pca_plot_df,aes(PC1,PC2,color=Group,label=Sample)) +
  geom_point(size=4) +
  geom_text_repel(size=3.2,show.legend=FALSE) +
  stat_ellipse(aes(group=Group),type="norm",linewidth=0.8,linetype=2,show.legend=FALSE) +
  labs(title="QAPA PAU principal component analysis",
       x=paste0("PC1 (",round(variance[1],1),"%)"),
       y=paste0("PC2 (",round(variance[2],1),"%)")) +
  theme_classic(base_size=13) +
  theme(plot.title=element_text(face="bold"),legend.title=element_blank())
ggsave(PCA_OUT,p_pca,width=7,height=6,device=cairo_pdf)

# ---------------- Top 50 APA heatmap ----------------
heat_sig <- df %>% filter(Significance!="NS") %>% arrange(FDR,desc(abs(DeltaPAU))) %>% distinct(APA_ID,.keep_all=TRUE)
heat_other <- df %>% filter(Expression_Pass,is.finite(P_value),!APA_ID %in% heat_sig$APA_ID) %>%
  arrange(P_value,desc(abs(DeltaPAU))) %>% distinct(APA_ID,.keep_all=TRUE)
heat_df <- bind_rows(heat_sig,heat_other) %>% slice_head(n=TOP_HEATMAP)
heat_mat <- as.matrix(heat_df[,all_pau,drop=FALSE])
rownames(heat_mat) <- make.unique(paste0(heat_df$Gene_Name," | ",heat_df$APA_Class))
heat_z <- t(scale(t(heat_mat))); heat_z[!is.finite(heat_z)] <- 0
annotation_col <- data.frame(Group=c(rep("Control",length(control_pau)),rep("MATR3 KO",length(case_pau))))
rownames(annotation_col) <- all_pau
pdf(HEATMAP_OUT,width=8.5,height=10,useDingbats=FALSE)
pheatmap(heat_z,cluster_rows=TRUE,cluster_cols=TRUE,annotation_col=annotation_col,
         show_rownames=TRUE,show_colnames=TRUE,fontsize_row=7,fontsize_col=9,border_color=NA,
         main=paste0("Top ",nrow(heat_z)," differential APA events"))
dev.off()

# ---------------- Gene-level weighted 3'UTR length ----------------
multi_df <- df %>% filter(Expression_Pass,Num_Events>=2,is.finite(Length),Length>0)

calc_weighted_length <- function(dat,pau_col) {
  pau <- dat[[pau_col]]
  good <- is.finite(pau) & is.finite(dat$Length)
  if (sum(good)==0) return(NA_real_)
  pau <- pau[good]; lengths <- dat$Length[good]
  total <- sum(pau)
  if (!is.finite(total) || total<=0) return(NA_real_)
  sum((pau/total)*lengths)
}

gene_list <- split(multi_df,multi_df$Gene)
length_list <- lapply(gene_list,function(g) {
  vals <- sapply(all_pau,function(s) calc_weighted_length(g,s))
  z <- data.frame(Gene=unique(g$Gene)[1],Gene_Name=unique(g$Gene_Name)[1],Num_APA_Sites=nrow(g),stringsAsFactors=FALSE)
  for (s in all_pau) z[[s]] <- vals[[s]]
  z
})
length_df <- bind_rows(length_list)
length_df$Control_mean_length <- rowMeans(length_df[,control_pau,drop=FALSE],na.rm=TRUE)
length_df$Case_mean_length <- rowMeans(length_df[,case_pau,drop=FALSE],na.rm=TRUE)
length_df$Delta_UTR_length <- length_df$Case_mean_length-length_df$Control_mean_length
length_df$P_value <- vapply(seq_len(nrow(length_df)),function(i) safe_ttest(length_df[i,case_pau],length_df[i,control_pau]),numeric(1))
length_df$FDR <- p.adjust(length_df$P_value,method="BH")
length_df$Direction <- case_when(
  is.finite(length_df$Delta_UTR_length) & length_df$Delta_UTR_length<=-LENGTH_CHANGE_CUTOFF ~ "Shortening",
  is.finite(length_df$Delta_UTR_length) & length_df$Delta_UTR_length>= LENGTH_CHANGE_CUTOFF ~ "Lengthening",
  TRUE ~ "Stable"
)
length_df$minus_log10_P <- -log10(pmax(length_df$P_value,1e-300))

top_short_genes <- length_df %>% filter(Direction=="Shortening",is.finite(P_value),is.finite(Delta_UTR_length)) %>%
  arrange(P_value,Delta_UTR_length) %>% distinct(Gene_Name,.keep_all=TRUE) %>% slice_head(n=TOP_LABELS)
top_long_genes <- length_df %>% filter(Direction=="Lengthening",is.finite(P_value),is.finite(Delta_UTR_length)) %>%
  arrange(P_value,desc(Delta_UTR_length)) %>% distinct(Gene_Name,.keep_all=TRUE) %>% slice_head(n=TOP_LABELS)

# ============================================================
# Publication-ready gene-level 3'UTR shortening/lengthening plot
#
# Shortening = blue
# Lengthening = red
# NS = gray
#
# Labels:
#   top 10 shortening genes
#   top 10 lengthening genes
#
# Improvements:
#   - cleaner typography
#   - lighter NS points
#   - balanced label placement
#   - explicit legend
#   - more breathing room around labels
#   - publication-friendly 8.6 x 6.4 inch layout
# ============================================================

length_df$PlotDirection <- factor(
  length_df$Direction,
  levels = c(
    "Shortening",
    "Lengthening",
    "Stable"
  ),
  labels = c(
    "Shortening",
    "Lengthening",
    "NS"
  )
)

# Symmetric x-range around zero, based on the observed data.
x_abs_max <- max(
  abs(length_df$Delta_UTR_length / 1000),
  na.rm = TRUE
)

if (!is.finite(x_abs_max) || x_abs_max <= 0) {
  x_abs_max <- 1
}

x_limit <- x_abs_max * 1.08

# Slightly different label nudges for the two directions.
top_short_genes$LabelX <- (
  top_short_genes$Delta_UTR_length / 1000
) - 0.10 * x_limit

top_long_genes$LabelX <- (
  top_long_genes$Delta_UTR_length / 1000
) + 0.10 * x_limit

p_sl <- ggplot(
  length_df,
  aes(
    x = Delta_UTR_length / 1000,
    y = minus_log10_P
  )
) +

  # ----------------------------------------------------------
  # NS points
  # ----------------------------------------------------------
  geom_point(
    data = subset(
      length_df,
      PlotDirection == "NS"
    ),
    aes(color = PlotDirection),
    size = 1.45,
    alpha = 0.32
  ) +

  # ----------------------------------------------------------
  # Significant shortening / lengthening
  # ----------------------------------------------------------
  geom_point(
    data = subset(
      length_df,
      PlotDirection != "NS"
    ),
    aes(color = PlotDirection),
    size = 2.15,
    alpha = 0.88
  ) +

  # ----------------------------------------------------------
  # Thresholds
  # ----------------------------------------------------------
  geom_vline(
    xintercept = c(
      -LENGTH_CHANGE_CUTOFF / 1000,
       LENGTH_CHANGE_CUTOFF / 1000
    ),
    linetype = "dashed",
    linewidth = 0.38,
    color = "grey40"
  ) +

  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.38,
    color = "grey40"
  ) +

  # Zero reference line
  geom_vline(
    xintercept = 0,
    linewidth = 0.35,
    color = "grey65"
  ) +

  # ----------------------------------------------------------
  # Top 10 shortening labels
  # ----------------------------------------------------------
  geom_label_repel(
    data = top_short_genes,
    aes(
      label = Gene_Name,
      x = Delta_UTR_length / 1000,
      y = minus_log10_P
    ),
    fill = "#DCECF5",
    color = "black",
    size = 2.8,
    label.size = 0.18,
    label.padding = unit(0.14, "lines"),
    box.padding = 0.55,
    point.padding = 0.28,
    force = 8,
    force_pull = 0.35,
    direction = "both",
    min.segment.length = 0,
    segment.color = "grey45",
    segment.size = 0.28,
    max.overlaps = Inf,
    seed = 101
  ) +

  # ----------------------------------------------------------
  # Top 10 lengthening labels
  # ----------------------------------------------------------
  geom_label_repel(
    data = top_long_genes,
    aes(
      label = Gene_Name,
      x = Delta_UTR_length / 1000,
      y = minus_log10_P
    ),
    fill = "#F7D8D6",
    color = "black",
    size = 2.8,
    label.size = 0.18,
    label.padding = unit(0.14, "lines"),
    box.padding = 0.55,
    point.padding = 0.28,
    force = 8,
    force_pull = 0.35,
    direction = "both",
    min.segment.length = 0,
    segment.color = "grey45",
    segment.size = 0.28,
    max.overlaps = Inf,
    seed = 202
  ) +

  # ----------------------------------------------------------
  # Manual colors
  # ----------------------------------------------------------
  scale_color_manual(
    values = c(
      "Shortening" = "#2C7FB8",
      "Lengthening" = "#D73027",
      "NS" = "grey72"
    ),
    breaks = c(
      "Shortening",
      "Lengthening",
      "NS"
    ),
    name = NULL
  ) +

  # ----------------------------------------------------------
  # Balanced x-axis
  # ----------------------------------------------------------
  scale_x_continuous(
    limits = c(
      -x_limit,
       x_limit
    ),
    breaks = scales::pretty_breaks(
      n = 6
    ),
    expand = expansion(
      mult = c(
        0.015,
        0.015
      )
    )
  ) +

  # ----------------------------------------------------------
  # Titles
  # ----------------------------------------------------------
  labs(
    title = "MATR3 KO alters 3'UTR length",
    subtitle = "Gene-level PAU-weighted 3'UTR usage",
    x = expression(
      Delta*"3'UTR length ("*"MATR3 KO"*" - Control, kb)"
    ),
    y = expression(
      -log[10]*"(P value)"
    )
  ) +

  # ----------------------------------------------------------
  # Publication theme
  # ----------------------------------------------------------
  theme_classic(
    base_size = 11
  ) +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 14,
      hjust = 0.5,
      margin = margin(
        b = 4
      )
    ),

    plot.subtitle = element_text(
      size = 10,
      hjust = 0.5,
      color = "grey25",
      margin = margin(
        b = 9
      )
    ),

    axis.title.x = element_text(
      size = 11,
      margin = margin(
        t = 7
      )
    ),

    axis.title.y = element_text(
      size = 11,
      margin = margin(
        r = 7
      )
    ),

    axis.text = element_text(
      size = 9.5,
      color = "black"
    ),

    axis.line = element_line(
      linewidth = 0.45,
      color = "black"
    ),

    axis.ticks = element_line(
      linewidth = 0.4,
      color = "black"
    ),

    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",

    legend.text = element_text(
      size = 9.5
    ),

    legend.key.width = unit(
      0.8,
      "lines"
    ),

    legend.spacing.x = unit(
      0.25,
      "cm"
    ),

    plot.margin = margin(
      t = 8,
      r = 18,
      b = 8,
      l = 8
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

ggsave(
  SL_OUT,
  p_sl,
  width = 8.6,
  height = 6.4,
  device = cairo_pdf
)

# Also save a 600-dpi PNG for quick viewing / figure assembly.
ggsave(
  sub(
    "\\.pdf$",
    ".png",
    SL_OUT
  ),
  p_sl,
  width = 8.6,
  height = 6.4,
  dpi = 600
)


# ---------------- Mechanistic P/D classification ----------------
mech <- df %>% filter(Expression_Pass,is.finite(P_value))
mech$Mechanistic_category <- "NS"
mech$Mechanistic_category[mech$Site=="P" & mech$P_value<NOMINAL_P_CUTOFF & mech$DeltaPAU>= DPAU_CUTOFF] <- "Proximal increased"
mech$Mechanistic_category[mech$Site=="D" & mech$P_value<NOMINAL_P_CUTOFF & mech$DeltaPAU<=-DPAU_CUTOFF] <- "Distal decreased"
mech$Mechanistic_category[mech$Site=="P" & mech$P_value<NOMINAL_P_CUTOFF & mech$DeltaPAU<=-DPAU_CUTOFF] <- "Proximal decreased"
mech$Mechanistic_category[mech$Site=="D" & mech$P_value<NOMINAL_P_CUTOFF & mech$DeltaPAU>= DPAU_CUTOFF] <- "Distal increased"
mech$Direction <- case_when(
  mech$Mechanistic_category %in% c("Proximal increased","Distal decreased") ~ "Shortening",
  mech$Mechanistic_category %in% c("Proximal decreased","Distal increased") ~ "Lengthening",
  TRUE ~ "NS"
)
mech$minus_log10_P <- -log10(pmax(mech$P_value,1e-300))

p_mech <- ggplot(mech,aes(DeltaPAU,minus_log10_P)) +
  geom_point(data=subset(mech,Direction=="NS"),color="grey78",size=1.4,alpha=0.5) +
  geom_point(data=subset(mech,Direction=="Shortening"),color="#377EB8",size=2,alpha=0.85) +
  geom_point(data=subset(mech,Direction=="Lengthening"),color="#E41A1C",size=2,alpha=0.85) +
  geom_vline(xintercept=c(-DPAU_CUTOFF,DPAU_CUTOFF),linetype="dashed",linewidth=0.45) +
  geom_hline(yintercept=-log10(NOMINAL_P_CUTOFF),linetype="dashed",linewidth=0.45) +
  labs(title="MATR3 KO: mechanistic APA-site changes",
       subtitle="Blue = 3'UTR shortening; red = 3'UTR lengthening",
       x="DeltaPAU (MATR3 KO - Control)",y="-log10(P value)") +
  theme_classic(base_size=13)
ggsave(MECH_VOLCANO_OUT,p_mech,width=8,height=6.5,device=cairo_pdf)

# ---------------- Differential Excel ----------------
wb <- createWorkbook()
addWorksheet(wb,"Summary"); writeData(wb,"Summary",summary_table)
addWorksheet(wb,"All_APA_events"); writeData(wb,"All_APA_events",df_out)
addWorksheet(wb,"Significant_APA"); writeData(wb,"Significant_APA",sig_df)
addWorksheet(wb,"UTR_length"); writeData(wb,"UTR_length",length_df %>% arrange(desc(abs(Delta_UTR_length))))
addWorksheet(wb,"Shortening"); writeData(wb,"Shortening",length_df %>% filter(Direction=="Shortening") %>% arrange(Delta_UTR_length))
addWorksheet(wb,"Lengthening"); writeData(wb,"Lengthening",length_df %>% filter(Direction=="Lengthening") %>% arrange(desc(Delta_UTR_length)))
addWorksheet(wb,"Mechanistic_APA"); writeData(wb,"Mechanistic_APA",mech %>% arrange(P_value,desc(abs(DeltaPAU))))
saveWorkbook(wb,XLSX_OUT,overwrite=TRUE)

# ---------------- GO ----------------
short_genes <- mech %>% filter(Direction=="Shortening") %>% pull(Gene_Name) %>% unique() %>% na.omit()
long_genes <- mech %>% filter(Direction=="Lengthening") %>% pull(Gene_Name) %>% unique() %>% na.omit()
discordant <- intersect(short_genes,long_genes)
short_genes_clean <- setdiff(short_genes,discordant)
long_genes_clean <- setdiff(long_genes,discordant)
background <- df %>% filter(Site %in% c("P","D")) %>% pull(Gene_Name) %>% unique() %>% na.omit()

wb_genes <- createWorkbook()
addWorksheet(wb_genes,"Shortening"); writeData(wb_genes,"Shortening",data.frame(Gene=short_genes_clean))
addWorksheet(wb_genes,"Lengthening"); writeData(wb_genes,"Lengthening",data.frame(Gene=long_genes_clean))
addWorksheet(wb_genes,"Discordant"); writeData(wb_genes,"Discordant",data.frame(Gene=discordant))
addWorksheet(wb_genes,"Background"); writeData(wb_genes,"Background",data.frame(Gene=background))
saveWorkbook(wb_genes,GO_GENES_XLSX,overwrite=TRUE)

run_GO <- function(genes) {
  if (length(genes)<2) return(NULL)
  gost(query=genes,organism="hsapiens",ordered_query=FALSE,custom_bg=background,
       correction_method="fdr",sources=c("GO:BP","GO:CC","GO:MF"),significant=FALSE)
}
GO_short <- run_GO(short_genes_clean)
GO_long <- run_GO(long_genes_clean)
short_result <- if (!is.null(GO_short) && !is.null(GO_short$result)) GO_short$result else data.frame()
long_result <- if (!is.null(GO_long) && !is.null(GO_long$result)) GO_long$result else data.frame()

wb_go <- createWorkbook()
for (s in c("Shortening_BP","Shortening_CC","Shortening_MF","Lengthening_BP","Lengthening_CC","Lengthening_MF")) addWorksheet(wb_go,s)
write_go <- function(wb,sheet,dat,source_name) {
  if (nrow(dat)==0) return(invisible(NULL))
  x <- dat %>% filter(source==source_name) %>% arrange(p_value)
  if (nrow(x)>0) writeData(wb,sheet,x)
}
write_go(wb_go,"Shortening_BP",short_result,"GO:BP")
write_go(wb_go,"Shortening_CC",short_result,"GO:CC")
write_go(wb_go,"Shortening_MF",short_result,"GO:MF")
write_go(wb_go,"Lengthening_BP",long_result,"GO:BP")
write_go(wb_go,"Lengthening_CC",long_result,"GO:CC")
write_go(wb_go,"Lengthening_MF",long_result,"GO:MF")
saveWorkbook(wb_go,GO_XLSX,overwrite=TRUE)

plot_GO <- function(dat,direction,ontology,outfile,top_n=TOP_GO) {
  if (nrow(dat)==0) return(NULL)
  source_name <- paste0("GO:",ontology)
  x <- dat %>% filter(source==source_name,is.finite(p_value),p_value>0) %>% arrange(p_value) %>% slice_head(n=top_n)
  if (nrow(x)==0) return(NULL)
  x <- x %>% mutate(GeneRatio=intersection_size/query_size,GeneCount=intersection_size,
                    GO_term=stringr::str_wrap(term_name,width=42)) %>% arrange(desc(p_value))
  x$GO_term <- factor(x$GO_term,levels=x$GO_term)
  ontology_text <- case_when(ontology=="BP"~"biological processes",ontology=="CC"~"cellular components",ontology=="MF"~"molecular functions",TRUE~ontology)
  subtitle_text <- if ("significant" %in% colnames(x) && any(x$significant==TRUE,na.rm=TRUE)) "FDR-significant GO terms" else "Top exploratory terms; none significant after FDR correction"
  xmax <- max(x$GeneRatio,na.rm=TRUE)*1.12
  p <- ggplot(x,aes(GeneRatio,GO_term)) +
    geom_segment(aes(x=0,xend=GeneRatio,y=GO_term,yend=GO_term),linewidth=0.35,color="grey50") +
    geom_point(aes(size=GeneCount,color=p_value),alpha=0.95) +
    scale_color_viridis_c(option="viridis",direction=-1,name="FDR adjusted\np-value") +
    scale_size_continuous(name="Gene count",range=c(2.4,7.2),breaks=scales::pretty_breaks(n=4)) +
    scale_x_continuous(limits=c(0,xmax),breaks=scales::pretty_breaks(n=5),labels=scales::label_number(accuracy=0.01),expand=expansion(mult=c(0,0))) +
    labs(title=paste0("GO ",ontology_text," associated with MATR3 KO 3'UTR ",tolower(direction)," genes"),
         subtitle=subtitle_text,x="GeneRatio",y=NULL) +
    theme_bw(base_size=8) +
    theme(panel.background=element_rect(fill="#EAF2F7",color=NA),panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(),plot.title=element_text(size=8,hjust=0.5),
          plot.subtitle=element_text(size=8,hjust=0.5),axis.text=element_text(size=8,color="black"),
          axis.title=element_text(size=8),legend.title=element_text(size=8),legend.text=element_text(size=8),
          plot.margin=margin(7,8,7,28))
  ggsave(outfile,p,width=8.8,height=5.8,device=cairo_pdf)
}

plot_GO(short_result,"Shortening","BP",file.path(GO_DIR,"MATR3_KO_shortening_GO_BP.pdf"))
plot_GO(short_result,"Shortening","CC",file.path(GO_DIR,"MATR3_KO_shortening_GO_CC.pdf"))
plot_GO(short_result,"Shortening","MF",file.path(GO_DIR,"MATR3_KO_shortening_GO_MF.pdf"))
plot_GO(long_result,"Lengthening","BP",file.path(GO_DIR,"MATR3_KO_lengthening_GO_BP.pdf"))
plot_GO(long_result,"Lengthening","CC",file.path(GO_DIR,"MATR3_KO_lengthening_GO_CC.pdf"))
plot_GO(long_result,"Lengthening","MF",file.path(GO_DIR,"MATR3_KO_lengthening_GO_MF.pdf"))

message("==============================================")
message("ANALYSIS COMPLETE: MATR3 KO")
message("==============================================")
print(summary_table)
message("Gene-level 3'UTR direction:")
print(table(length_df$Direction))
message("GO workbook: ",GO_XLSX)
