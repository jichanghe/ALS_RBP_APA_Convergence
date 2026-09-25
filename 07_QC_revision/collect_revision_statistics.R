
# ============================================================
# Neuroscience revision:
# Quantitative summary of all QAPA datasets
#
# Legacy discovery/QC utility.
# This script recursively scans a complete QAPA analysis tree.
#
# Usage:
# Rscript collect_revision_statistics.R /path/to/QAPA_root
#
# Output:
# <QAPA_root>/Revision_quantitative_summary/
# ============================================================

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop(
    "\nUsage:\n",
    "Rscript collect_revision_statistics.R /path/to/QAPA_root\n"
  )
}

ROOT <- normalizePath(
  args[1],
  mustWork = TRUE
)

OUT <- file.path(
  ROOT,
  "Revision_quantitative_summary"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Packages
# ============================================================

need <- c("data.table")

for (pkg in need) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "\nMissing package: ", pkg,
      "\nInstall with:\n",
      "install.packages('", pkg, "')\n"
    )
  }
}

has_readxl <- requireNamespace("readxl", quietly = TRUE)

# ============================================================
# Dataset definitions
#
# More-specific mutant names MUST be tested before KO names.
# ============================================================

dataset_defs <- list(

  TDP43_K263E = c(
    "k263e",
    "tdp.?43.*k263e",
    "tardbp.*k263e"
  ),

  FUS_P525L = c(
    "p525l",
    "fus.*p525l"
  ),

  MATR3_KO = c(
    "matr3.*ko",
    "matr3.*kd",
    "matr3[_-]?ko",
    "matr3[_-]?kd"
  ),

  TDP43_KO = c(
    "tdp.?43.*ko",
    "tdp.?43.*kd",
    "tardbp.*ko",
    "tardbp.*kd",
    "tdp43[_-]?ko",
    "tdp43[_-]?kd"
  ),

  FUS_KO = c(
    "fus.*ko",
    "fus.*kd",
    "fus[_-]?ko",
    "fus[_-]?kd"
  )
)

# ============================================================
# Utilities
# ============================================================

clean_name <- function(x) {
  gsub("[^a-z0-9]+", "", tolower(x))
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

match_any <- function(x, patterns) {

  x <- tolower(x)

  any(
    vapply(
      patterns,
      function(p) grepl(p, x, ignore.case = TRUE, perl = TRUE),
      logical(1)
    )
  )
}

dataset_from_path <- function(path) {

  p <- tolower(path)

  for (nm in names(dataset_defs)) {

    if (match_any(p, dataset_defs[[nm]])) {
      return(nm)
    }
  }

  return(NA_character_)
}

# ============================================================
# Read table safely
# ============================================================

read_table_safe <- function(f, nrows = Inf) {

  ext <- tolower(tools::file_ext(f))

  tryCatch({

    if (ext %in% c("xlsx", "xls")) {

      if (!has_readxl) {
        return(NULL)
      }

      d <- readxl::read_excel(
        f,
        sheet = 1,
        n_max = if (is.finite(nrows)) nrows else Inf
      )

      return(as.data.frame(d))
    }

    d <- data.table::fread(
      f,
      nrows = nrows,
      data.table = FALSE,
      showProgress = FALSE,
      fill = TRUE,
      check.names = FALSE
    )

    as.data.frame(d)

  }, error = function(e) {

    NULL

  })
}

# ============================================================
# Find columns
# ============================================================

find_col <- function(df, patterns, exclude = NULL) {

  if (is.null(df) || ncol(df) == 0)
    return(NA_character_)

  nn <- clean_name(names(df))

  hit <- integer()

  for (p in patterns) {

    x <- grep(p, nn, perl = TRUE)

    if (length(x) > 0) {
      hit <- c(hit, x)
    }
  }

  hit <- unique(hit)

  if (!is.null(exclude)) {

    bad <- integer()

    for (p in exclude) {
      bad <- c(bad, grep(p, nn, perl = TRUE))
    }

    hit <- setdiff(hit, bad)
  }

  if (length(hit) == 0)
    return(NA_character_)

  names(df)[hit[1]]
}

get_gene_col <- function(df) {

  find_col(
    df,
    c(
      "^genename$",
      "^genesymbol$",
      "^symbol$",
      "^gene$",
      "^geneid$",
      "gene"
    )
  )
}

get_event_col <- function(df) {

  find_col(
    df,
    c(
      "^eventid$",
      "^event$",
      "^siteid$",
      "^site$",
      "^transcriptid$",
      "^transcript$",
      "polyasite",
      "pasid"
    )
  )
}

get_p_col <- function(df) {

  find_col(
    df,
    c(
      "^p$",
      "^pvalue$",
      "^pval$",
      "^pvaluewelch$",
      "^welchp$",
      "^ttestp$",
      "^nominalp$",
      "pvalue",
      "pval"
    ),
    exclude = c(
      "adjust",
      "padj",
      "fdr",
      "qvalue"
    )
  )
}

get_fdr_col <- function(df) {

  find_col(
    df,
    c(
      "^fdr$",
      "^padj$",
      "^adjustedp$",
      "^adjustedpvalue$",
      "^adjp$",
      "^qvalue$",
      "fdr",
      "padj",
      "qvalue",
      "adjustedp"
    )
  )
}

get_delta_col <- function(df) {

  find_col(
    df,
    c(
      "^deltapau$",
      "^dpau$",
      "deltapau",
      "delta.*pau",
      "paudelta",
      "mean.*difference",
      "difference.*pau"
    )
  )
}

get_direction_col <- function(df) {

  find_col(
    df,
    c(
      "^direction$",
      "^apadirection$",
      "^classification$",
      "^mechanism$",
      "lengthdirection",
      "utr.*direction",
      "apaclass"
    )
  )
}

get_length_delta_col <- function(df) {

  find_col(
    df,
    c(
      "delt.*weighted.*length",
      "weighted.*length.*delta",
      "deltalength",
      "delta3utr",
      "utr.*length.*change",
      "lengthchange"
    )
  )
}

# ============================================================
# Detect sample columns in a PAU matrix
# ============================================================

detect_sample_columns <- function(df) {

  if (is.null(df) || ncol(df) < 3)
    return(character())

  annotation_regex <- paste(
    c(
      "^gene$",
      "geneid",
      "genename",
      "symbol",
      "chrom",
      "^chr$",
      "strand",
      "start",
      "end",
      "position",
      "coordinate",
      "site",
      "transcript",
      "length",
      "numsite",
      "numberofsite",
      "event",
      "pas",
      "utr",
      "pvalue",
      "pval",
      "padj",
      "fdr",
      "qvalue",
      "delta",
      "direction",
      "mean",
      "median",
      "log2",
      "fold",
      "ratio",
      "status",
      "class"
    ),
    collapse = "|"
  )

  out <- character()

  for (nm in names(df)) {

    cn <- clean_name(nm)

    if (grepl(annotation_regex, cn, perl = TRUE))
      next

    x <- safe_numeric(df[[nm]])

    valid <- is.finite(x)

    if (mean(valid) < 0.70)
      next

    xx <- x[valid]

    if (length(xx) < 3)
      next

    # PAU is normally 0-1 or 0-100.
    # This helps avoid genomic coordinates.
    frac_reasonable <- mean(xx >= -0.001 & xx <= 100.001)

    if (frac_reasonable >= 0.95) {
      out <- c(out, nm)
    }
  }

  unique(out)
}

# ============================================================
# Infer Control / Experimental labels
# ============================================================

infer_group <- function(sample_name, dataset = NULL) {

  x <- tolower(sample_name)

  # Remove .PAU suffix safely
  x <- sub("[.]pau$", "", x, ignore.case = TRUE)

  # ==========================================================
  # Dataset-specific rules
  # ==========================================================

  # TDP-43 K263E
  if (!is.null(dataset) && dataset == "TDP43_K263E") {

    if (grepl("^wt", x))
      return("Control")

    if (grepl("k263e", x))
      return("Experimental")
  }

  # FUS P525L
  if (!is.null(dataset) && dataset == "FUS_P525L") {

    if (grepl("^wt", x))
      return("Control")

    if (grepl("p525l", x))
      return("Experimental")
  }

  # MATR3 KO
  if (!is.null(dataset) && dataset == "MATR3_KO") {

    if (grepl("^wt", x))
      return("Control")

    if (grepl("^matr3_?ko", x))
      return("Experimental")
  }

  # TDP-43 KO
  # Cont_S1-S4 = control
  # TDP43_S5-S7 = experimental
  if (!is.null(dataset) && dataset == "TDP43_KO") {

    if (grepl("^cont", x))
      return("Control")

    if (grepl("^tdp43", x))
      return("Experimental")
  }

  # FUS KO
  # WT = control
  # A4 and A5 = FUS KO clones
  if (!is.null(dataset) && dataset == "FUS_KO") {

    if (grepl("^wt", x))
      return("Control")

    if (grepl("^a4", x) || grepl("^a5", x))
      return("Experimental")
  }

  # ==========================================================
  # Generic fallback
  # ==========================================================

  if (
    grepl(
      "(^|[_.-])(ctrl|control|cont|con|ref|wt|wildtype|normal|mock|vehicle)([_.-]|$)",
      x,
      perl = TRUE
    )
  ) {
    return("Control")
  }

  if (
    grepl(
      "^(ctrl|control|cont|con|ref|wt|normal)[_-]?[0-9]",
      x,
      perl = TRUE
    )
  ) {
    return("Control")
  }

  if (
    grepl(
      "(ko|kd|knockout|knockdown|k263e|p525l|mut|mutant|patient|case)",
      x,
      perl = TRUE
    )
  ) {
    return("Experimental")
  }

  NA_character_
}

# ============================================================
# Direction classification
# ============================================================

classify_direction <- function(df) {

  n <- nrow(df)

  ans <- rep(NA_character_, n)

  dcol <- get_direction_col(df)

  if (!is.na(dcol)) {

    x <- tolower(as.character(df[[dcol]]))

    ans[grepl("short", x)] <- "Shortening"

    ans[
      grepl("length|long", x)
    ] <- "Lengthening"

    return(ans)
  }

  # If there is a PAU-weighted 3'UTR length delta,
  # positive = lengthening
  # negative = shortening
  lcol <- get_length_delta_col(df)

  if (!is.na(lcol)) {

    x <- safe_numeric(df[[lcol]])

    ans[x > 0] <- "Lengthening"
    ans[x < 0] <- "Shortening"

    return(ans)
  }

  # IMPORTANT:
  # Do NOT classify shortening/lengthening from deltaPAU sign alone.
  # Proximal and distal sites have opposite biological meaning.

  ans
}

# ============================================================
# Score differential-statistics candidate files
# ============================================================

score_stats_file <- function(f) {

  d <- read_table_safe(f, nrows = 100)

  if (is.null(d))
    return(-Inf)

  s <- 0

  if (!is.na(get_p_col(d)))         s <- s + 8
  if (!is.na(get_fdr_col(d)))       s <- s + 6
  if (!is.na(get_delta_col(d)))     s <- s + 8
  if (!is.na(get_gene_col(d)))      s <- s + 4
  if (!is.na(get_event_col(d)))     s <- s + 2
  if (!is.na(get_direction_col(d))) s <- s + 4
  if (!is.na(get_length_delta_col(d))) s <- s + 4

  base <- tolower(basename(f))

  if (grepl("result|stat|differ|volcano|mechanistic|apa", base))
    s <- s + 3

  if (grepl("pca", base))
    s <- s - 10

  s
}

# ============================================================
# Score PAU matrix candidate
# ============================================================

score_pau_file <- function(f) {

  d <- read_table_safe(f, nrows = 250)

  if (is.null(d))
    return(-Inf)

  sample_cols <- detect_sample_columns(d)

  s <- length(sample_cols)

  base <- tolower(basename(f))

  if (grepl("pau", base))
    s <- s + 10

  if (grepl("quant", base))
    s <- s + 4

  if (grepl("matrix", base))
    s <- s + 3

  if (grepl("pca", base))
    s <- s - 10

  s
}

# ============================================================
# PCA calculation
# ============================================================

calculate_pca <- function(df, sample_cols) {

  ans <- list(
    PC1_scaled   = NA_real_,
    PC2_scaled   = NA_real_,
    PC1_unscaled = NA_real_,
    PC2_unscaled = NA_real_,
    n_events_PCA = NA_integer_
  )

  if (
    is.null(df) ||
    length(sample_cols) < 3
  ) {
    return(ans)
  }

  mat <- sapply(
    sample_cols,
    function(x) safe_numeric(df[[x]])
  )

  mat <- as.matrix(mat)

  # Rows = APA events
  # Columns = samples

  keep <- apply(
    mat,
    1,
    function(x) {
      sum(is.finite(x)) >= max(2, ceiling(ncol(mat) / 2))
    }
  )

  mat <- mat[keep, , drop = FALSE]

  if (nrow(mat) < 3)
    return(ans)

  # Row-mean imputation
  for (i in seq_len(nrow(mat))) {

    miss <- !is.finite(mat[i, ])

    if (any(miss)) {

      m <- mean(mat[i, !miss], na.rm = TRUE)

      mat[i, miss] <- m
    }
  }

  # Remove zero-variance APA events
  sds <- apply(mat, 1, sd, na.rm = TRUE)

  mat <- mat[
    is.finite(sds) & sds > 0,
    ,
    drop = FALSE
  ]

  ans$n_events_PCA <- nrow(mat)

  if (nrow(mat) < 3)
    return(ans)

  # --------------------------------------------
  # Scaled PCA
  # --------------------------------------------

  pca1 <- tryCatch(
    prcomp(
      t(mat),
      center = TRUE,
      scale. = TRUE
    ),
    error = function(e) NULL
  )

  if (!is.null(pca1)) {

    vv <- pca1$sdev^2 / sum(pca1$sdev^2) * 100

    if (length(vv) >= 1)
      ans$PC1_scaled <- vv[1]

    if (length(vv) >= 2)
      ans$PC2_scaled <- vv[2]
  }

  # --------------------------------------------
  # Unscaled PCA
  # --------------------------------------------

  pca2 <- tryCatch(
    prcomp(
      t(mat),
      center = TRUE,
      scale. = FALSE
    ),
    error = function(e) NULL
  )

  if (!is.null(pca2)) {

    vv <- pca2$sdev^2 / sum(pca2$sdev^2) * 100

    if (length(vv) >= 1)
      ans$PC1_unscaled <- vv[1]

    if (length(vv) >= 2)
      ans$PC2_unscaled <- vv[2]
  }

  ans
}

# ============================================================
# Recursively find all tables
# ============================================================

cat("\n============================================================\n")
cat("Scanning:\n", ROOT, "\n")
cat("============================================================\n\n")

all_files <- list.files(
  ROOT,
  recursive = TRUE,
  full.names = TRUE
)

all_files <- all_files[
  grepl(
    "\\.(csv|tsv|txt|xlsx|xls)$",
    all_files,
    ignore.case = TRUE
  )
]

# Do not scan our own new output
all_files <- all_files[
  !grepl(
    "Revision_quantitative_summary",
    all_files,
    fixed = TRUE
  )
]

cat("Tabular files found:", length(all_files), "\n\n")

# ============================================================
# Assign files to datasets
# ============================================================

inventory <- data.frame(
  file = all_files,
  dataset = vapply(
    all_files,
    dataset_from_path,
    character(1)
  ),
  stringsAsFactors = FALSE
)

# ============================================================
# Score files
# ============================================================

inventory$stats_score <- NA_real_
inventory$pau_score   <- NA_real_

for (i in seq_len(nrow(inventory))) {

  f <- inventory$file[i]

  if (is.na(inventory$dataset[i]))
    next

  inventory$stats_score[i] <- score_stats_file(f)
  inventory$pau_score[i]   <- score_pau_file(f)
}

data.table::fwrite(
  inventory,
  file.path(
    OUT,
    "00_all_candidate_files.csv"
  )
)

# ============================================================
# Output objects
# ============================================================

summary_list <- list()
top_list     <- list()
sample_list  <- list()
warning_list <- character()

# ============================================================
# Dataset loop
# ============================================================

for (dataset in names(dataset_defs)) {

  cat("\n============================================================\n")
  cat("DATASET:", dataset, "\n")
  cat("============================================================\n")

  inv <- inventory[
    inventory$dataset == dataset &
      !is.na(inventory$dataset),
    ,
    drop = FALSE
  ]

  if (nrow(inv) == 0) {

    cat("No files detected.\n")

    warning_list <- c(
      warning_list,
      paste(dataset, ": no matching files detected")
    )

    next
  }

  # ----------------------------------------------------------
  # Best statistics file
  # ----------------------------------------------------------

  stats_file <- NA_character_
  stats_df   <- NULL

  ok_stats <- is.finite(inv$stats_score)

  if (any(ok_stats)) {

    z <- inv[ok_stats, , drop = FALSE]

    z <- z[
      order(z$stats_score, decreasing = TRUE),
      ,
      drop = FALSE
    ]

    if (
      nrow(z) > 0 &&
      z$stats_score[1] >= 10
    ) {

      stats_file <- z$file[1]

      stats_df <- read_table_safe(stats_file)

      cat(
        "Statistics file:\n  ",
        stats_file,
        "\n"
      )

      cat(
        "Stats score:",
        z$stats_score[1],
        "\n"
      )
    }
  }

  # ----------------------------------------------------------
  # Best PAU matrix
  # ----------------------------------------------------------

  pau_file <- NA_character_
  pau_df   <- NULL
  sample_cols <- character()

  ok_pau <- is.finite(inv$pau_score)

  if (any(ok_pau)) {

    z <- inv[ok_pau, , drop = FALSE]

    z <- z[
      order(z$pau_score, decreasing = TRUE),
      ,
      drop = FALSE
    ]

    if (
      nrow(z) > 0 &&
      z$pau_score[1] >= 3
    ) {

      pau_file <- z$file[1]

      pau_df <- read_table_safe(pau_file)

      sample_cols <- detect_sample_columns(pau_df)

      cat(
        "PAU matrix:\n  ",
        pau_file,
        "\n"
      )

      cat(
        "Detected sample columns:",
        length(sample_cols),
        "\n"
      )
    }
  }

  # ----------------------------------------------------------
  # Sample counts
  # ----------------------------------------------------------

  n_control      <- NA_integer_
  n_experimental <- NA_integer_
  n_samples      <- NA_integer_

  if (length(sample_cols) > 0) {

    groups <- vapply(
      sample_cols,
      function(x) infer_group(x, dataset),
      character(1)
    )

    sample_table <- data.frame(
      Dataset = dataset,
      Sample = sample_cols,
      Inferred_group = groups,
      PAU_file = pau_file,
      stringsAsFactors = FALSE
    )

    sample_list[[dataset]] <- sample_table

    n_samples <- length(sample_cols)

    if (any(groups == "Control", na.rm = TRUE)) {
      n_control <- sum(groups == "Control", na.rm = TRUE)
    }

    if (any(groups == "Experimental", na.rm = TRUE)) {
      n_experimental <- sum(
        groups == "Experimental",
        na.rm = TRUE
      )
    }

    cat("\nSamples:\n")

    print(sample_table[, c(
      "Sample",
      "Inferred_group"
    )])

    cat(
      "\nControl n =",
      n_control,
      "\nExperimental n =",
      n_experimental,
      "\n"
    )

    if (any(is.na(groups))) {

      warning_list <- c(
        warning_list,
        paste(
          dataset,
          ": some sample groups could not be inferred:",
          paste(
            sample_cols[is.na(groups)],
            collapse = ", "
          )
        )
      )
    }
  }

  # ----------------------------------------------------------
  # Detected APA events / genes
  # Prefer PAU matrix for total detected counts
  # ----------------------------------------------------------

  detected_events <- NA_integer_
  detected_genes  <- NA_integer_

  gene_col_pau <- NA_character_

  if (!is.null(pau_df)) {

    detected_events <- nrow(pau_df)

    gene_col_pau <- get_gene_col(pau_df)

    if (!is.na(gene_col_pau)) {

      g <- as.character(
        pau_df[[gene_col_pau]]
      )

      g <- g[
        !is.na(g) &
        nzchar(g)
      ]

      detected_genes <- length(unique(g))
    }
  }

  # ----------------------------------------------------------
  # Differential stats
  # ----------------------------------------------------------

  n_tested_events <- NA_integer_
  n_tested_genes  <- NA_integer_

  nominal_events <- NA_integer_
  nominal_genes  <- NA_integer_

  fdr_events <- NA_integer_
  fdr_genes  <- NA_integer_

  nominal_short <- NA_integer_
  nominal_long  <- NA_integer_

  fdr_short <- NA_integer_
  fdr_long  <- NA_integer_

  median_dpau_all <- NA_real_
  median_abs_dpau_all <- NA_real_

  median_dpau_nominal <- NA_real_
  median_abs_dpau_nominal <- NA_real_

  median_dpau_fdr <- NA_real_
  median_abs_dpau_fdr <- NA_real_

  max_abs_dpau <- NA_real_

  pcol     <- NA_character_
  fdrcol   <- NA_character_
  dcol     <- NA_character_
  genecol  <- NA_character_
  eventcol <- NA_character_

  if (!is.null(stats_df)) {

    n_tested_events <- nrow(stats_df)

    pcol     <- get_p_col(stats_df)
    fdrcol   <- get_fdr_col(stats_df)
    dcol     <- get_delta_col(stats_df)
    genecol  <- get_gene_col(stats_df)
    eventcol <- get_event_col(stats_df)

    # --------------------------------------------------------
    # If FDR missing but P exists, calculate BH FDR
    # --------------------------------------------------------

    if (
      is.na(fdrcol) &&
      !is.na(pcol)
    ) {

      pp <- safe_numeric(stats_df[[pcol]])

      stats_df$FDR_BH_calculated <- p.adjust(
        pp,
        method = "BH"
      )

      fdrcol <- "FDR_BH_calculated"

      cat(
        "FDR column absent; calculated BH-adjusted FDR.\n"
      )
    }

    if (!is.na(genecol)) {

      gg <- as.character(
        stats_df[[genecol]]
      )

      gg <- gg[
        !is.na(gg) &
        nzchar(gg)
      ]

      n_tested_genes <- length(unique(gg))
    }

    # --------------------------------------------------------
    # P / FDR vectors
    # --------------------------------------------------------

    pvec <- rep(NA_real_, nrow(stats_df))
    qvec <- rep(NA_real_, nrow(stats_df))
    dvec <- rep(NA_real_, nrow(stats_df))

    if (!is.na(pcol))
      pvec <- safe_numeric(stats_df[[pcol]])

    if (!is.na(fdrcol))
      qvec <- safe_numeric(stats_df[[fdrcol]])

    if (!is.na(dcol))
      dvec <- safe_numeric(stats_df[[dcol]])

    nominal_hit <- is.finite(pvec) & pvec < 0.05
    fdr_hit     <- is.finite(qvec) & qvec < 0.05

    if (!is.na(pcol)) {

      nominal_events <- sum(nominal_hit)

      if (!is.na(genecol)) {

        nominal_genes <- length(
          unique(
            as.character(
              stats_df[[genecol]][nominal_hit]
            )
          )
        )
      }
    }

    if (!is.na(fdrcol)) {

      fdr_events <- sum(fdr_hit)

      if (!is.na(genecol)) {

        fdr_genes <- length(
          unique(
            as.character(
              stats_df[[genecol]][fdr_hit]
            )
          )
        )
      }
    }

    # --------------------------------------------------------
    # Delta PAU
    # --------------------------------------------------------

    if (!is.na(dcol)) {

      valid <- is.finite(dvec)

      median_dpau_all <- median(
        dvec[valid],
        na.rm = TRUE
      )

      median_abs_dpau_all <- median(
        abs(dvec[valid]),
        na.rm = TRUE
      )

      if (any(nominal_hit & valid)) {

        median_dpau_nominal <- median(
          dvec[nominal_hit & valid],
          na.rm = TRUE
        )

        median_abs_dpau_nominal <- median(
          abs(dvec[nominal_hit & valid]),
          na.rm = TRUE
        )
      }

      if (any(fdr_hit & valid)) {

        median_dpau_fdr <- median(
          dvec[fdr_hit & valid],
          na.rm = TRUE
        )

        median_abs_dpau_fdr <- median(
          abs(dvec[fdr_hit & valid]),
          na.rm = TRUE
        )
      }

      if (any(valid)) {
        max_abs_dpau <- max(
          abs(dvec[valid]),
          na.rm = TRUE
        )
      }
    }

    # --------------------------------------------------------
    # Shortening / lengthening
    # --------------------------------------------------------

    direction <- classify_direction(stats_df)

    if (any(!is.na(direction))) {

      if (!is.na(pcol)) {

        nominal_short <- sum(
          nominal_hit &
          direction == "Shortening",
          na.rm = TRUE
        )

        nominal_long <- sum(
          nominal_hit &
          direction == "Lengthening",
          na.rm = TRUE
        )
      }

      if (!is.na(fdrcol)) {

        fdr_short <- sum(
          fdr_hit &
          direction == "Shortening",
          na.rm = TRUE
        )

        fdr_long <- sum(
          fdr_hit &
          direction == "Lengthening",
          na.rm = TRUE
        )
      }

    } else {

      warning_list <- c(
        warning_list,
        paste(
          dataset,
          ": shortening/lengthening could NOT be derived safely.",
          "No direction or weighted-3UTR-length-change column was found.",
          "DeltaPAU sign alone was deliberately NOT used."
        )
      )
    }

    # --------------------------------------------------------
    # Top representative events
    # --------------------------------------------------------

    if (!is.na(dcol)) {

      top <- stats_df

      top$.__delta__ <- dvec
      top$.__abs_delta__ <- abs(dvec)
      top$.__p__ <- pvec
      top$.__fdr__ <- qvec
      top$.__direction__ <- direction

      top <- top[
        is.finite(top$.__abs_delta__),
        ,
        drop = FALSE
      ]

      top <- top[
        order(
          top$.__abs_delta__,
          decreasing = TRUE
        ),
        ,
        drop = FALSE
      ]

      top <- head(top, 20)

      keep <- c()

      if (!is.na(genecol))
        keep <- c(keep, genecol)

      if (!is.na(eventcol))
        keep <- c(keep, eventcol)

      keep <- unique(c(
        keep,
        ".__delta__",
        ".__abs_delta__",
        ".__p__",
        ".__fdr__",
        ".__direction__"
      ))

      top2 <- top[, keep, drop = FALSE]

      names(top2)[names(top2) == ".__delta__"] <-
        "DeltaPAU"

      names(top2)[names(top2) == ".__abs_delta__"] <-
        "Abs_DeltaPAU"

      names(top2)[names(top2) == ".__p__"] <-
        "P_value"

      names(top2)[names(top2) == ".__fdr__"] <-
        "FDR"

      names(top2)[names(top2) == ".__direction__"] <-
        "APA_direction"

      top2$Dataset <- dataset

      top_list[[dataset]] <- top2
    }
  }

  # ----------------------------------------------------------
  # PCA
  # ----------------------------------------------------------

  pca <- calculate_pca(
    pau_df,
    sample_cols
  )

  cat(
    "\nPCA scaled:\n",
    "PC1 =",
    round(pca$PC1_scaled, 2),
    "%\n",
    "PC2 =",
    round(pca$PC2_scaled, 2),
    "%\n"
  )

  cat(
    "\nPCA unscaled:\n",
    "PC1 =",
    round(pca$PC1_unscaled, 2),
    "%\n",
    "PC2 =",
    round(pca$PC2_unscaled, 2),
    "%\n"
  )

  # ----------------------------------------------------------
  # Final summary row
  # ----------------------------------------------------------

  summary_list[[dataset]] <- data.frame(

    Dataset = dataset,

    Control_n = n_control,
    Experimental_n = n_experimental,
    Total_samples = n_samples,

    APA_events_detected = detected_events,
    APA_genes_detected = detected_genes,

    Events_tested = n_tested_events,
    Genes_tested = n_tested_genes,

    Nominal_P_lt_0.05_events = nominal_events,
    Nominal_P_lt_0.05_genes = nominal_genes,

    FDR_lt_0.05_events = fdr_events,
    FDR_lt_0.05_genes = fdr_genes,

    Nominal_shortening = nominal_short,
    Nominal_lengthening = nominal_long,

    FDR_shortening = fdr_short,
    FDR_lengthening = fdr_long,

    Median_DeltaPAU_all = median_dpau_all,
    Median_abs_DeltaPAU_all = median_abs_dpau_all,

    Median_DeltaPAU_nominal = median_dpau_nominal,
    Median_abs_DeltaPAU_nominal =
      median_abs_dpau_nominal,

    Median_DeltaPAU_FDR = median_dpau_fdr,
    Median_abs_DeltaPAU_FDR =
      median_abs_dpau_fdr,

    Max_abs_DeltaPAU = max_abs_dpau,

    PCA_PC1_scaled_percent =
      pca$PC1_scaled,

    PCA_PC2_scaled_percent =
      pca$PC2_scaled,

    PCA_PC1_unscaled_percent =
      pca$PC1_unscaled,

    PCA_PC2_unscaled_percent =
      pca$PC2_unscaled,

    PCA_events_used =
      pca$n_events_PCA,

    Statistics_file =
      ifelse(
        is.na(stats_file),
        "",
        stats_file
      ),

    PAU_file =
      ifelse(
        is.na(pau_file),
        "",
        pau_file
      ),

    P_column =
      ifelse(
        is.na(pcol),
        "",
        pcol
      ),

    FDR_column =
      ifelse(
        is.na(fdrcol),
        "",
        fdrcol
      ),

    DeltaPAU_column =
      ifelse(
        is.na(dcol),
        "",
        dcol
      ),

    stringsAsFactors = FALSE
  )
}

# ============================================================
# Combine outputs
# ============================================================

if (length(summary_list) > 0) {

  summary_df <- data.table::rbindlist(
    summary_list,
    fill = TRUE
  )

  summary_df <- as.data.frame(summary_df)

  data.table::fwrite(
    summary_df,
    file.path(
      OUT,
      "01_revision_quantitative_summary.csv"
    )
  )

  cat("\n\n")
  cat("============================================================\n")
  cat("FINAL QUANTITATIVE SUMMARY\n")
  cat("============================================================\n\n")

  print(
    summary_df,
    row.names = FALSE
  )
}

# ============================================================
# Samples
# ============================================================

if (length(sample_list) > 0) {

  sample_df <- data.table::rbindlist(
    sample_list,
    fill = TRUE
  )

  data.table::fwrite(
    sample_df,
    file.path(
      OUT,
      "02_detected_samples_and_groups.csv"
    )
  )
}

# ============================================================
# Top representative events
# ============================================================

if (length(top_list) > 0) {

  top_df <- data.table::rbindlist(
    top_list,
    fill = TRUE
  )

  data.table::fwrite(
    top_df,
    file.path(
      OUT,
      "03_top_representative_DeltaPAU_events.csv"
    )
  )
}

# ============================================================
# Warnings
# ============================================================

if (length(warning_list) == 0) {

  warning_list <- "No warnings."

} else {

  warning_list <- unique(warning_list)
}

writeLines(
  warning_list,
  file.path(
    OUT,
    "04_warnings.txt"
  )
)

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")

cat(
  "\nOutput directory:\n",
  OUT,
  "\n\n"
)

cat(
  "Main file:\n",
  file.path(
    OUT,
    "01_revision_quantitative_summary.csv"
  ),
  "\n\n"
)

