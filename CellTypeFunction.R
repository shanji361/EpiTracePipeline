# ============================================================================
# EpiTrace Cell-Type Rebalancing — Unified Three-Regime Perturbation Operator
#
# Purpose : Controlled perturbation of class distributions to measure
#           EpiTrace robustness — NOT biological correction.
#
# Regimes :
#   "down"  — Regime A  Pure downsampling. Majority classes lose cells.
#                        Rare classes are untouched. Total N shrinks.
#   "mixed" — Regime B  Convergence toward a geometric mean target.
#                        Large classes shrink; small classes are duplicated
#                        (with replacement) up to a moderate ceiling.
#                        Total N stays close to original.
#   "up"    — Regime C  Stress-test. All classes are upsampled toward the
#                        largest class size (with replacement). Rare classes
#                        are heavily duplicated. Total N grows.
#
# Core formula :
#   target_c = f(n_c, alpha, mode)
#   alpha in [0, 1] — strength of perturbation
#   alpha = 0 always returns the original distribution regardless of mode
#
# new_to_orig convention (used throughout):
#   names(new_to_orig) = new barcode (unique, used as Seurat cell name)
#   unname(new_to_orig) = original barcode (used to slice assay matrices)
# ============================================================================


# ── Internal utilities ───────────────────────────────────────────────────────

.gini <- function(counts) {
  v <- sort(counts[counts > 0])
  n <- length(v)
  2 * sum(seq_len(n) * v) / (n * sum(v)) - (n + 1) / n
}

.print_summary <- function(alpha, mode, n_c, targets, classes, seed) {
  summary_df <- data.frame(
    celltype   = classes,
    n_before   = n_c[classes],
    n_after    = targets[classes],
    pct_change = round(100 * (targets[classes] / n_c[classes] - 1), 1),
    row.names  = NULL
  )
  summary_df <- summary_df[order(-summary_df$n_before), ]
  
  message(sprintf(
    "\n── resample_cells  alpha=%.2f  mode='%s'  seed=%d ──",
    alpha, mode, seed))
  message(sprintf("  Total cells : %d -> %d",
                  sum(n_c), sum(targets)))
  message(sprintf("  Gini  before: %.3f   after: %.3f",
                  .gini(n_c), .gini(targets)))
  message(sprintf("  Ratio before: %.1fx  after: %.1fx",
                  max(n_c) / min(n_c),
                  max(targets) / min(targets)))
  print(summary_df, row.names = FALSE)
}


# ── Target-size calculators, one per regime ──────────────────────────────────

.targets_down <- function(n_c, alpha) {
  n_min   <- min(n_c)
  targets <- round(n_c^(1 - alpha) * n_min^alpha)
  pmin(targets, n_c)
}

.targets_mixed <- function(n_c, alpha, cap_multiplier = 3) {
  G       <- exp(mean(log(n_c)))
  targets <- round(n_c^(1 - alpha) * G^alpha)
  cap     <- round(min(n_c) * cap_multiplier)
  ifelse(targets > n_c, pmin(targets, cap), targets)
}

.targets_up <- function(n_c, alpha) {
  n_max   <- max(n_c)
  targets <- round(n_c^(1 - alpha) * n_max^alpha)
  pmax(targets, n_c)
}


# ── Main function ─────────────────────────────────────────────────────────────

resample_cells <- function(seurat_obj,
                           alpha          = 0.5,
                           mode           = c("down", "mixed", "up"),
                           celltype_col   = "celltype",
                           cap_multiplier = 3,
                           seed           = 1234,
                           verbose        = TRUE) {
  
  mode <- match.arg(mode)
  stopifnot(
    is.numeric(alpha), length(alpha) == 1, alpha >= 0, alpha <= 1,
    celltype_col %in% colnames(seurat_obj@meta.data),
    is.numeric(cap_multiplier), cap_multiplier >= 1
  )
  
  set.seed(seed)
  
  meta      <- seurat_obj@meta.data
  celltypes <- as.character(meta[[celltype_col]])
  all_bcs   <- rownames(meta)
  
  ct_tab  <- table(celltypes)
  classes <- names(ct_tab)
  n_c     <- as.integer(ct_tab)
  names(n_c) <- classes
  
  targets <- switch(mode,
                    down  = .targets_down(n_c, alpha),
                    mixed = .targets_mixed(n_c, alpha, cap_multiplier),
                    up    = .targets_up(n_c, alpha)
  )
  names(targets) <- classes
  
  # new_to_orig: names = new unique barcode, values = original source barcode
  # For unmodified cells both are identical.
  # For duplicates: name = "BARCODE_dup1", value = "BARCODE"
  new_to_orig <- character(0)
  
  for (ct in classes) {
    pool <- all_bcs[celltypes == ct]
    tgt  <- targets[ct]
    
    if (tgt <= length(pool)) {
      # Downsampling — no duplication, names and values are the same barcode
      drawn       <- sample(pool, size = tgt, replace = FALSE)
      new_to_orig <- c(new_to_orig, setNames(drawn, drawn))
      
    } else {
      # Upsampling — keep every original cell, then add duplicates
      gap   <- tgt - length(pool)
      extra <- sample(pool, size = gap, replace = TRUE)
      
      # FIX 1: dup_counter is a plain list — no rlang %||% needed.
      # FIX 2: counter starts at 1L directly; avoids 0L + 1L arithmetic
      #         that produced NA when the list lookup returned NULL.
      dup_counter <- list()
      extra_new   <- character(gap)  # new unique barcodes for duplicates
      extra_orig  <- character(gap)  # corresponding original barcodes
      
      for (j in seq_len(gap)) {
        src  <- extra[j]
        cnt  <- dup_counter[[src]]
        cnt  <- if (is.null(cnt) || is.na(cnt)) 1L else cnt + 1L
        dup_counter[[src]] <- cnt
        extra_new[j]  <- sprintf("%s_dup%d", src, cnt)  # new unique name
        extra_orig[j] <- src                              # original source
      }
      
      # FIX 3: setNames(value, name) — value must be the original barcode
      #         so the matrix slicer can find it; name must be the new unique
      #         barcode so Seurat rownames stay unique.
      #         Old code had these backwards: setNames(extra, extra_new)
      #         which put the original barcodes as names and new as values,
      #         causing the matrix slice to use _dupN as a column index
      #         (which doesn't exist) and rownames to get the original
      #         barcode (causing duplicates).
      new_to_orig <- c(new_to_orig,
                       setNames(pool,       pool),       # identity
                       setNames(extra_orig, extra_new))  # value=orig, name=new
    }
  }
  
  # ── Fast path: mode="down", no duplicates ────────────────────────────────
  # subset() propagates ALL assays atomically — RNA and ATAC are guaranteed
  # to cover the same balanced cell set.
  needs_dup <- any(duplicated(names(new_to_orig)))
  
  if (!needs_dup) {
    result_obj <- subset(seurat_obj, cells = names(new_to_orig))
    result_obj$original_cell  <- names(new_to_orig)[
      match(rownames(result_obj@meta.data), names(new_to_orig))]
    result_obj$resample_alpha <- alpha
    result_obj$resample_mode  <- mode
    if (verbose) .print_summary(alpha, mode, n_c, targets, classes, seed)
    return(result_obj)
  }
  
  # ── Slow path: mode="mixed" or "up", duplicated cells ───────────────────
  # Both ATAC and RNA are sliced with the same new_to_orig vector so the
  # two modalities cover exactly the same balanced cell set.
  
  orig_unique <- unique(unname(new_to_orig))       # original barcodes needed
  base_obj    <- subset(seurat_obj, cells = orig_unique)
  
  # Build metadata for the new (possibly duplicated) cell set.
  # Index orig_meta by unname(new_to_orig) (original barcodes) to get the
  # right row for each new cell, then rename rows to names(new_to_orig).
  orig_meta          <- base_obj@meta.data
  new_meta           <- orig_meta[unname(new_to_orig), , drop = FALSE]
  rownames(new_meta) <- names(new_to_orig)         # unique new barcodes
  new_meta$original_cell  <- unname(new_to_orig)
  new_meta$resample_alpha <- alpha
  new_meta$resample_mode  <- mode
  
  # Transfer assays — skip any that have no counts layer (e.g. chromvar)
  assay_names <- Seurat::Assays(base_obj)
  
  has_counts <- vapply(assay_names, function(an) {
    mtx <- tryCatch(
      Seurat::GetAssayData(base_obj, assay = an, layer = "counts"),
      error = function(e) NULL
    )
    !is.null(mtx) && nrow(mtx) > 0 && ncol(mtx) > 0
  }, logical(1))
  
  assay_names_ok <- assay_names[has_counts]
  
  if (any(!has_counts)) {
    message("  Note: skipping assays with no counts layer (not transferred): ",
            paste(assay_names[!has_counts], collapse = ", "))
  }
  
  # Slice every valid assay using unname(new_to_orig) as column indices
  # (original barcodes), then rename columns to names(new_to_orig) (new
  # unique barcodes). RNA and ATAC get identical column sets — aligned.
  new_assay_list <- lapply(assay_names_ok, function(an) {
    mtx     <- Seurat::GetAssayData(base_obj, assay = an, layer = "counts")
    new_mtx <- mtx[, unname(new_to_orig), drop = FALSE]  # slice by orig barcode
    colnames(new_mtx) <- names(new_to_orig)               # rename to new barcode
    Seurat::CreateAssayObject(counts = new_mtx,
                              min.cells = 0, min.features = 0,
                              check.matrix = FALSE)
  })
  names(new_assay_list) <- assay_names_ok
  
  result_obj <- Seurat::CreateSeuratObject(
    counts    = Seurat::GetAssayData(new_assay_list[[1]], layer = "counts"),
    meta.data = new_meta,
    min.cells = 0, min.features = 0
  )
  if (length(assay_names_ok) > 1) {
    for (an in assay_names_ok[-1]) result_obj[[an]] <- new_assay_list[[an]]
  }
  
  if (length(base_obj@reductions) > 0) {
    message("  Note: reductions not transferred — ",
            "duplicated cells have no unique embedding.")
  }
  
  if (verbose) .print_summary(alpha, mode, n_c, targets, classes, seed)
  return(result_obj)
}
