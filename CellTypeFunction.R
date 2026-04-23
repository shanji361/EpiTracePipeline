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
  message(sprintf("  Total cells : %d → %d",
                  sum(n_c), sum(targets)))
  message(sprintf("  Gini  before: %.3f   after: %.3f",
                  .gini(n_c), .gini(targets)))
  message(sprintf("  Ratio before: %.1fx  after: %.1fx",
                  max(n_c) / min(n_c),
                  max(targets) / min(targets)))
  print(summary_df, row.names = FALSE)
}


# ── Target-size calculators, one per regime ──────────────────────────────────

# Regime A — Downsampling only
# target_c = n_c^(1-alpha) * n_min^alpha
# alpha=0 → target_c = n_c  (no change)
# alpha=1 → target_c = n_min for all (everything floored to the rarest)
.targets_down <- function(n_c, alpha) {
  n_min   <- min(n_c)
  targets <- round(n_c^(1 - alpha) * n_min^alpha)
  pmin(targets, n_c)    # strict downsampling: never exceed original
}

# Regime B — Mixed (convergence to geometric mean)
# G = geometric mean of all class sizes
# target_c = n_c^(1-alpha) * G^alpha
# alpha=0 → target_c = n_c
# alpha=1 → target_c = G for all classes
# Classes above G shrink; classes below G grow (capped at cap_multiplier * n_min)
.targets_mixed <- function(n_c, alpha, cap_multiplier = 3) {
  G       <- exp(mean(log(n_c)))
  targets <- round(n_c^(1 - alpha) * G^alpha)
  cap     <- round(min(n_c) * cap_multiplier)
  # Upsampled classes are capped; downsampled classes keep their computed value
  ifelse(targets > n_c, pmin(targets, cap), targets)
}

# Regime C — Full upsampling (stress test)
# target_c = n_c^(1-alpha) * n_max^alpha
# alpha=0 → target_c = n_c
# alpha=1 → target_c = n_max for all (inflate everything to largest class)
.targets_up <- function(n_c, alpha) {
  n_max   <- max(n_c)
  targets <- round(n_c^(1 - alpha) * n_max^alpha)
  pmax(targets, n_c)    # strict upsampling: never shrink
}


# ── Main function ─────────────────────────────────────────────────────────────

#' Resample a Seurat object under a chosen imbalance regime
#'
#' @param seurat_obj      A Seurat object with cell-type labels.
#' @param alpha           Perturbation strength in [0, 1].
#'                        0 = no change (identity); 1 = maximum perturbation.
#' @param mode            Sampling regime:
#'                        "down"  — Regime A, pure downsampling (default).
#'                        "mixed" — Regime B, mild upsample + mild downsample.
#'                        "up"    — Regime C, full upsampling / duplication.
#' @param celltype_col    Metadata column holding cell-type labels
#'                        (default: "celltype").
#' @param cap_multiplier  Only used when mode="mixed". Maximum fold-increase
#'                        allowed for rare classes relative to n_min (default: 3).
#' @param seed            Integer random seed (default: 1234).
#' @param verbose         Print per-class summary (default: TRUE).
#'
#' @return A Seurat object tagged with metadata columns:
#'         resample_alpha  — alpha value used
#'         resample_mode   — mode used ("down", "mixed", or "up")
#'         original_cell   — source barcode (tracks duplicated cells)
#'
#' @details
#'   Regime A ("down") samples WITHOUT replacement only.
#'   Regimes B ("mixed") and C ("up") keep all original cells once, then
#'   fill the gap to the target by sampling WITH replacement. Duplicated
#'   cells receive suffixed barcodes (_dup1, _dup2, ...) so Seurat cell
#'   names remain unique. Reductions are NOT transferred for objects with
#'   duplicated cells since duplicates have no unique embedding.
#'
#' @examples
#'   # Regime A — remove cells from large classes
#'   obj_down  <- resample_cells(epitrace_obj, alpha = 0.7, mode = "down")
#'
#'   # Regime B — converge toward geometric mean
#'   obj_mixed <- resample_cells(epitrace_obj, alpha = 0.5, mode = "mixed")
#'
#'   # Regime C — inflate rare classes to stress-test EpiTrace
#'   obj_up    <- resample_cells(epitrace_obj, alpha = 1.0, mode = "up")
#'
#'   # Drop-in replacement for the old cap=350 block:
#'   epitrace_balanced <- resample_cells(
#'     epitrace_obj_age_estimated_multiome,
#'     alpha = 0.7,
#'     mode  = "down"
#'   )
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
  
  new_to_orig <- character(0)
  
  for (ct in classes) {
    pool <- all_bcs[celltypes == ct]
    tgt  <- targets[ct]
    
    if (tgt <= length(pool)) {
      drawn       <- sample(pool, size = tgt, replace = FALSE)
      new_to_orig <- c(new_to_orig, setNames(drawn, drawn))
    } else {
      gap         <- tgt - length(pool)
      extra       <- sample(pool, size = gap, replace = TRUE)
      dup_counter <- integer(0)
      extra_new   <- character(gap)
      for (j in seq_len(gap)) {
        src            <- extra[j]
        dup_counter[src] <- (dup_counter[src] %||% 0L) + 1L
        extra_new[j]   <- sprintf("%s_dup%d", src, dup_counter[src])
      }
      new_to_orig <- c(new_to_orig,
                       setNames(pool,  pool),
                       setNames(extra, extra_new))
    }
  }
  
  # ── Fast path: mode="down", no duplicates — subset keeps ALL assays ─────────
  needs_dup <- any(duplicated(unname(new_to_orig)))
  
  if (!needs_dup) {
    result_obj <- subset(seurat_obj, cells = unname(new_to_orig))
    result_obj$original_cell  <- unname(new_to_orig)[
      match(rownames(result_obj@meta.data), names(new_to_orig))]
    result_obj$resample_alpha <- alpha
    result_obj$resample_mode  <- mode
    if (verbose) .print_summary(alpha, mode, n_c, targets, classes, seed)
    return(result_obj)
  }
  
  # ── Slow path: mode="mixed" or "up", duplicated cells ───────────────────────
  orig_unique <- unique(unname(new_to_orig))
  base_obj    <- subset(seurat_obj, cells = orig_unique)
  
  orig_meta           <- base_obj@meta.data
  new_meta            <- orig_meta[new_to_orig, , drop = FALSE]
  rownames(new_meta)  <- names(new_to_orig)
  new_meta$original_cell  <- unname(new_to_orig)
  new_meta$resample_alpha <- alpha
  new_meta$resample_mode  <- mode
  
  assay_names <- Seurat::Assays(base_obj)
  
  # FIX: skip assays that have no counts layer (e.g. chromvar stores only
  # a 'data' layer; trying GetAssayData(..., layer="counts") on it returns
  # an empty matrix and CreateAssayObject() then errors.
  has_counts <- vapply(assay_names, function(an) {
    mtx <- tryCatch(
      Seurat::GetAssayData(base_obj, assay = an, layer = "counts"),
      error = function(e) NULL
    )
    !is.null(mtx) && nrow(mtx) > 0 && ncol(mtx) > 0
  }, logical(1))
  
  assay_names_ok <- assay_names[has_counts]
  
  if (any(!has_counts)) {
    message("  Note: skipping assays with no counts layer (will not be ",
            "transferred): ", paste(assay_names[!has_counts], collapse = ", "))
  }
  
  new_assay_list <- lapply(assay_names_ok, function(an) {
    mtx     <- Seurat::GetAssayData(base_obj, assay = an, layer = "counts")
    new_mtx <- mtx[, new_to_orig, drop = FALSE]
    colnames(new_mtx) <- names(new_to_orig)
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
