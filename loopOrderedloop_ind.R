#' Combined LOOP imputation: ordered ranger K-fold with ind. random forest fallback
#'
#' This function estimates arm-specific conditional means
#'   \eqn{E[Y(1) \mid Z]} and \eqn{E[Y(0) \mid Z]}
#' using:
#'
#' 1. Primary path: ranger with
#'    \code{respect.unordered.factors = "order"} and a global K-fold
#'    cross-fitting scheme that is assignment-agnostic (folds do not depend
#'    on treatment).
#'
#' 2. Per-fold, per-arm, and ind. fallback: for a given fold \eqn{k} and arm
#'    (treated & control), if a categorical level appears in that fold-arm's
#'    test set but has zero training rows for that fold-arm, we avoid relying
#'    on ranger for those specific rows (which just guess on
#'    unseen levels). For those rows only, we use predictions from a
#'    numeric-only \code{loop_rf} fit (all categorical covariates dropped
#'    from Z). Other rows in that fold-arm still use ranger.
#'
#' Thus:
#'    "Safe" rows (all categorical levels seen in that arm's training for
#'     the fold-arm) use ordered \code{ranger} K-fold predictions.
#'    "Unsafe" rows (at least one unseen categorical level for that arm/fold)
#'     inherit predictions from a single global numeric-only \code{loop_rf}
#'     model.
#'
#' The cross-fitting is always outcome-agnostic and treatment-assignment-agnostic
#' in how folds are constructed.
#'
#' @param Y Numeric outcome (can be 0/1). No NAs.
#' @param Tr 0/1 treatment indicator. No NAs.
#' @param Z Data.frame or matrix of pre-treatment covariates.
#' @param K Number of global folds (default 30).
#' @param num.trees Number of trees per ranger model (default 500).
#' @param seed RNG seed for fold assignment (default 1).
#' @param fallback_to_loop_rf Logical; if TRUE (default), rows with unseen
#'   categorical levels (per arm/fold) use numeric-only \code{loop_rf}
#'   predictions instead of ranger.
#' @param dropobs Passed through to \code{loop_rf} in the fallback case.
#' @param verbose_warnings Logical; if FALSE (default), print warnings about
#'   unseen levels, fallback usage, and NA predictions. If FALSE, these are
#'   suppressed.
#' @param ... Additional arguments passed to \code{ranger::ranger()} for the
#'   K-fold ranger fits.
#'
#' @return An \eqn{n \times 2} numeric matrix with columns:
#'   \itemize{
#'     \item \code{that} = \eqn{E[Y(1)\mid Z]}
#'     \item \code{chat} = \eqn{E[Y(0)\mid Z]}
#'   }
#'
#' In the ranger-K-fold path, the return value has an attribute \code{"forests"}
#' containing:
#'   \itemize{
#'     \item \code{forest1_fits}: length-K list of treated-arm ranger models
#'           (some entries may be \code{NULL} when a fold has too few training rows),
#'     \item \code{forest0_fits}: analogous list for the control arm,
#'     \item \code{fold_id}: length-n integer vector of global fold assignments,
#'     \item \code{spec}: list of main tuning hyperparameters.
#'   }
#'
#' If \code{fallback_to_loop_rf = FALSE}, ranger is used whenever possible, and
#' no \code{loop_rf} call is made.
#'
#' @export
loopOrderedloop <- function(
    Y, Tr, Z,
    K = 30,
    num.trees = 500,
    seed = 1,
    fallback_to_loop_rf = TRUE,
    dropobs = NULL,
    verbose_warnings = FALSE,
    ...
) {
  # ---- basic type coercions ----
  Y  <- as.numeric(Y)      # handles 1-col matrices
  Tr <- as.integer(Tr)     # TRUE/FALSE -> 1/0, numeric -> integer
  
  # ---- input checks ----
  if (anyNA(Y))  stop("Error: Missing data in dependent variable Y.")
  if (anyNA(Tr)) stop("Error: Missing data in treatment indicator Tr.")
  if (length(Y) != length(Tr)) stop("Y and Tr lengths differ.")
  
  Zdf <- as.data.frame(Z)
  if (nrow(Zdf) != length(Y)) stop("Row count of Z must match length of Y.")
  if (!all(Tr %in% c(0L, 1L))) stop("Tr must be 0/1.")
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package 'ranger' is required.")
  }
  
  set.seed(seed)
  
  n     <- length(Y)
  idx_t <- which(Tr == 1L)
  idx_c <- which(Tr == 0L)
  
  # ---- helper: coerce ranger predictions to numeric ----
  # Regression -> numeric vector.
  # Classification (probability=TRUE) -> choose column corresponding to
  # "1" or "TRUE" if possible, else 2nd or last column as fallback.
  .rg_pred_numeric <- function(pred_obj) {
    p <- pred_obj$predictions
    if (!is.matrix(p)) return(as.numeric(p))
    
    cn <- colnames(p)
    if (!is.null(cn)) {
      idx <- NULL
      if ("1" %in% cn) {
        idx <- match("1", cn)
      } else if ("TRUE" %in% cn) {
        idx <- match("TRUE", cn)
      } else if (all(c("0", "1") %in% cn)) {
        idx <- match("1", cn)
      } else if (all(c("FALSE", "TRUE") %in% cn)) {
        idx <- match("TRUE", cn)
      } else if (ncol(p) == 2L) {
        idx <- 2L
      } else {
        idx <- ncol(p)
      }
      return(as.numeric(p[, idx]))
    }
    
    # no colnames
    if (ncol(p) == 2L) return(as.numeric(p[, 2L]))
    as.numeric(p[, ncol(p)])
  }
  
  # ---- Z preprocessing: factors + outcome-free NA imputation ----
  is_cat <- vapply(Zdf, function(v) is.character(v) || is.factor(v), logical(1L))
  
  # Categorical: cast to factor, tag missing explicitly as "___MISSING___"
  if (any(is_cat)) {
    for (nm in names(Zdf)[is_cat]) {
      f <- as.factor(Zdf[[nm]])
      f <- base::addNA(f)
      if (any(is.na(levels(f)))) {
        lv <- levels(f)
        lv[is.na(lv)] <- "___MISSING___"
        f <- factor(f, levels = lv)
      }
      Zdf[[nm]] <- f
    }
  }
  
  # Numeric: median-impute NA (outcome-free)
  is_num <- vapply(Zdf, is.numeric, logical(1L))
  if (any(is_num)) {
    for (nm in names(Zdf)[is_num]) {
      if (anyNA(Zdf[[nm]])) {
        med <- stats::median(Zdf[[nm]], na.rm = TRUE)
        if (is.na(med)) {
          stop(sprintf("All values NA in numeric column '%s'.", nm))
        }
        Zdf[[nm]][is.na(Zdf[[nm]])] <- med
      }
    }
  }
  
  # ---- global K-fold assignment (assignment-agnostic) ----
  if (K < 2L) stop("K must be at least 2.")
  if (K > n)  stop("K cannot exceed sample size.")
  fold_id <- sample(rep_len(seq_len(K), n))
  
  # ---- per-arm, per-fold row-level zero-train-level checks ----
  bad_fold_treated <- rep(FALSE, K)
  bad_fold_control <- rep(FALSE, K)
  bad_rows_treated <- vector("list", K)  # each element: integer indices of unsafe rows
  bad_rows_control <- vector("list", K)
  
  if (any(is_cat)) {
    cat_names <- names(Zdf)[is_cat]
    
    for (arm in c(0L, 1L)) {
      idx_arm <- which(Tr == arm)
      if (!length(idx_arm)) next
      
      for (k in seq_len(K)) {
        train_idx <- idx_arm[fold_id[idx_arm] != k]
        test_idx  <- which(fold_id == k)
        if (!length(train_idx) || !length(test_idx)) next
        
        all_missing_lvls <- character(0)
        unsafe_rows_k    <- integer(0)
        
        for (nm in cat_names) {
          f <- Zdf[[nm]]
          
          test_lvls  <- unique(f[test_idx])
          train_lvls <- unique(f[train_idx])
          
          missing_lvls <- setdiff(test_lvls, train_lvls)
          missing_lvls <- missing_lvls[!is.na(missing_lvls)]
          
          if (length(missing_lvls)) {
            all_missing_lvls <- union(all_missing_lvls, as.character(missing_lvls))
            new_bad_rows <- test_idx[f[test_idx] %in% missing_lvls]
            unsafe_rows_k <- union(unsafe_rows_k, new_bad_rows)
          }
        }
        
        if (length(unsafe_rows_k)) {
          if (arm == 1L) {
            bad_fold_treated[k] <- TRUE
            bad_rows_treated[[k]] <- unsafe_rows_k
          } else {
            bad_fold_control[k] <- TRUE
            bad_rows_control[[k]] <- unsafe_rows_k
          }
          if (verbose_warnings) {
            warning(sprintf(
              "Detected levels in categorical covariates that appear in fold %d (arm=%d) but have zero training rows for that arm/fold: %s.",
              k, arm, paste(all_missing_lvls, collapse = ", ")
            ))
          }
        }
      }
    }
  }
  
  # ---- prepare numeric-only Z and (optional) random forest predictions ----
  loop_that <- loop_chat <- NULL
  if (fallback_to_loop_rf && (any(bad_fold_treated) || any(bad_fold_control))) {
    if (!exists("loop_rf", mode = "function")) {
      stop("fallback_to_loop_rf = TRUE but loop_rf() is not available in the current environment.")
    }
    
    Z_num <- Zdf[, !is_cat, drop = FALSE]
    if (ncol(Z_num) == 0L) {
      if (verbose_warnings) {
        warning("All covariates are categorical; in fallback path, creating an intercept-only numeric covariate.")
      }
      Z_num <- matrix(1, nrow = nrow(Zdf), ncol = 1)
      colnames(Z_num) <- "intercept_only"
    }
    
    if (verbose_warnings) {
      warning("Using loop_rf with numeric-only covariates for rows where categorical levels are unseen in training.")
    }
    loop_out <- loop_rf(Y = Y, Tr = Tr, Z = Z_num, dropobs = dropobs)
    loop_that <- as.numeric(loop_out[, "that"])
    loop_chat <- as.numeric(loop_out[, "chat"])
  }
  
  # ---- storage for K-fold predictions ----
  that <- rep(NA_real_, n)  # E[Y(1)|Z]
  chat <- rep(NA_real_, n)  # E[Y(0)|Z]
  forest1_fits <- vector("list", K)  # treated-arm ranger models
  forest0_fits <- vector("list", K)  # control-arm ranger models
  
  # helper: fit one arm with row-level fallback
  .fit_arm_models <- function(idx_arm, bad_rows_list, arm_label) {
    pred <- rep(NA_real_, n)
    fits <- vector("list", K)
    
    if (!length(idx_arm)) {
      return(list(pred = pred, fits = fits))
    }
    
    for (k in seq_len(K)) {
      test_idx_k  <- which(fold_id == k)
      if (!length(test_idx_k)) next
      
      train_idx_k <- idx_arm[fold_id[idx_arm] != k]
      
      # No training data for this arm in this fold => we can't fit ranger.
      if (!length(train_idx_k)) {
        fits[[k]] <- NULL
        next
      }
      
      # Always fit ranger on training rows for this arm/fold
      fit <- ranger::ranger(
        Y ~ .,
        data   = data.frame(Y = Y[train_idx_k], Zdf[train_idx_k, , drop = FALSE]),
        respect.unordered.factors = "order",
        num.trees = num.trees,
        oob.error = FALSE,
        ...
      )
      fits[[k]] <- fit
      
      # Base predictions for all rows in this fold
      fold_pred <- .rg_pred_numeric(
        predict(fit, data = Zdf[test_idx_k, , drop = FALSE])
      )
      pred[test_idx_k] <- fold_pred
      
      # Row-level override for "unsafe" rows if fallback is enabled
      unsafe_k <- NULL
      if (!is.null(bad_rows_list) && length(bad_rows_list) >= k) {
        unsafe_k <- bad_rows_list[[k]]
      }
      if (!is.null(unsafe_k) && length(unsafe_k) &&
          fallback_to_loop_rf && !is.null(loop_that)) {
        if (arm_label == 1L) {
          pred[unsafe_k] <- loop_that[unsafe_k]
        } else {
          pred[unsafe_k] <- loop_chat[unsafe_k]
        }
      }
    }
    
    list(pred = pred, fits = fits)
  }
  
  # ---- fit treated and control arms ----
  treated_res <- .fit_arm_models(idx_t, bad_rows_treated, arm_label = 1L)
  control_res <- .fit_arm_models(idx_c, bad_rows_control, arm_label = 0L)
  
  that <- treated_res$pred  # E[Y(1)|Z]
  chat <- control_res$pred  # E[Y(0)|Z]
  
  forest1_fits <- treated_res$fits
  forest0_fits <- control_res$fits
  
  # ---- sanity checks ----
  if (anyNA(that) && verbose_warnings) {
    warning("Some E[Y(1)|Z] predictions are NA (e.g., empty treated training fold). Consider increasing sample size or adjusting K.")
  }
  if (anyNA(chat) && verbose_warnings) {
    warning("Some E[Y(0)|Z] predictions are NA (e.g., empty control training fold). Consider increasing sample size or adjusting K.")
  }
  
  # ---- output ----
  out <- cbind(that = that, chat = chat)
  attr(out, "forests") <- list(
    forest1_fits = forest1_fits,
    forest0_fits = forest0_fits,
    fold_id      = fold_id,
    spec         = list(
      K         = K,
      num.trees = num.trees,
      respect.unordered.factors = "order",
      seed      = seed,
      fallback_to_loop_rf = fallback_to_loop_rf,
      verbose_warnings    = verbose_warnings
    )
  )
  out
}