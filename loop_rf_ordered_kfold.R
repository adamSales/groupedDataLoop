#' Random Forest Imputation (ranger, ordered factors; GLOBAL K-fold cross-fit, per arm)
#'
#' Produces no-outcome leakage, treatment-assignment-agnostic predictions:
#' - Uses ONE global K-fold partition for ALL units.
#' - For each arm (Tr=1 and Tr=0) and each fold k, trains a model on that arm's units
#'   with fold != k; then predicts ALL units in fold k with that model.
#'   => For every unit i, both \code{that[i]} and \code{chat[i]} are computed by models
#'      that do NOT use Y_i, and selection does not depend on Tr_i.
#'
#' @param Y Numeric outcome (can be 0/1). No NAs.
#' @param Tr 0/1 treatment. No NAs.
#' @param Z Data.frame/matrix of pre-treatment covariates.
#' @param K Number of global folds (default 10).
#' @param rare_threshold Collapse factor levels with global count <= rare_threshold into "___RARE___" (default 1).
#' @param num.trees Trees per fold model (default 500).
#' @param seed RNG seed for fold assignment (default 1).
#' @return An n x 2 matrix with columns:
#'   \itemize{
#'     \item \code{that} = E[Y(1)|Z]
#'     \item \code{chat} = E[Y(0)|Z]
#'   }
#'         The returned object has an attribute \code{"forests"} with fitted fold models
#'         and the global \code{fold_id}.
#' @export

loop_rf_ordered_kfold <- function(Y, Tr, Z, K = 10, rare_threshold = 1,
                                  num.trees = 500, seed = 1, ...) {
  # ---- coerce types up front ----
  Y  <- as.numeric(Y)      # handles 1-col matrices too
  Tr <- as.integer(Tr)     # TRUE/FALSE -> 1/0, numeric stays numeric
  # -------- input checks --------
  if (anyNA(Y))  stop("Error: Missing data in dependent variable Y.")
  if (anyNA(Tr)) stop("Error: Missing data in treatment indicator Tr.")
  if (length(Y) != length(Tr)) stop("Y and Tr lengths differ.")
  Zdf <- as.data.frame(Z)
  if (nrow(Zdf) != length(Y)) stop("Row count of Z must match length of Y.")
  if (!all(Tr %in% c(0, 1))) stop("Tr must be 0/1.")
  if (!requireNamespace("ranger", quietly = TRUE)) stop("Package 'ranger' is required.")
  set.seed(seed)
  
  n <- length(Y)
  idx_t <- which(Tr == 1L)
  idx_c <- which(Tr == 0L)
  
  # -------- helper: coerce ranger predictions to numeric --------
  # Regression -> numeric; Classification (probability=TRUE) -> pick P(class-2) by name if possible, else 2nd/last col fallback.
  .rg_pred_numeric <- function(pred_obj) {
    p <- pred_obj$predictions
    if (!is.matrix(p)) return(as.numeric(p))
    cn <- colnames(p)
    if (!is.null(cn)) {
      if ("1" %in% cn)    return(as.numeric(p[, "1"]))
      if ("TRUE" %in% cn) return(as.numeric(p[, "TRUE"]))
      if (all(c("0","1") %in% cn))       return(as.numeric(p[, "1"]))
      if (all(c("FALSE","TRUE") %in% cn)) return(as.numeric(p[, "TRUE"]))
      if (ncol(p) == 2L) return(as.numeric(p[, 2]))
      return(as.numeric(p[, ncol(p)]))
    } else {
      if (ncol(p) == 2L) return(as.numeric(p[, 2]))
      return(as.numeric(p[, ncol(p)]))
    }
  }
  
  # -------- Z preprocessing: factors + rare levels + simple NA handling (outcome-free) --------
  is_cat <- vapply(Zdf, function(v) is.character(v) || is.factor(v), logical(1))
  # Categorical: cast to factor, collapse rare levels, tag missing
  if (any(is_cat)) {
    for (nm in names(Zdf)[is_cat]) {
      f <- as.factor(Zdf[[nm]])
      # mark missing explicitly
      f <- addNA(f)
      if (any(is.na(levels(f)))) {
        lv <- levels(f)
        lv[is.na(lv)] <- "___MISSING___"
        f <- factor(f, levels = lv)
      }
      # collapse rare levels by global counts (excluding NAs which are now a level)
      cnt <- table(f, useNA = "no")
      rare_lvls <- names(cnt)[cnt <= rare_threshold & names(cnt) != "___MISSING___"]
      if (length(rare_lvls)) {
        f_chr <- as.character(f)
        f_chr[f_chr %in% rare_lvls] <- "___RARE___"
        f <- factor(f_chr)
      }
      # lock: keep the collapsed/marked factor
      Zdf[[nm]] <- f
    }
  }
  # Numeric: median-impute NA (outcome-free)
  is_num <- vapply(Zdf, is.numeric, logical(1))
  if (any(is_num)) {
    for (nm in names(Zdf)[is_num]) {
      if (anyNA(Zdf[[nm]])) {
        med <- stats::median(Zdf[[nm]], na.rm = TRUE)
        if (is.na(med)) stop(sprintf("All values NA in numeric column '%s'.", nm))
        Zdf[[nm]][is.na(Zdf[[nm]])] <- med
      }
    }
  }
  
  # -------- global K-fold assignment (assignment-agnostic) --------
  if (K < 2L) stop("K must be at least 2.")
  if (K > n)  stop("K cannot exceed sample size.")
  fold_id <- sample(rep_len(seq_len(K), n))  # one global fold for ALL units
  
  # storage
  that <- rep(NA_real_, n)  # E[Y(1)|Z]
  chat <- rep(NA_real_, n)  # E[Y(0)|Z]
  forest1_fits <- vector("list", K)  # treated-arm models, one per fold (leave that fold out)
  forest0_fits <- vector("list", K)  # control-arm models, one per fold (leave that fold out)
  
  # -------- fit K treated-arm models (each leaves out treated units in fold k) --------
  if (length(idx_t)) {
    for (k in seq_len(K)) {
      tr_idx <- idx_t[ fold_id[idx_t] != k ]
      if (length(tr_idx) < 2L) {
        forest1_fits[[k]] <- NULL
        next
      }
      fit_t <- ranger::ranger(
        Y ~ ., data = data.frame(Y = Y[tr_idx], Zdf[tr_idx, , drop = FALSE]),
        respect.unordered.factors = "order",
        num.trees = num.trees,
        oob.error = FALSE
      )
      forest1_fits[[k]] <- fit_t
      # predict ALL units in fold k (treated AND control) with this treated-arm model
      ii <- which(fold_id == k)
      if (length(ii)) {
        that[ii] <- .rg_pred_numeric(predict(fit_t, data = Zdf[ii, , drop = FALSE]))
      }
    }
  }
  
  # -------- fit K control-arm models (each leaves out control units in fold k) --------
  if (length(idx_c)) {
    for (k in seq_len(K)) {
      tr_idx <- idx_c[ fold_id[idx_c] != k ]
      if (length(tr_idx) < 2L) {
        forest0_fits[[k]] <- NULL
        next
      }
      fit_c <- ranger::ranger(
        Y ~ ., data = data.frame(Y = Y[tr_idx], Zdf[tr_idx, , drop = FALSE]),
        respect.unordered.factors = "order",
        num.trees = num.trees,
        oob.error = FALSE
      )
      forest0_fits[[k]] <- fit_c
      # predict ALL units in fold k with this control-arm model
      ii <- which(fold_id == k)
      if (length(ii)) {
        chat[ii] <- .rg_pred_numeric(predict(fit_c, data = Zdf[ii, , drop = FALSE]))
      }
    }
  }
  
  # -------- sanity checks --------
  if (anyNA(that)) warning("Some E[Y(1)|Z] predictions are NA (e.g., empty treated training fold). Consider increasing K.")
  if (anyNA(chat)) warning("Some E[Y(0)|Z] predictions are NA (e.g., empty control training fold). Consider increasing K.")
  
  # output
  out <- cbind(that = that, chat = chat)
  attr(out, "forests") <- list(
    forest1_fits = forest1_fits,
    forest0_fits = forest0_fits,
    fold_id = fold_id,
    spec = list(
      K = K,
      rare_threshold = rare_threshold,
      num.trees = num.trees,
      respect.unordered.factors = "order",
      seed = seed
    )
  )
  out
}
