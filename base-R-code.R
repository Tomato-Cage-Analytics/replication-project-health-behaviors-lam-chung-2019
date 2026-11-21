# =============================================================================
# Model results — Meta-Analysis Studio 4.6.0, Two-Arm module
# Updated to follow the Methods: "Two-Arm Meta-Analysis with Robust Variance
# Estimation" (nested random-effects model, CR2 robust inference).
# Data: demo-dataset.csv
#
# Contents
#   Step 1  Effect sizes (Hedges' g; means/SDs first, other statistics only
#           for rows without them; direction so that + = improvement with IPSRT)
#   Step 2  Main analysis, run separately for each outcome group (domain):
#           A. Nested model: outcomes (es_id) within trials, REML, rho = 0.5,
#              CR2 + Satterthwaite robust inference, model-based results alongside
#              Q, I2 (total, between-trial, within-trial)
#           B. One inverse-variance-weighted composite effect per trial, pooled
#   Step 3  Sensitivity analyses
#           C1. rho = 0 to 0.9
#           C2. leave-one-trial-out (all outcomes of the trial removed)
#           C3. composite-per-trial pooling (Section B)
#           C4. restricted to confidence = "standard"
#
# Not covered here: Egger test / funnel plot / trim-and-fill (Section D) and
# the risk-of-bias module.
#
# Packages: metafor, clubSandwich
# Put this file and demo-dataset.csv in the same folder, then
# Session > Set Working Directory > To Source File Location.
# =============================================================================

library(metafor)
library(clubSandwich)

# ---- Settings (values from the Methods) -------------------------------------
DATA_FILE     <- "demo-dataset.csv"
RHO           <- 0.5                  # assumed correlation between outcomes of one trial
RHO_GRID      <- seq(0, 0.9, by = 0.1) # sensitivity analysis C1
PRE_POST_R    <- 0.5                  # pre-post r for change-score SDs (used if the sheet has no pre_post_r)
MIN_TRIALS    <- 3                    # an outcome group is pooled only with >= 3 independent trials
MIN_ROBUST_DF <- 4                    # robust df below this are flagged as unreliable (Tipton, 2015)
COMPOSITE_TEST <- "knha"              # test for pooling composites (Section B / C3); "z" for plain REML
# Cohort rule from the Methods: main analysis = cohort "original"; extended analysis =
# "original" + "new". If the file has no cohort called "original" (as in this demo),
# all rows are analysed together.
COHORT_SETS <- list("Original trials" = "original",
                    "Original + new"  = c("original", "new"))

dat <- read.csv(DATA_FILE, check.names = FALSE, stringsAsFactors = FALSE)

# =============================================================================
# Step 1. Effect sizes
# =============================================================================
col <- function(d, nm) if (nm %in% names(d)) d[[nm]] else rep(NA, nrow(d))
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
pnum <- function(x) num(sub("^\\s*[<>=]+\\s*", "", as.character(x)))   # "<0.001" -> 0.001
first_val <- function(...) {                 # first non-missing value across several columns
  m <- cbind(...); apply(m, 1, function(z) { z <- z[!is.na(z)]; if (length(z)) z[1] else NA_real_ })
}
Jc <- function(df) 1 - 3 / (4 * df - 1)       # Hedges' small-sample correction
# d (unsigned or signed) -> Hedges' g and its variance
from_d <- function(d, n1, n2) {
  j <- Jc(n1 + n2 - 2)
  c(yi = j * d, vi = j^2 * ((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2))))
}
# log odds ratio -> Hedges' g via Chinn (2000): d = lnOR * sqrt(3)/pi, Var = Var(lnOR) * 3/pi^2
from_logor <- function(lor, se, n1, n2) {
  j <- Jc(n1 + n2 - 2)
  c(yi = j * lor * sqrt(3) / pi, vi = j^2 * se^2 * 3 / pi^2)
}

n1 <- num(dat$n1); n2 <- num(dat$n2)
m1 <- num(dat$m1); sd1 <- num(dat$sd1); m2 <- num(dat$m2); sd2 <- num(dat$sd2)
pp_r <- first_val(num(col(dat, "pre_post_r")), rep(PRE_POST_R, nrow(dat)))

# Change scores where post means/SDs are missing: change = post - pre,
# SD(change) = sqrt(sd_pre^2 + sd_post^2 - 2 r sd_pre sd_post)
fill_change <- function(m, s, pre_m, pre_s, post_m, post_s) {
  need <- is.na(m) | is.na(s)
  ok <- need & !is.na(pre_m) & !is.na(pre_s) & !is.na(post_m) & !is.na(post_s)
  m[ok] <- (post_m - pre_m)[ok]
  s[ok] <- sqrt(pre_s^2 + post_s^2 - 2 * pp_r * pre_s * post_s)[ok]
  list(m = m, s = s, filled = ok)
}
a1 <- fill_change(m1, sd1, num(col(dat, "pre_m1_EXP")),  num(col(dat, "pre_sd1_EXP")),
                  num(col(dat, "post_m1_EXP")), num(col(dat, "post_sd1_EXP")))
a2 <- fill_change(m2, sd2, num(col(dat, "pre_m2_CONT")),  num(col(dat, "pre_sd2_CONT")),
                  num(col(dat, "post_m2_CONT")), num(col(dat, "post_sd2_CONT")))
m1 <- a1$m; sd1 <- a1$s; m2 <- a2$m; sd2 <- a2$s
dat$es_method <- NA_character_
dat$es_method[a1$filled | a2$filled] <- "change_score_from_pre_post"

# Direction: positive = improvement with IPSRT. Use es_sign; if it is empty, read
# scale_direction (means/SDs, events) or favours (statistic-only rows).
read_sign <- function(es_sign, scale_dir, favours, stat_only) {
  if (!is.na(es_sign)) return(es_sign)
  sd_txt <- tolower(trimws(as.character(scale_dir)))
  fv_txt <- tolower(trimws(as.character(favours)))
  if (!stat_only && !is.na(sd_txt)) {
    if (grepl("^(higher_better|event_good)", sd_txt)) return(1)
    if (grepl("^(lower_better|event_bad)", sd_txt))  return(-1)
  }
  if (!is.na(fv_txt)) {
    if (grepl("^(ipsrt|srt|intervention|exp)", fv_txt)) return(1)
    if (grepl("^(control|tau|cont)", fv_txt))           return(-1)
  }
  NA_real_
}

has_means <- !is.na(m1) & !is.na(sd1) & !is.na(n1) & !is.na(m2) & !is.na(sd2) & !is.na(n2)
dat$yi <- NA_real_; dat$vi <- NA_real_

# (1) means and SDs: Hedges' g
if (any(has_means)) {
  es <- escalc(measure = "SMD", m1i = m1[has_means], sd1i = sd1[has_means], n1i = n1[has_means],
               m2i = m2[has_means], sd2i = sd2[has_means], n2i = n2[has_means])
  dat$yi[has_means] <- es$yi; dat$vi[has_means] <- es$vi
  dat$es_method[has_means & is.na(dat$es_method)] <- "means_sd"
}

# (2) other statistics, only for rows without means/SDs, in the order of the Methods
stat_es <- function(i) {
  g_rep <- num(col(dat, "hedges_g"))[i];  d_rep <- num(col(dat, "cohens_d"))[i]
  chi   <- first_val(num(col(dat, "chi_sq")))[i]
  eta2  <- num(col(dat, "eta2"))[i]
  tt    <- first_val(num(col(dat, "t_stat")), num(col(dat, "t")))[i]
  ff    <- first_val(num(col(dat, "F_stat")), num(col(dat, "F")))[i]
  te    <- num(col(dat, "generic_TE"))[i];  se <- num(col(dat, "generic_seTE"))[i]
  pv    <- first_val(num(col(dat, "p_value")), pnum(col(dat, "P")))[i]
  rpb   <- num(col(dat, "r_pb"))[i]
  e1 <- num(col(dat, "e1"))[i]; e2 <- num(col(dat, "e2"))[i]
  lor <- num(col(dat, "log_or"))[i]; lse <- num(col(dat, "log_or_se"))[i]
  orr <- num(col(dat, "or"))[i]; olo <- num(col(dat, "or_lo"))[i]; ohi <- num(col(dat, "or_hi"))[i]
  op  <- num(col(dat, "or_p"))[i]
  N1 <- n1[i]; N2 <- n2[i]; okn <- is.finite(N1) && is.finite(N2)
  src <- as.character(col(dat, "es_source_file"))[i]
  # Statistic-based magnitudes are unsigned; the sign comes from es_sign (direction step).
  if (okn && !is.na(g_rep))  return(c(from_d(abs(g_rep) / Jc(N1 + N2 - 2), N1, N2), m = "reported_g"))
  if (okn && !is.na(d_rep))  return(c(from_d(abs(d_rep), N1, N2),
                                      m = ifelse(!is.na(src) && src == "partial_eta2", "partial_eta2 (via d)", "reported_d")))
  if (okn && !is.na(chi))    { phi <- sqrt(chi / (N1 + N2)); return(c(from_d(2 * phi / sqrt(1 - phi^2), N1, N2), m = "chi_square_phi")) }
  if (okn && !is.na(eta2))   return(c(from_d(2 * sqrt(eta2 / (1 - eta2)), N1, N2), m = "partial_eta2"))
  if (okn && !is.na(tt))     return(c(from_d(abs(tt) * sqrt(1 / N1 + 1 / N2), N1, N2), m = "t"))
  if (okn && !is.na(ff))     return(c(from_d(sqrt(ff) * sqrt(1 / N1 + 1 / N2), N1, N2), m = "F_1df"))
  if (!is.na(te) && !is.na(se)) return(c(yi = abs(te), vi = se^2, m = "precalculated_es_se"))
  if (okn && !is.na(pv))     { tq <- qt(1 - pv / 2, N1 + N2 - 2); return(c(from_d(tq * sqrt(1 / N1 + 1 / N2), N1, N2), m = "p_value")) }
  if (okn && !is.na(rpb))    return(c(from_d(2 * abs(rpb) / sqrt(1 - rpb^2), N1, N2), m = "r_pb"))
  if (okn && !is.na(e1) && !is.na(e2)) {                       # event counts (signed by direction later)
    a <- e1; b <- N1 - e1; c_ <- e2; d_ <- N2 - e2
    if (any(c(a, b, c_, d_) == 0)) { a <- a + .5; b <- b + .5; c_ <- c_ + .5; d_ <- d_ + .5 }
    return(c(from_logor(log(a * d_ / (b * c_)), sqrt(1/a + 1/b + 1/c_ + 1/d_), N1, N2), m = "events_chinn"))
  }
  if (okn && is.na(lor) && !is.na(orr) && !is.na(olo) && !is.na(ohi)) { lor <- log(orr); lse <- (log(ohi) - log(olo)) / (2 * qnorm(.975)) }
  if (okn && is.na(lor) && !is.na(orr) && !is.na(op)) { lor <- log(orr); lse <- abs(lor) / qnorm(1 - op / 2) }
  if (okn && !is.na(lor) && !is.na(lse)) return(c(from_logor(lor, lse, N1, N2), m = "odds_ratio_chinn"))
  c(yi = NA, vi = NA, m = NA)
}
for (i in which(!has_means)) {
  r <- stat_es(i)
  dat$yi[i] <- as.numeric(r["yi"]); dat$vi[i] <- as.numeric(r["vi"]); dat$es_method[i] <- r["m"]
}

# (3) direction
dat$stat_only <- !has_means & !(dat$es_method %in% c("events_chinn", "odds_ratio_chinn"))
dat$sgn <- mapply(read_sign, num(col(dat, "es_sign")), col(dat, "scale_direction"),
                  col(dat, "favours"), dat$stat_only)
if (anyNA(dat$sgn[!is.na(dat$yi)]))
  stop("Direction could not be read for: ", paste(dat$study[is.na(dat$sgn) & !is.na(dat$yi)], collapse = "; "),
       ". Set es_sign to 1 or -1 for these rows.")
# Statistic-based magnitudes (everything except means/SDs and event/OR rows) are unsigned:
unsigned <- !has_means & !(dat$es_method %in% c("events_chinn", "odds_ratio_chinn"))
dat$yi <- ifelse(unsigned, abs(dat$yi), dat$yi) * dat$sgn

# Rows without a usable effect are excluded (reported, not silently dropped)
bad <- is.na(dat$yi) | is.na(dat$vi) | !is.finite(dat$yi) | !is.finite(dat$vi)
if (any(bad)) { warning("Excluded (no usable effect size): ", paste(dat$study[bad], collapse = "; ")); dat <- dat[!bad, ] }

# Cluster and effect-size identifiers: trial = level 3, es_id = level 2
dat$trial <- if ("trial" %in% names(dat)) ifelse(is.na(dat$trial) | dat$trial == "", dat$study, dat$trial) else dat$study
dat$es_id <- if ("es_id" %in% names(dat)) dat$es_id else seq_len(nrow(dat))
dat <- dat[order(dat$trial, dat$es_id), ]
dat$domain <- if ("domain" %in% names(dat)) dat$domain else "All outcomes"

cat("\n==== Effect sizes (Hedges' g; positive = improvement with IPSRT) ====\n")
print(data.frame(study = dat$study, trial = dat$trial, es_id = dat$es_id, method = dat$es_method,
                 sign = dat$sgn, g = round(dat$yi, 3), se = round(sqrt(dat$vi), 3)), row.names = FALSE)
if (any(dat$es_method == "change_score_from_pre_post", na.rm = TRUE))
  cat("Note: change-score SDs used r =", PRE_POST_R, "unless the sheet gave pre_post_r.\n")
if (any(grepl("eta2", dat$es_method)))
  cat("Note: partial eta-squared effects are approximate; check them in leave-one-out.\n")

# =============================================================================
# Model functions
# =============================================================================
cohen_label <- function(g) { a <- abs(g); ifelse(a > 0.8, "large", ifelse(a >= 0.5, "medium", ifelse(a >= 0.2, "small", "negligible"))) }

# Nested model: outcomes (es_id) within trials, REML, working correlation rho between
# sampling errors of outcomes in the same trial; CR2 + Satterthwaite with trials as clusters.
fit_nested <- function(d, rho) {
  V <- vcalc(vi = vi, cluster = trial, obs = es_id, rho = rho, data = d)
  m <- rma.mv(yi, V, random = ~ 1 | trial/es_id, data = d, method = "REML")
  ct <- coef_test(m, vcov = "CR2", cluster = d$trial)
  ci <- conf_int(m, vcov = "CR2", cluster = d$trial)
  list(model = m, V = V, data = d, ct = ct, ci = ci, rho = rho)
}
tidy_fit <- function(f, label = NA) {
  m <- f$model; ct <- f$ct; ci <- f$ci
  data.frame(label = label, trials = length(unique(f$data$trial)), effects = nrow(f$data),
             g = unname(m$beta[1]),
             robust_lo = ci[["CI_L"]], robust_hi = ci[["CI_U"]], robust_df = ci[["df"]],
             robust_p = ct[["p_Satt"]],
             robust_df_flag = ifelse(ci[["df"]] < MIN_ROBUST_DF, "df < 4: unreliable", ""),
             model_lo = m$ci.lb, model_hi = m$ci.ub, model_p = m$pval,
             sigma2_trial = m$sigma2[1], sigma2_effect = m$sigma2[2], row.names = NULL)
}
safe_fit <- function(d, rho, label) {
  tryCatch(tidy_fit(fit_nested(d, rho), label),
           error = function(e) data.frame(label = label, trials = length(unique(d$trial)), effects = nrow(d),
                                          g = NA, robust_lo = NA, robust_hi = NA, robust_df = NA, robust_p = NA,
                                          robust_df_flag = paste("not estimable:", conditionMessage(e)),
                                          model_lo = NA, model_hi = NA, model_p = NA,
                                          sigma2_trial = NA, sigma2_effect = NA))
}
show <- function(x, digits = 3) { is_num <- vapply(x, is.numeric, logical(1)); x[is_num] <- lapply(x[is_num], round, digits); print(x, row.names = FALSE) }

# Section B: one inverse-variance-weighted composite per trial. Its variance uses the
# within-trial covariances implied by rho: Var = w'Vw / (sum w)^2, w = 1/vi.
composites <- function(d, rho) {
  V <- as.matrix(vcalc(vi = vi, cluster = trial, obs = es_id, rho = rho, data = d))
  do.call(rbind, lapply(unique(d$trial), function(tr) {
    ix <- which(d$trial == tr); w <- 1 / d$vi[ix]
    data.frame(trial = tr, k_effects = length(ix), yi = sum(w * d$yi[ix]) / sum(w),
               vi = drop(t(w) %*% V[ix, ix, drop = FALSE] %*% w) / sum(w)^2)
  }))
}
pool_composites <- function(d, rho) {
  cp <- composites(d, rho)
  list(composites = cp, fit = rma(yi, vi, data = cp, method = "REML", test = COMPOSITE_TEST))
}
tidy_comp <- function(pc, label) {
  f <- pc$fit
  data.frame(label = label, trials = f$k, g = unname(f$beta[1]), ci_lo = f$ci.lb, ci_hi = f$ci.ub, p = f$pval,
             tau2 = f$tau2, I2 = f$I2, row.names = NULL)
}

# =============================================================================
# Steps 2-3. Run for each cohort set and each outcome group
# =============================================================================
if (!"cohort" %in% names(dat) || !any(tolower(dat$cohort) == "original", na.rm = TRUE)) {
  cat("\nNo cohort called 'original' in this file: analysing all rows together.\n")
  cohort_sets <- list("All rows in file" = NULL)
} else cohort_sets <- COHORT_SETS

for (cs in names(cohort_sets)) {
  dc <- if (is.null(cohort_sets[[cs]])) dat else dat[tolower(dat$cohort) %in% cohort_sets[[cs]], ]
  for (dom in unique(dc$domain)) {
    d <- dc[dc$domain == dom, ]
    d <- d[order(d$trial, d$es_id), ]
    hdr <- paste0(cs, " | outcome group: ", dom)
    cat("\n\n", strrep("=", 78), "\n", hdr, "\n", strrep("=", 78), "\n", sep = "")
    n_tr <- length(unique(d$trial))
    if (n_tr < MIN_TRIALS) {
      cat("Only", n_tr, "independent trial(s); the Methods require at least", MIN_TRIALS, "-- not pooled.\n"); next
    }
    
    # ---- A. Nested model ------------------------------------------------------
    f <- fit_nested(d, RHO); m <- f$model; res <- tidy_fit(f, "nested, rho = 0.5")
    cat("\n-- A. Nested random-effects model (REML, rho =", RHO, ") --\n")
    cat(sprintf("Trials = %d, effects = %d\n", res$trials, res$effects))
    cat(sprintf("Pooled g = %.3f (%s)\n", res$g, cohen_label(res$g)))
    cat(sprintf("Robust (CR2): 95%% CI [%.3f, %.3f], df = %.2f, p = %.4f %s\n",
                res$robust_lo, res$robust_hi, res$robust_df, res$robust_p,
                ifelse(res$robust_df < MIN_ROBUST_DF, "  <-- df < 4, robust inference unreliable", "")))
    cat(sprintf("Model-based : 95%% CI [%.3f, %.3f], p = %.4f\n", res$model_lo, res$model_hi, res$model_p))
    sig_rob <- sign(res$robust_lo) == sign(res$robust_hi)      # CI excludes 0
    sig_mb  <- sign(res$model_lo) == sign(res$model_hi)
    if (!is.na(sig_rob) && !is.na(sig_mb) && sig_rob != sig_mb)
      cat("Robust and model-based results disagree on significance -- say so in the paper.\n")
    # Heterogeneity: Q, I2 total / between trials / within trials (Cheung, 2014)
    w <- 1 / d$vi; k <- nrow(d)
    s2 <- (k - 1) * sum(w) / (sum(w)^2 - sum(w^2))              # typical sampling variance
    tot <- sum(m$sigma2) + s2
    I2_b <- 100 * m$sigma2[1] / tot; I2_w <- 100 * m$sigma2[2] / tot
    cat(sprintf("Q(%d) = %.2f, p = %.4f\n", m$k - m$p, m$QE, m$QEp))
    cat(sprintf("I2 total = %.1f%%  (between trials = %.1f%%, within trials = %.1f%%)\n", I2_b + I2_w, I2_b, I2_w))
    cat(sprintf("tau2 between trials = %.4f, within trials = %.4f\n", m$sigma2[1], m$sigma2[2]))
    cat("Heterogeneity:", ifelse(m$QEp > 0.05 && (I2_b + I2_w) <= 50, "low (Q p > .05 and I2 <= 50%)", "not low"), "\n")
    
    # ---- B. One composite per trial ------------------------------------------
    pcmp <- pool_composites(d, RHO)
    cat("\n-- B. One weighted composite effect per trial --\n")
    show(pcmp$composites)
    cat("Pooled composites (REML, test =", COMPOSITE_TEST, "):\n"); show(tidy_comp(pcmp, "composite per trial"))
    
    # ---- C. Sensitivity analyses ---------------------------------------------
    cat("\n-- C1. Sensitivity to rho (0 to 0.9) --\n")
    show(do.call(rbind, lapply(RHO_GRID, function(r) safe_fit(d, r, paste0("rho = ", r)))))
    
    cat("\n-- C2. Leave-one-trial-out (all outcomes of the trial removed; rho =", RHO, ") --\n")
    show(do.call(rbind, lapply(unique(d$trial), function(tr)
      safe_fit(d[d$trial != tr, ], RHO, paste("without", tr)))))
    cat("(With few trials, leaving one out leaves very few clusters; robust df are then small or not estimable.)\n")
    
    cat("\n-- C3. Composite-per-trial pooling vs all outcomes --\n")
    show(rbind(data.frame(label = "all outcomes (nested, robust)", trials = res$trials, g = res$g,
                          ci_lo = res$robust_lo, ci_hi = res$robust_hi, p = res$robust_p, tau2 = NA, I2 = NA),
               tidy_comp(pcmp, "composite per trial")))
    
    if ("confidence" %in% names(d) && any(d$confidence == "standard", na.rm = TRUE) &&
        any(d$confidence != "standard", na.rm = TRUE)) {
      ds <- d[d$confidence == "standard", ]
      cat("\n-- C4. Restricted to confidence = 'standard' --\n")
      if (length(unique(ds$trial)) < MIN_TRIALS)
        cat("Only", length(unique(ds$trial)), "trial(s) remain after restriction; fewer than", MIN_TRIALS, "-- reported for reference only.\n")
      show(safe_fit(ds, RHO, "confidence = standard"))
    }
  }
}

sessionInfo()