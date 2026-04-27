#!/usr/bin/env Rscript
# Fit one of the multivariate hormone-imputation models (BC / FC / LH) on the
# BioCycle data. Designed to be launched on a single RStudio node so that the
# three models can run in parallel across nodes.
#
# Usage: Rscript fit_mv.R {bc|fc|lh}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1 || !args[[1]] %in% c("bc", "fc", "lh")) {
  stop("Usage: Rscript fit_mv.R {bc|fc|lh}")
}
which_model <- args[[1]]
cat(sprintf("[fit_mv] starting model: %s\n", which_model))

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(rio)
})

options(mc.cores = parallel::detectCores(),
        brms.backend = "cmdstanr",
        brms.file_refit = "on_change")

biocycle <- readRDS("biocycle.rds") %>% filter(between(cycle_length, 20, 35))
biocycle$logE   <- log(biocycle$total_estradiol)
biocycle$logFE  <- log(biocycle$estradiol)
biocycle$logP   <- log(biocycle$progesterone)
biocycle$logT   <- log(biocycle$Testo)
biocycle$logFT  <- log(biocycle$`Free T (ng/dL)`)
biocycle$logLH  <- log(biocycle$LH)
biocycle$logFSH <- log(biocycle$FSH)

if (which_model == "bc") {
  day_var <- "bc_day"
  data_used <- biocycle
  newdata <- data.frame(
    bc_day = c(-28:-1, -29:-40),
    prc_stirn_bc = c(.01, .01, .02, .03, .05, .09, .16, .27, .38, .48, .56, .58, .55, .48, .38, .28, .20, .14, .10, .07, .06, .04, .03, .02, .01, .01, .01, .01, rep(.01, times = 12)),
    prc_wcx_bc   = c(.000, .000, .001, .002, .004, .009, .018, .032, .050, .069, .085, .094, .093, .085, .073, .059, .047, .036, .028, .021, .016, .013, .010, .008, .007, .006, .005, .005, rep(.005, times = 12))
  ) %>% arrange(bc_day)
  model_file <- "models/mod_mv_bc"
  out_stub   <- "merge_files/bc_days_mv"
  suffix     <- "_bc"
} else if (which_model == "fc") {
  day_var <- "fc_day"
  data_used <- biocycle
  newdata <- data.frame(
    fc_day = c(0:39),
    prc_stirn_fc = c(.01, .01, .02, .03, .05, .09, .16, .27, .38, .48, .56, .58, .55, .48, .38, .28, .20, .14, .10, .07, .06, .04, .03, .02, .01, .01, .01, .01, rep(.01, times = 12)),
    prc_wcx_fc   = c(.000, .000, .001, .002, .004, .009, .018, .032, .050, .069, .085, .094, .093, .085, .073, .059, .047, .036, .028, .021, .016, .013, .010, .008, .007, .006, .005, .005, rep(.005, times = 12))
  )
  model_file <- "models/mod_mv_fc"
  out_stub   <- "merge_files/fc_days_mv"
  suffix     <- "_fc"
} else {
  day_var <- "lh_day"
  data_used <- biocycle %>% filter(between(lh_day, -20, 14))
  newdata <- tibble(
    conception_risk_lh = c(rep(0,8), 0.00, 0.01, 0.02, 0.06, 0.16, 0.20, 0.25, 0.24, 0.10, 0.02, 0.02, rep(0,12)),
    lh_day = -15:15
  ) %>% mutate(fertile_lh = conception_risk_lh / max(conception_risk_lh) * 0.25 / 0.31)
  model_file <- "models/mod_mv_lh"
  out_stub   <- "merge_files/lh_days_mv"
  suffix     <- "_lh"
}

mk_bf <- function(resp) {
  bf(as.formula(sprintf("%s ~ s(%s) + (1 | id) + (1|id:cycle)", resp, day_var)))
}

bf_set <- mk_bf("logE") + mk_bf("logFE") + mk_bf("logP") + mk_bf("logT") +
          mk_bf("logFT") + mk_bf("logLH") + mk_bf("logFSH") + set_rescor(TRUE)

mod <- brm(bf_set,
           data = data_used,
           control = list(adapt_delta = 0.95, max_treedepth = 20),
           file = model_file)
print(mod)

preds <- predict(mod, newdata = newdata, re_formula = NA)

newdata[[paste0("mv_est_estradiol",         suffix)]] <- round(preds[, "Estimate", "logE"],   2)
newdata[[paste0("mv_est_free_estradiol",    suffix)]] <- round(preds[, "Estimate", "logFE"],  2)
newdata[[paste0("mv_est_progesterone",      suffix)]] <- round(preds[, "Estimate", "logP"],   2)
newdata[[paste0("mv_est_testosterone",      suffix)]] <- round(preds[, "Estimate", "logT"],   2)
newdata[[paste0("mv_est_free_testosterone", suffix)]] <- round(preds[, "Estimate", "logFT"],  2)
newdata[[paste0("mv_est_luteinising",       suffix)]] <- round(preds[, "Estimate", "logLH"],  2)
newdata[[paste0("mv_est_fsh",               suffix)]] <- round(preds[, "Estimate", "logFSH"], 2)

rio::export(newdata, paste0(out_stub, ".rds"))
rio::export(newdata, paste0(out_stub, ".tsv"))
rio::export(newdata, paste0(out_stub, ".sav"))

cat(sprintf("[fit_mv] done: %s -> %s.{rds,tsv,sav}\n", which_model, out_stub))
