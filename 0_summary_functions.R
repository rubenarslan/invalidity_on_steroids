library(tidyverse)
library(brms)
library(gt)

options(mc.cores = parallel::detectCores(), 
        brms.backend = "cmdstanr",
        brms.file_refit = "on_change")

n_nonmissing <- function(x) { sum(!is.na(x)) }

rmse <- function(error) { sqrt(mean(error^2, na.rm = TRUE)) }
rmse_brms <- function(model, re_formula = NA) {
  resid <- residuals(model, summary = F, re_formula = re_formula, 
                     method = "posterior_predict")
  rmse_samples <- apply(resid, 1, rmse)
  rmse_ci <- rstantools::posterior_interval(as.matrix(rmse_samples), prob = 0.95)
  sprintf("%.2f [%.2f;%.2f]", mean(rmse_samples), rmse_ci[1], rmse_ci[2]) 
}

compute_factor <- function(day, imp_horm, real_horm, threshold) {
  median_imp <- median(exp(imp_horm[which(between(day, threshold[1], threshold[2]))]), na.rm = T)
  median_real <- median(real_horm[which(between(day, threshold[1], threshold[2]))], na.rm = T)
  factor <- median_imp/median_real
}

cycle_phase_plot <- function(summary, kind, custom_limits = NULL, size = 25) {
  model <- summary[[kind]]
  LOD <- summary$`Limit of detection`
  LOQ <- summary$LOQ
  Hormone_label <- summary$Hormone
  n_days <- nrow(model$data)
  n_women <- n_distinct(model$data$id)
  
  
  df <- model$data
  point_alpha <- 1 / nrow(df) * 150
  point_alpha <- if_else(point_alpha > 1, 1, point_alpha)
  
  if(str_detect(Hormone_label, "Progesterone")) {
    Hormone_label <- "Progesterone"
    color <- '#587f6f05'
    labels <- c(0.2, 1, 5, 10, 20, 50, 100, 200, 500, 1000, 1500, 10000, 25000)
    breaks <- log(labels)
    limits <- log(c(0.09, 28000))
  } else if(str_detect(Hormone_label, "Estradiol")) {
    color <- '#b85a5305'
    labels <- c(0.1, 0.5, 1, 5, 10, 20, 30)
    breaks <- log(labels)
    limits <- log(c(0.09, 32))
  } else if(str_detect(Hormone_label, "Testosterone")) {
    color <- '#b85a5305'
    labels <- c(0.1, 0.5, 1, 5, 10, 20, 30, 100, 150, 200, 250)
    breaks <- log(labels)
    limits <- log(c(0.09, 260))
  }
  
  if(!is.null(custom_limits)) {
    limits <- custom_limits
  }
  if(kind == "lh_day_model") {
    x_axis <- scale_x_continuous("Day relative to LH-surge", 
                                 breaks = seq(-15, 15, by = 3))
    coords <- coord_cartesian(xlim = c(-15, 15), ylim = limits)
  } else if(kind == "bc_day_model") {
    x_axis <- scale_x_continuous("Days until next menstrual onset", 
                                 breaks = seq(-30, 0, by = 3))
    coords <- coord_cartesian(xlim = c(-30, 0), ylim = limits) 
  } else  if(kind == "fc_day_model") {
    x_axis <- scale_x_continuous("Days since last menstrual onset", 
                                 breaks = seq(0, 30, by = 3))
    coords <-coord_cartesian(xlim = c(0, 30), ylim = limits)
  }
  
  
  plot(conditional_effects(model, 
                           conditions = tibble(hormone_cens = c("none", "left")), 
                           spaghetti = TRUE, ndraws = 200), 
       points = TRUE, 
       point_args = list(alpha = point_alpha),
       spaghetti_args = list(alpha = 0.01, colour = color), 
       line_args = list(size = 0), 
       plot = FALSE)[[1]] +
    facet_null() +
    scale_y_continuous(paste0("Log(",Hormone_label,")"),
                       breaks = breaks, labels = labels) +
    x_axis +
    coords +
    geom_hline(yintercept = log(as.numeric(LOQ)),
               linetype = 'dashed') +
    geom_hline(yintercept = log(as.numeric(LOD)),
               linetype = 'solid') +
    theme_minimal(base_size = size) +
    ggtitle(summary$Dataset,subtitle = paste0("n = ", n_women, ", days = ",n_days))
}

cycle_phase_plot_ep <- function(model, Dataset, custom_limits = NULL, size = 20) {
  Hormone_label <- "Estradiol"
  kind = "ep_day_model"
  n_days <- nrow(model$data)
  n_women <- n_distinct(model$data$id)
  
  
  df <- model$data
  point_alpha <- 1 / nrow(df) * 150
  point_alpha <- if_else(point_alpha > 1, 1, point_alpha)
  
  if(str_detect(Hormone_label, "Progesterone")) {
    Hormone_label <- "Progesterone"
    color <- '#587f6f05'
    labels <- c(0.2, 1, 5, 10, 20, 50, 100, 200, 500, 1000, 1500, 10000, 25000)
    breaks <- log(labels)
    limits <- log(c(0.09, 28000))
  } else {
    color <- '#b85a5305'
    labels <- c(0.1, 0.5, 1, 5, 10, 20, 30)
    breaks <- log(labels)
    limits <- log(c(0.09, 32))
  }
  
  if(!is.null(custom_limits)) {
    limits <- custom_limits
  }
  if(kind == "ep_day_model") {
    x_axis <- scale_x_continuous("Day relative to estradiol peak", 
                                 breaks = seq(-15, 15, by = 3))
    coords <- coord_cartesian(xlim = c(-15, 15), ylim = limits)
  } else if(kind == "bc_day_model") {
    x_axis <- scale_x_continuous("Days until next menstrual onset", 
                                 breaks = seq(-30, 0, by = 3))
    coords <- coord_cartesian(xlim = c(-30, 0), ylim = limits) 
  } else  if(kind == "fc_day_model") {
    x_axis <- scale_x_continuous("Days since last menstrual onset", 
                                 breaks = seq(0, 30, by = 3))
    coords <-coord_cartesian(xlim = c(0, 30), ylim = limits)
  }
  
  r2 <- sqrt(loo_R2(model, re_formula = NA))

  plot(conditional_effects(model, 
                           conditions = tibble(hormone_cens = c("none", "left")), 
                           spaghetti = TRUE, ndraws = 200), 
       points = TRUE, 
       point_args = list(alpha = point_alpha),
       spaghetti_args = list(alpha = 0.01, colour = color), 
       line_args = list(size = 0), 
       plot = FALSE)[[1]] +
    facet_null() +
    scale_y_continuous(paste0("Log(",Hormone_label,")"),
                       breaks = breaks, labels = labels) +
    x_axis +
    coords +
    theme_minimal(base_size = size) +
    ggtitle(Dataset,subtitle = paste0("n = ", n_women, ", days = ",n_days, ", LOO-R = ", 
                                      sprintf("%.2f [%.2f;%.2f]", r2[,"Estimate"], r2[,"Q2.5"], r2[,"Q97.5"])))
}


summarise_hormone <- function(df, 
                              Dataset, 
                              Hormone, 
                              Method = "Salivary Immunoassay",
                              Citation = "",
                              LOD = NA,
                              LOQ = NA,
                              CV_intra = NA,
                              CV_inter = NA,
                              Procedure = "",
                              Scheduling = "",
                              LH_test = "") {
  # df
  # df$id
  # df$cycle
  # df$bc_day
  # df$fc_day
  # df$cycle_length
  # df$hormone
  # df$hormone_cens
  
  results <- list()
  results$Dataset <- Dataset
  results$Citation <- Citation
  results$Hormone <- Hormone
  results$Method <- Method
  results$`Limit of detection` <- LOD
  results$LOQ <- LOQ
  results$`Intraassay CV` = CV_intra
  results$`Interassay CV` = CV_inter
  results$Procedure <- Procedure
  results$Scheduling <- Scheduling
  results$LH_test <- LH_test
  results$distribution <- ggplot(df, aes(hormone, fill = hormone_cens)) + geom_histogram() + scale_x_log10() + xlab(Hormone)
  
  # Summary statistics
  results$mean <- mean(df$hormone, na.rm = T)
  results$logmean <- mean(log(df$hormone), na.rm = T)
  results$logsd <- sd(log(df$hormone), na.rm = T)
  results$median <- median(df$hormone, na.rm = T)
  results$sd <- sd(df$hormone, na.rm = T)
  results$mad <- mad(df$hormone, na.rm = T)
  rng <- range(df$hormone, na.rm = T)
  results$range <- sprintf("%.2f, %.2f", rng[1], rng[2])
  results$missing <- sum(is.na(df$hormone))
  df$outlier <- log(df$hormone) > results$logmean + results$logsd * 3
  results$outliers <- sum(df$outlier, na.rm = T)
  results$censored <- sum(df$hormone_cens != "none", na.rm = TRUE)
  
  results$n_women <- n_distinct(df$id)
  results$n_cycles <- n_distinct(paste(df$id, df$cycle))
  results$n_days <- nrow(df)
  
  minna <- function(x) {
    if(all(is.na(x))) {
      NA_real_
    } else {
      min(x, na.rm = T)
    }
  }
  
  dfsumm <- df %>% select(id, cycle, cycle_length) %>% group_by(id, cycle) %>% summarise(cycle_length = minna(cycle_length)) %>% ungroup()
  
  if(exists("age", df)) {
    dfage <- df %>% select(id, age) %>% group_by(id) %>% summarise(age = minna(age))
    
    results$age <- sprintf("%.1f±%.2f", mean(dfage$age, na.rm = T), sd(dfage$age, na.rm = T))
  } else {
    results$age <- NA
  }
  if(exists("partner", df)) {
    dfage <- df %>% select(id, partner, age) %>% arrange(age) %>% group_by(id) %>% summarise(partner = first(as.numeric(partner)))
    
    results$in_relationship <- sprintf("%.0f%%", 100*mean(dfage$partner, na.rm = T))
  } else {
    results$in_relationship <- NA
  }
  
  if(min(df$fc_day, na.rm = T) < 0) {
    stop("FC day min wrong", paste(range(df$fc_day, na.rm = T),collapse = ";"))
  }
  if(any(!is.na(df$bc_day)) && max(df$bc_day, na.rm = T) > -1) {
    stop("BC day min wrong ", paste(range(df$bc_day, na.rm = T),collapse = ";"))
  }
  
  df <- df %>% filter(!is.na(hormone), 
                      is.na(cycle_length) | between(cycle_length, 20,35))
  
  results$cycle_length <- sprintf("%.1f±%.2f", mean(dfsumm$cycle_length, na.rm = T), sd(dfsumm$cycle_length, na.rm = T))
  
  stopifnot(all(df$fc_day < 35, na.rm = TRUE))
  stopifnot(all(df$bc_day > -36, na.rm = TRUE))

  df <- df %>% group_by(id, cycle) %>% 
    mutate(hormone_diff = hormone - mean(hormone),
           log_hormone_diff = log(hormone) - mean(log(hormone))) %>% 
    ungroup()
  results$usable_n <- nrow(df)
  results$usable_n_women <- n_distinct(df$id)
  
  bc_days <- readRDS("merge_files/bc_days.rds")
  fc_days <- readRDS("merge_files/fc_days.rds")
  
  if(any(!is.na(df$bc_day)) && max(df$bc_day, na.rm = T) != -1) {
    warning("Check BC day range ", paste(range(df$bc_day, na.rm = T),collapse = ";"))
  }
  df = left_join(df, bc_days)
  
  if(min(df$fc_day, na.rm = T) != 0) {
    warning("Check FC day range ", paste(range(df$fc_day, na.rm = T),collapse = ";"))
  }
  df = left_join(df, fc_days, by = "fc_day")
  
  
  # from Jünger/Stern et al. 2018 Supplementary Material
  # Day relative to ovulation	Schwartz et al., (1980)	Wilcox et al., (1998)	Colombo & Masarotto (2000)	Weighted average
  lh_days <- readRDS("merge_files/lh_days.rds")
  
  if(is.null(df$lh_day) || all(is.na(df$lh_day))) { 
    df$lh_day <- NA_real_
    results$no_lh_surge_woman <- "NA"
    results$no_lh_surge_cycle <- "NA"
  } else {
    results$no_lh_surge_woman <- df %>% group_by(id) %>% 
      filter(!is.na(hormone), !is.na(fc_day)) %>% 
      summarise(no_lh_surge = all(is.na(lh_day))) %>% 
      summarise(count = sum(no_lh_surge), percent = mean(no_lh_surge),
                cycles = n()) %>% 
      { sprintf("%.0f/%.0f (%.0f%%)", .$count, .$cycles, .$percent*100) }
    results$no_lh_surge_cycle <- df %>% group_by(id, cycle) %>% 
      filter(!is.na(hormone), !is.na(fc_day)) %>% 
      summarise(no_lh_surge = all(is.na(lh_day))) %>% 
      ungroup() %>% 
      summarise(count = sum(no_lh_surge), percent = mean(no_lh_surge),
                cycles = n()) %>% 
      { sprintf("%.0f/%.0f (%.0f%%)", .$count, .$cycles, .$percent*100) }
  }
    
  df = df %>% left_join(lh_days, by = "lh_day")
  
  several_cycles <- { df %>% group_by(id) %>% filter(n_distinct(cycle) > 1) %>% nrow() } > 0
  several_timepoints <- { df %>% group_by(id) %>% filter(n() > 1) %>% nrow() } > 0
  
  if (several_cycles) {
    (icc_model_id <- brm(log(hormone) | cens(hormone_cens) ~ (1|id), data = df, file_refit = "on_change",
                      control = list(adapt_delta = 0.99),
                      file = paste0("models/m_", Dataset, "_", Hormone, "_icc_id"))) %>% 
      add_criterion("loo_R2")

    (icc_model_id_cycle <- brm(log(hormone) | cens(hormone_cens) ~ (1|id) + (1|id:cycle), data = df, file_refit = "on_change",
                         control = list(adapt_delta = 0.99),
                         file = paste0("models/m_", Dataset, "_", Hormone, "_icc"))) %>% 
      add_criterion("loo_R2")

    r2 <- loo_R2(icc_model_id)
    results$var_id_loo <- sprintf("%.2f [%.2f;%.2f]", r2[,"Estimate"], r2[,"Q2.5"], r2[,"Q97.5"])
    r2 <- loo_R2(icc_model_id_cycle)
    results$var_cycle_loo <- sprintf("%.2f [%.2f;%.2f]", r2[,"Estimate"], r2[,"Q2.5"], r2[,"Q97.5"])
    
    vc <- VarCorr(icc_model_id_cycle)
    r2_id <- bayes_R2(icc_model_id_cycle, re_formula = ~ (1|id))
    r2_cycle <- bayes_R2(icc_model_id_cycle, re_formula = ~ (1|id:cycle))
    
    results$var_id <- sprintf("%.2f [%.2f;%.2f] (%.0f%%)", vc$id$sd[,"Estimate"], vc$id$sd[,"Q2.5"], vc$id$sd[,"Q97.5"], 100* r2_id[,"Estimate"]) 
    results$var_cycle <- sprintf("%.2f [%.2f;%.2f] (%.0f%%)", vc$`id:cycle`$sd[,"Estimate"], vc$`id:cycle`$sd[,"Q2.5"], vc$`id:cycle`$sd[,"Q97.5"], 100* r2_cycle[,"Estimate"])
    results$var_resid <- sprintf("%.2f [%.2f;%.2f] (%.0f%%)", vc$residual$sd[,"Estimate"], vc$residual$sd[,"Q2.5"], vc$residual$sd[,"Q97.5"], (1 - r2_cycle[,"Estimate"] - r2_id[,"Estimate"]) * 100)
    
    
    if(n_nonmissing(df$bc_day) > 20) {
      (bc_day_model <- brm(log(hormone) | cens(hormone_cens) ~ s(bc_day) + (1|id) + (1|id:cycle) , data = df , file_refit = "on_change",
                           control = list(adapt_delta = 0.99),
                           file = paste0("models/m_", Dataset, "_", Hormone, "_bc")) %>% 
         add_criterion("loo_R2", re_formula = NA) %>%
         add_criterion("bayes_R2", re_formula = NA))
    }
    
    if(n_nonmissing(df$fc_day) > 20) {
      (fc_day_model <- brm(log(hormone) | cens(hormone_cens) ~ s(fc_day) + (1|id) + (1|id:cycle) , data = df , file_refit = "on_change",
                           control = list(adapt_delta = 0.99),
                           file = paste0("models/m_", Dataset, "_", Hormone, "_fc")) %>% 
         add_criterion("loo_R2", re_formula = NA) %>%
         add_criterion("bayes_R2", re_formula = NA))
    }
    
    if(n_nonmissing(df$lh_day)>20) {
      (lh_day_model <- brm(log(hormone) | cens(hormone_cens) ~ s(lh_day) + (1|id) + (1|id:cycle) , data = df %>% filter(between(lh_day, -15, 15)) ,
                           control = list(adapt_delta = 0.99),
                           file_refit = "on_change", file = paste0("models/m_", Dataset, "_", Hormone, "_lh")) %>% 
         add_criterion("loo_R2", re_formula = NA) %>%
         add_criterion("bayes_R2", re_formula = NA))
    }
    
  } else if(several_timepoints) 
  {
    (icc_model <- brm(log(hormone) | cens(hormone_cens) ~ (1|id), data = df, file_refit = "on_change",
                      control = list(adapt_delta = 0.99),
                      file = paste0("models/m_", Dataset, "_", Hormone, "_icc")))

    r2 <- loo_R2(icc_model)
    results$var_id_loo <- sprintf("%.2f [%.2f;%.2f]", r2[,"Estimate"], r2[,"Q2.5"], r2[,"Q97.5"])
    
    vc <- VarCorr(icc_model)
    r2_id <- bayes_R2(icc_model, re_formula = ~ (1|id))

    results$var_id <- sprintf("%.2f [%.2f;%.2f] (%.0f%%)", vc$id$sd[,"Estimate"], vc$id$sd[,"Q2.5"], vc$id$sd[,"Q97.5"], 100* r2_id[,"Estimate"]) 
    results$var_resid <- sprintf("%.2f [%.2f;%.2f] (%.0f%%)", vc$residual$sd[,"Estimate"], vc$residual$sd[,"Q2.5"], vc$residual$sd[,"Q97.5"], (1 - r2_id[,"Estimate"]) * 100)
    
    
    if(n_nonmissing(df$bc_day)>20) {
      (bc_day_model <- brm(log(hormone) | cens(hormone_cens) ~ s(bc_day) + (1|id), data = df , file_refit = "on_change",
                           control = list(adapt_delta = 0.99),
                           file = paste0("models/m_", Dataset, "_", Hormone, "_bc")) %>% 
         add_criterion("loo_R2", re_formula = NA) %>%
         add_criterion("bayes_R2", re_formula = NA))
    }
    
    if(n_nonmissing(df$fc_day)>20) {
      (fc_day_model <- brm(log(hormone) | cens(hormone_cens) ~ s(fc_day) + (1|id), data = df , file_refit = "on_change",
                           control = list(adapt_delta = 0.99),
                           file = paste0("models/m_", Dataset, "_", Hormone, "_fc")) %>% 
         add_criterion("loo_R2", re_formula = NA) %>%
         add_criterion("bayes_R2", re_formula = NA))
    }
    
    if(n_nonmissing(df$lh_day)>20) {
      (lh_day_model <- brm(log(hormone) | cens(hormone_cens) ~ s(lh_day) + (1|id) , data = df %>% filter(between(lh_day, -15, 15)) ,
                           control = list(adapt_delta = 0.99),
                           file_refit = "on_change", file = paste0("models/m_", Dataset, "_", Hormone, "_lh")) %>% 
         add_criterion("loo_R2", re_formula = NA) %>%
         add_criterion("bayes_R2", re_formula = NA))
    }
  } else 
  {
    
    if(n_nonmissing(df$bc_day)>20) {
      (bc_day_model <- brm(log(hormone) | cens(hormone_cens) ~ s(bc_day), data = df , file_refit = "on_change",
                           control = list(adapt_delta = 0.99),
                           file = paste0("models/m_", Dataset, "_", Hormone, "_bc")) %>% 
         add_criterion("loo_R2", re_formula = NA) %>%
         add_criterion("bayes_R2", re_formula = NA))
    }
    
    if(n_nonmissing(df$fc_day)>20) {
      (fc_day_model <- brm(log(hormone) | cens(hormone_cens) ~ s(fc_day), data = df , file_refit = "on_change",
                           control = list(adapt_delta = 0.99),
                           file = paste0("models/m_", Dataset, "_", Hormone, "_fc")) %>% 
         add_criterion("loo_R2", re_formula = NA) %>%
         add_criterion("bayes_R2", re_formula = NA))
    }
  }
  
  Hormone_name <- case_when(
    str_detect(Hormone, "Estradiol") ~ "free_estradiol", 
    str_detect(Hormone, "Progesterone") ~ "progesterone", 
    str_detect(Hormone, "Testosterone") ~ "free_testosterone",
    TRUE ~ "free_estradiol")
  
  if(n_nonmissing(df$bc_day)>20) {
    results$bc_day_model <- bc_day_model
    r2 <- bayes_R2(bc_day_model, re_formula = NA)
    results$r_bc <- sprintf("%.2f [%.2f;%.2f]", sqrt(r2[,"Estimate"]), sqrt(r2[,"Q2.5"]), sqrt(r2[,"Q97.5"])) 
    
    results$rmse_bc <- rmse_brms(bc_day_model)
    
    r2 <- loo_R2(bc_day_model, re_formula = NA)
    results$r_loo_bc <- sprintf("%.2f [%.2f;%.2f]", sqrt(r2[,"Estimate"]), sqrt(r2[,"Q2.5"]), sqrt(r2[,"Q97.5"])) 
    
    cor_stirn <- broom::tidy(cor.test(df$prc_stirn_bc, df$hormone))
    results$r_bc_stirn <- sprintf("% .2f [% .2f;% .2f]", cor_stirn["estimate"], cor_stirn["conf.low"], cor_stirn["conf.high"])

    cor_stirn <- broom::tidy(cor.test(df$prc_stirn_bc, log(df$hormone)))
    results$r_log_bc_stirn <- sprintf("% .2f [% .2f;% .2f]", cor_stirn["estimate"], cor_stirn["conf.low"], cor_stirn["conf.high"])

    cor_stirn <- broom::tidy(cor.test(df$prc_stirn_bc, df$hormone_diff))
    results$r_diff_bc_stirn <- sprintf("% .2f [% .2f;% .2f]", cor_stirn["estimate"], cor_stirn["conf.low"], cor_stirn["conf.high"])

    cor_stirn <- broom::tidy(cor.test(df$prc_stirn_bc, df$log_hormone_diff))
    results$r_log_diff_bc_stirn <- sprintf("% .2f [% .2f;% .2f]", cor_stirn["estimate"], cor_stirn["conf.low"], cor_stirn["conf.high"])
    
    cor_imp <- broom::tidy(cor.test(df[[paste0("est_", Hormone_name, "_bc")]], log(df$hormone)))
    results$r_bc_imputed <- sprintf("% .2f [% .2f;% .2f]", cor_imp["estimate"], cor_imp["conf.low"], cor_imp["conf.high"])
    
    cor_imp <- broom::tidy(cor.test(df[[paste0("est_", Hormone_name, "_bc")]], df$log_hormone_diff))
    results$r_diff_bc_imputed <- sprintf("% .2f [% .2f;% .2f]", cor_imp["estimate"], cor_imp["conf.low"], cor_imp["conf.high"])
    sd_unrestricted <- sd(bc_days[which(bc_days$bc_day > -30), paste0("est_", Hormone_name, "_bc")])
    sd_observed <- sd(df[[paste0("est_", Hormone_name, "_bc")]], na.rm = T)
    results$sd_bc_imputed <- sprintf("%.2f (%.0f%%)", sd_observed, 100*sd_observed^2/sd_unrestricted^2)
    rr_cor <- psychometric::cRRr(cor_imp[,c("estimate", "conf.low", "conf.high")],
                                 sd_observed, sd_unrestricted)
    results$r_diff_bc_imputed_rr <- sprintf("% .2f [% .2f;% .2f]", rr_cor["unrestricted.estimate"], rr_cor["unrestricted.conf.low"], rr_cor["unrestricted.conf.high"])
    
    
    factor <- compute_factor(df$bc_day, df[[paste0("est_", Hormone_name, "_bc")]], df$hormone, threshold = c(-10, -2))
    results$imputed_bc_vs_measured_graph <- ggplot(df, aes_string("bc_day",
                                            y = paste0("est_", Hormone_name, "_bc"),
                                            ymin = paste0("est_", Hormone_name, "_bc_low"),
                                            ymax = paste0("est_", Hormone_name, "_bc_high"))) +
      geom_smooth(stat = "identity", color = "#6F906F", fill = "#6F906F") +
      geom_point(aes(y = log(hormone*factor)), alpha = 0.2) +
      ylab(paste0("log(", Hormone_name," × ",sprintf("%.1f", factor),")")) +
      scale_x_continuous("Days until next menstrual onset", limits = c(-30, 0))
    
    
    
    results$rmse_bc_imputed <- sprintf("% .2f", rmse(log(df$hormone * factor) - df[[paste0("est_", Hormone_name, "_bc")]]))
    
  }
  
  if(n_nonmissing(df$fc_day)>20) {
    results$fc_day_model <- fc_day_model
    r2 <- bayes_R2(fc_day_model, re_formula = NA)
    results$r_fc <- sprintf("%.2f [%.2f;%.2f]", sqrt(r2[,"Estimate"]), sqrt(r2[,"Q2.5"]), sqrt(r2[,"Q97.5"])) 
    
    results$rmse_fc <- rmse_brms(fc_day_model)
    
    
    r2 <- loo_R2(fc_day_model, re_formula = NA)
    results$r_loo_fc <- sprintf("%.2f [%.2f;%.2f]", sqrt(r2[,"Estimate"]), sqrt(r2[,"Q2.5"]), sqrt(r2[,"Q97.5"])) 
    
    cor_stirn <- broom::tidy(cor.test(df$prc_stirn_fc, df$hormone))
    results$r_fc_stirn <- sprintf("% .2f [% .2f;% .2f]", cor_stirn["estimate"], cor_stirn["conf.low"], cor_stirn["conf.high"])

    cor_stirn <- broom::tidy(cor.test(df$prc_stirn_fc, log(df$hormone)))
    results$r_log_fc_stirn <- sprintf("% .2f [% .2f;% .2f]", cor_stirn["estimate"], cor_stirn["conf.low"], cor_stirn["conf.high"])

    cor_stirn <- broom::tidy(cor.test(df$prc_stirn_fc, df$hormone_diff))
    results$r_diff_fc_stirn <- sprintf("% .2f [% .2f;% .2f]", cor_stirn["estimate"], cor_stirn["conf.low"], cor_stirn["conf.high"])
    
    cor_stirn <- broom::tidy(cor.test(df$prc_stirn_fc, df$log_hormone_diff))
    results$r_log_diff_fc_stirn <- sprintf("% .2f [% .2f;% .2f]", cor_stirn["estimate"], cor_stirn["conf.low"], cor_stirn["conf.high"])
    
    cor_imp <- broom::tidy(cor.test(df[[paste0("est_", Hormone_name, "_fc")]], log(df$hormone)))
    results$r_fc_imputed <- sprintf("% .2f [% .2f;% .2f]", cor_imp["estimate"], cor_imp["conf.low"], cor_imp["conf.high"])
    
    cor_imp <- broom::tidy(cor.test(df[[paste0("est_", Hormone_name, "_fc")]], df$log_hormone_diff))
    results$r_diff_fc_imputed <- sprintf("% .2f [% .2f;% .2f]", cor_imp["estimate"], cor_imp["conf.low"], cor_imp["conf.high"])
    sd_unrestricted <- sd(fc_days[which(fc_days$fc_day < 29), paste0("est_", Hormone_name, "_fc")])
    sd_observed <- sd(df[[paste0("est_", Hormone_name, "_fc")]], na.rm = T)
    results$sd_fc_imputed <- sprintf("%.2f (%.0f%%)", sd_observed, 100*sd_observed^2/sd_unrestricted^2)
    rr_cor <- psychometric::cRRr(cor_imp[,c("estimate", "conf.low", "conf.high")],
                                 sd_observed, sd_unrestricted)
    results$r_diff_fc_imputed_rr <- sprintf("% .2f [% .2f;% .2f]", rr_cor["unrestricted.estimate"], rr_cor["unrestricted.conf.low"], rr_cor["unrestricted.conf.high"])
    
    factor <- compute_factor(df$fc_day, df[[paste0("est_", Hormone_name, "_fc")]], df$hormone, threshold = c(20, 28))
    
    results$imputed_fc_vs_measured_graph <- ggplot(df, aes_string("fc_day",
                                                                  y = paste0("est_", Hormone_name, "_fc"),
                                                                  ymin = paste0("est_", Hormone_name, "_fc_low"),
                                                                  ymax = paste0("est_", Hormone_name, "_fc_high"))) +
      geom_smooth(stat = "identity", color = "#6F906F", fill = "#6F906F") +
      geom_point(aes(y = log(hormone*factor)), alpha = 0.2) +
      ylab(paste0("log(", Hormone_name," × ",sprintf("%.0f", factor),")")) +
      scale_x_continuous("Days since last menstrual onset", limits = c(0, 30))
    
    
    
    results$rmse_fc_imputed <- sprintf("% .2f", rmse(log(df$hormone * factor) - df[[paste0("est_", Hormone_name, "_fc")]]))
    
  }
  
  if(n_nonmissing(df$lh_day)>20) {
    results$lh_day_model <- lh_day_model
    r2 <- bayes_R2(lh_day_model, re_formula = NA)
    results$r_lh <- sprintf("%.2f [%.2f;%.2f]", sqrt(r2[,"Estimate"]), sqrt(r2[,"Q2.5"]), sqrt(r2[,"Q97.5"])) 
    
    results$rmse_lh <- rmse_brms(lh_day_model)
    
    r2 <- loo_R2(lh_day_model, re_formula = NA)
    results$r_loo_lh <- sprintf("%.2f [%.2f;%.2f]", sqrt(r2[,"Estimate"]), sqrt(r2[,"Q2.5"]), sqrt(r2[,"Q97.5"])) 
    
    cor_stirn <- broom::tidy(cor.test(df$fertile_lh, df$hormone))
    results$r_prob_lh <- sprintf("% .2f [% .2f;% .2f]", cor_stirn["estimate"], cor_stirn["conf.low"], cor_stirn["conf.high"])

    cor_stirn <- broom::tidy(cor.test(df$fertile_lh, log(df$hormone)))
    results$r_log_prob_lh <- sprintf("% .2f [% .2f;% .2f]", cor_stirn["estimate"], cor_stirn["conf.low"], cor_stirn["conf.high"])
    
    cor_stirn <- broom::tidy(cor.test(df$fertile_lh, df$hormone_diff))
    results$r_diff_prob_lh <- sprintf("% .2f [% .2f;% .2f]", cor_stirn["estimate"], cor_stirn["conf.low"], cor_stirn["conf.high"])
    
    cor_stirn <- broom::tidy(cor.test(df$fertile_lh, df$log_hormone_diff))
    results$r_log_diff_prob_lh <- sprintf("% .2f [% .2f;% .2f]", cor_stirn["estimate"], cor_stirn["conf.low"], cor_stirn["conf.high"])
    
    
    cor_imp <- broom::tidy(cor.test(df[[paste0("est_", Hormone_name, "_lh")]], log(df$hormone)))
    results$r_lh_imputed <- sprintf("% .2f [% .2f;% .2f]", cor_imp["estimate"], cor_imp["conf.low"], cor_imp["conf.high"])
    
    cor_imp <- broom::tidy(cor.test(df[[paste0("est_", Hormone_name, "_lh")]], df$log_hormone_diff))
    results$r_diff_lh_imputed <- sprintf("% .2f [% .2f;% .2f]", cor_imp["estimate"], cor_imp["conf.low"], cor_imp["conf.high"])
    sd_unrestricted <- sd(lh_days[which(between(lh_days$lh_day, -15, 13)),paste0("est_", Hormone_name, "_lh")])
    sd_observed <- sd(df[[paste0("est_", Hormone_name, "_lh")]], na.rm = T)
    results$sd_lh_imputed <- sprintf("%.2f (%.0f%%)", sd_observed, 100*sd_observed^2/sd_unrestricted^2)
    rr_cor <- psychometric::cRRr(cor_imp[,c("estimate", "conf.low", "conf.high")],
                                 sd_observed, sd_unrestricted)
    results$r_diff_lh_imputed_rr <- sprintf("% .2f [% .2f;% .2f]", rr_cor["unrestricted.estimate"], rr_cor["unrestricted.conf.low"], rr_cor["unrestricted.conf.high"])
    
    factor <- compute_factor(df$lh_day, df[[paste0("est_", Hormone_name, "_lh")]], df$hormone, threshold = c(2, 10))
    
    results$imputed_lh_vs_measured_graph <- ggplot(df, aes_string("lh_day",
                                                                  y = paste0("est_", Hormone_name, "_lh"),
                                                                  ymin = paste0("est_", Hormone_name, "_lh_low"),
                                                                  ymax = paste0("est_", Hormone_name, "_lh_high"))) +
      geom_smooth(stat = "identity", color = "#6F906F", fill = "#6F906F") +
      geom_point(aes(y = log(hormone*factor)), alpha = 0.2) +
      ylab(paste0("log(", Hormone_name," × ",sprintf("%.0f", factor),")")) +
      scale_x_continuous("Day relative to LH-surge", limits = c(-15, 15))
    
    results$rmse_lh_imputed <- sprintf("% .2f", rmse(log(df$hormone * factor) - df[[paste0("est_", Hormone_name, "_lh")]]))
  }
  
  rio::export(results, paste0("summaries/table_", Dataset, "_", Hormone, ".rds"))
  results
}


summarise_hormones <- function(df, Dataset, Method = "Salivary Immunoassay") {
  # df
  # df$id
  # df$cycle
  # df$bc_day
  # df$fc_day
  # df$cycle_length
  # df$estradiol
  # df$progesterone
  
  results <- list()
  results$Dataset <- Dataset
  results$Hormone <- "Estradiol & Progesterone"
  results$Method <- Method
  results$scatterplot <- ggplot(df, aes(estradiol, progesterone)) + geom_point(alpha = 0.1)  + scale_x_log10() + scale_y_log10() + xlab("log(Estradiol)") + ylab("log(Progesterone)")
  results$n_women <- n_distinct(df$id)
  results$n_cycles <- n_distinct(paste(df$id, df$cycle))
  results$n_days <- nrow(df)
  
  df$ratio <- df$estradiol/df$progesterone
  
  cor_hormones <- broom::tidy(cor.test(df$estradiol, df$progesterone))
  results$r_ep <- sprintf("% .2f [% .2f;% .2f]", cor_hormones["estimate"], cor_hormones["conf.low"], cor_hormones["conf.high"])
  
  cor_hormones <- broom::tidy(cor.test(log(df$estradiol), log(df$progesterone)))
  results$r_log_ep <- sprintf("% .2f [% .2f;% .2f]", cor_hormones["estimate"], cor_hormones["conf.low"], cor_hormones["conf.high"])
  
  cor_hormones <- broom::tidy(cor.test((df$estradiol), (df$ratio)))
  results$r_e_ratio <- sprintf("% .2f [% .2f;% .2f]", cor_hormones["estimate"], cor_hormones["conf.low"], cor_hormones["conf.high"])
  
  
  cor_hormones <- broom::tidy(cor.test(log(df$estradiol), log(df$ratio)))
  results$r_log_e_ratio <- sprintf("% .2f [% .2f;% .2f]", cor_hormones["estimate"], cor_hormones["conf.low"], cor_hormones["conf.high"])
  
  cor_hormones <- broom::tidy(cor.test((df$progesterone), (df$ratio)))
  results$r_p_ratio <- sprintf("% .2f [% .2f;% .2f]", cor_hormones["estimate"], cor_hormones["conf.low"], cor_hormones["conf.high"])
  
  
  cor_hormones <- broom::tidy(cor.test(log(df$progesterone), log(df$ratio)))
  results$r_log_p_ratio <- sprintf("% .2f [% .2f;% .2f]", cor_hormones["estimate"], cor_hormones["conf.low"], cor_hormones["conf.high"])
  
  df <- df %>% filter(!is.na(estradiol), !is.na(progesterone), 
                      is.na(cycle_length) | between(cycle_length, 20,35))
  results$usable_n <- nrow(df)
  
  bc_days <- readRDS("merge_files/bc_days.rds")
  fc_days <- readRDS("merge_files/fc_days.rds")
  
  bc_days = bc_days %>% select(bc_day, prc_stirn_bc)
  df = left_join(df, bc_days)
  
  fc_days = fc_days %>% select(fc_day, prc_stirn_fc)
  df = left_join(df, fc_days, by = "fc_day")
  
  # from Jünger/Stern et al. 2018 Supplementary Material
  # Day relative to ovulation	Schwartz et al., (1980)	Wilcox et al., (1998)	Colombo & Masarotto (2000)	Weighted average
  lh_days <- readRDS("merge_files/lh_days.rds")
  
  if(is.null(df$lh_day)) df$lh_day <- NA_real_
  df = df %>% left_join(lh_days, by = "lh_day")
  
  several_cycles <- { df %>% group_by(id) %>% filter(n_distinct(cycle) > 1) %>% nrow() } > 0
  several_timepoints <- { df %>% group_by(id) %>% filter(n() > 1) %>% nrow() } > 0
  
  
  if(n_nonmissing(df$bc_day)>20) {
    
    cor_hormones <- broom::tidy(cor.test(df$prc_stirn_bc, (df$ratio)))
    results$r_bc_ratio <- sprintf("% .2f [% .2f;% .2f]", cor_hormones["estimate"], cor_hormones["conf.low"], cor_hormones["conf.high"])
    
    cor_hormones <- broom::tidy(cor.test(df$prc_stirn_bc, log(df$ratio)))
    results$r_bc_log_ratio <- sprintf("% .2f [% .2f;% .2f]", cor_hormones["estimate"], cor_hormones["conf.low"], cor_hormones["conf.high"])
    
    (m_bc_prob_rat <- brm(prc_stirn_bc ~ log(estradiol) + log(ratio), data = df , file_refit = "on_change",
                      control = list(adapt_delta = 0.99),
                      file = paste0("models/m_", Dataset, "_Er_bc")) %>% 
        add_criterion("loo_R2", re_formula = NA) %>%
        add_criterion("bayes_R2", re_formula = NA))
    
    (m_bc_prob <- brm(prc_stirn_bc ~ s(log(estradiol), log(progesterone)), data = df , file_refit = "on_change",
                      control = list(adapt_delta = 0.99),
                      file = paste0("models/m_", Dataset, "_ExP_bc")) %>% 
        add_criterion("loo_R2", re_formula = NA) %>%
        add_criterion("bayes_R2", re_formula = NA))
    
    results$m_bc_prob_rat <- m_bc_prob_rat
    r2 <- loo_R2(m_bc_prob_rat, re_formula = NA)
    results$r_erat_bc <- sprintf("%.2f [%.2f;%.2f]", sqrt(r2[,"Estimate"]), sqrt(r2[,"Q2.5"]), sqrt(r2[,"Q97.5"])) 
    
    results$m_bc_prob <- m_bc_prob

    r2 <- loo_R2(m_bc_prob, re_formula = NA)
    results$r_loo_bc <- sprintf("%.2f [%.2f;%.2f]", sqrt(r2[,"Estimate"]), sqrt(r2[,"Q2.5"]), sqrt(r2[,"Q97.5"])) 
  }
  
  if(n_nonmissing(df$fc_day)>20) {
    
    cor_hormones <- broom::tidy(cor.test(df$prc_stirn_fc, (df$ratio)))
    results$r_fc_ratio <- sprintf("% .2f [% .2f;% .2f]", cor_hormones["estimate"], cor_hormones["conf.low"], cor_hormones["conf.high"])
    
    cor_hormones <- broom::tidy(cor.test(df$prc_stirn_fc, log(df$ratio)))
    results$r_fc_log_ratio <- sprintf("% .2f [% .2f;% .2f]", cor_hormones["estimate"], cor_hormones["conf.low"], cor_hormones["conf.high"])
    
    
    (m_fc_prob_rat <- brm(prc_stirn_fc ~ log(estradiol) + log(ratio), data = df , file_refit = "on_change",
                          control = list(adapt_delta = 0.99),
                          file = paste0("models/m_", Dataset, "_Er_fc")) %>% 
        add_criterion("loo_R2", re_formula = NA) %>%
        add_criterion("bayes_R2", re_formula = NA))
    
    (m_fc_prob <- brm(prc_stirn_fc ~ s(log(estradiol), log(progesterone)), data = df , file_refit = "on_change",
                      control = list(adapt_delta = 0.99),
                      file = paste0("models/m_", Dataset, "_ExP_fc")) %>% 
        add_criterion("loo_R2", re_formula = NA) %>%
        add_criterion("bayes_R2", re_formula = NA))
    
    results$m_fc_prob_rat <- m_fc_prob_rat
    r2 <- loo_R2(m_fc_prob_rat, re_formula = NA)
    results$r_erat_fc <- sprintf("%.2f [%.2f;%.2f]", sqrt(r2[,"Estimate"]), sqrt(r2[,"Q2.5"]), sqrt(r2[,"Q97.5"])) 
    
    results$m_fc_prob <- m_fc_prob
    
    r2 <- loo_R2(m_fc_prob, re_formula = NA)
    results$r_loo_fc <- sprintf("%.2f [%.2f;%.2f]", sqrt(r2[,"Estimate"]), sqrt(r2[,"Q2.5"]), sqrt(r2[,"Q97.5"])) 
  }
  
  if(n_nonmissing(df$lh_day)>20) {
    
    cor_hormones <- broom::tidy(cor.test(df$fertile_lh, (df$ratio)))
    results$r_lh_ratio <- sprintf("% .2f [% .2f;% .2f]", cor_hormones["estimate"], cor_hormones["conf.low"], cor_hormones["conf.high"])
    
    cor_hormones <- broom::tidy(cor.test(df$fertile_lh, log(df$ratio)))
    results$r_lh_log_ratio <- sprintf("% .2f [% .2f;% .2f]", cor_hormones["estimate"], cor_hormones["conf.low"], cor_hormones["conf.high"])  
    
    (m_lh_prob_rat <- brm(fertile_lh ~ log(estradiol) + log(ratio), data = df, file_refit = "on_change",
                          control = list(adapt_delta = 0.99),
                          file = paste0("models/m_", Dataset, "_Er_lh")) %>% 
        add_criterion("loo_R2", re_formula = NA) %>%
        add_criterion("bayes_R2", re_formula = NA))
    
    results$m_lh_prob_rat <- m_lh_prob_rat
    r2 <- loo_R2(m_lh_prob_rat, re_formula = NA)
    results$r_erat_lh <- sprintf("%.2f [%.2f;%.2f]", sqrt(r2[,"Estimate"]), sqrt(r2[,"Q2.5"]), sqrt(r2[,"Q97.5"])) 
    
    (m_lh_prob <- brm(fertile_lh ~ s(log(estradiol),log(progesterone)), data = df , file_refit = "on_change",
                      control = list(adapt_delta = 0.99),
                      file = paste0("models/m_", Dataset, "_ExP_lh")) %>% 
        add_criterion("loo_R2", re_formula = NA) %>%
        add_criterion("bayes_R2", re_formula = NA))
    
    # (m_lh_prob_gp <- brm(fertile_lh ~ gp(log(estradiol), log(progesterone)), data = df , file_refit = "on_change",
    #  file = paste0("models/m_", Dataset, "_gpEP_lh")) %>% 
    #     add_criterion("loo_R2", re_formula = NA) %>%
    #     add_criterion("bayes_R2", re_formula = NA))
    
    results$m_lh_prob <- m_lh_prob

    r2 <- loo_R2(m_lh_prob, re_formula = NA)
    results$r_loo_lh <- sprintf("%.2f [%.2f;%.2f]", sqrt(r2[,"Estimate"]), sqrt(r2[,"Q2.5"]), sqrt(r2[,"Q97.5"]))
  }
  
  
  rio::export(results, paste0("summaries/table_", Dataset, "_ratio.rds"))
  results
}


export_anon <- function(filename, data) {
  set.seed(05102019)
  shareable <- data %>% select(id, cycle, cycle_length, fc_day, one_of("bc_day", "lh_day"), estradiol, estradiol_cens, progesterone, progesterone_cens) %>% 
    ungroup() %>% 
    arrange(rnorm(n())) %>% 
    mutate(id = factor(as.numeric(forcats::fct_inorder(factor(id))))) %>% 
    arrange(id, cycle, fc_day)
  
  
  var_nam <- list(
    id = "Random Participant ID",
    cycle = "Cycle number, incrementing within-subject",
    cycle_length = "Cycle length in days",
    fc_day = "Forward-counted cycle day. Counted from last (recalled) menstrual onset, starting at zero.",
    bc_day = "Backward-counted cycle day. Counted from next observed menstrual onset (0), beginning with -1.",
    lh_day = "Day relative to LH surge",
    estradiol = "Salivary estradiol in pg/ml",
    estradiol_cens = "Whether estradiol reached the limit of detection (left) or not (none)",
    progesterone = "Salivary progesterone in pg/ml",
    progesterone_cens = "Whether progesterone reached the limit of detection (left) or not (none)"
  )
  labelled::var_label(shareable) <- var_nam[names(shareable)]
  
  rio::export(shareable, paste0("osf_data/", filename, ".tsv"))
  rio::export(shareable, paste0("osf_data/", filename, ".rds"))
  rio::export(shareable, paste0("osf_data/", filename, ".sav"))
}
