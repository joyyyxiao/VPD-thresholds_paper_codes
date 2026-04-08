library(mgcv)
library(raster)
library(ggplot2)
library(car)
#load data
data <- read.csv("suppl_data/data_multiyearmean.csv",row.names = 1)


#load eCO2 effect data
total_eco2 <- read.csv("results_mr1/FLUXCOM_gpp_total_CO2_effect.csv")
total_eco2 <- read.csv("results_mr1/GLASS_gpp_total_CO2_effect.csv")
total_eco2 <- read.csv("results_mr1/MODIS_gpp_total_CO2_effect.csv")

clean_outliers <- function(x, probs = c(0.01, 0.99)) {
  
  x[!is.finite(x)] <- NA
  qs <- quantile(x, probs = probs, na.rm = TRUE)
  x[x < qs[1] | x > qs[2]] <- NA
  
  return(x)
}

data$eco2_re <- clean_outliers(total_eco2$RelativeChange)
data$eco2_ab <- clean_outliers(total_eco2$AbsoluteChange)

#set factors
data$Biome <- factor(data$Biome)
data$Myc <- factor(data$Myc)
#data00 <- na.omit(data)

# Load land-use change masks at different percentage thresholds 
# (the percentage indicates the proportion of each pixel that experienced land-use change during 2000–2022)
mask_20 <- raster("suppl_data/lulcc_0022_20perc.tif")
mask_20[mask_20 == 1] <- NA

data$lccmask <- as.numeric(as.matrix(mask_20))
data20 <- na.omit(data)

#smooth check
check_smooth <- function(x, name, top_n = 5) {
  x2 <- x[is.finite(x)]
  n  <- length(x2)
  u  <- length(unique(x2))
  prop_unique <- u / n
  
  qs <- quantile(x2, probs = seq(0.01, 0.99, by = 0.01), na.rm = TRUE)
  dup_q <- sum(duplicated(qs))
  
  tab <- sort(table(x2), decreasing = TRUE)
  top_share <- sum(tab[1:min(top_n, length(tab))]) / n
  
  data.frame(
    var = name,
    n = n,
    n_unique = u,
    prop_unique = prop_unique,
    duplicated_quantiles = dup_q,
    top_share = top_share,
    stringsAsFactors = FALSE
  )
}

vars <- c("MAP","TMP","PET","AET","MI","VPD","SM","PAR","CNr","P")
diag_tab <- do.call(rbind, lapply(vars, function(v) check_smooth(data20[[v]], v)))
diag_tab[order(diag_tab$prop_unique), ]

#family
hist(data20$eco2_re)
hist(data20$eco2_ab)

#Multicollinearity
cor_mat <- cor(data[,c("MAP","P","TMP","PET","AET","MI","VPD","SM","PAR","CNr")], use = "complete.obs")
corrplot::corrplot(cor_mat, method="color")


#lm_tmp <- lm(eco2 ~ MAP + TMP + PET + AET + MI + VPD + SM + PAR + CNr, data=data20)
#vif(lm_tmp)

###function###
vif_subsets <- function(response, predictors, data,
                                     require = c("VPD","SM"),
                                     vif_th = 5,
                                     min_extra = 0,
                                     max_extra = NULL) {
  if (!requireNamespace("car", quietly = TRUE)) {
    stop("Please install the 'car' package first.")
  }
  stopifnot(is.character(response), length(response) == 1)
  
  needed <- unique(c(response, predictors, require))
  miss <- setdiff(needed, names(data))
  if (length(miss) > 0) stop("Missing columns in data: ", paste(miss, collapse = ", "))
  
 
  predictors <- unique(predictors)
  if (!all(require %in% predictors)) {
    stop("All `require` variables must be included in `predictors`.")
  }
  
  extras <- setdiff(predictors, require)
  if (is.null(max_extra)) max_extra <- length(extras)
  max_extra <- min(max_extra, length(extras))
  min_extra <- max(0, min_extra)
  if (min_extra > max_extra) stop("min_extra must be <= max_extra")
  
  out <- list()
  idx <- 0L
  
  for (k in seq(min_extra, max_extra)) {
    combs <- if (k == 0) list(character(0)) else combn(extras, k, simplify = FALSE)
    
    for (add_vars in combs) {
      vars <- c(require, add_vars)
      fml <- as.formula(paste(response, "~", paste(vars, collapse = " + ")))
      
      m <- tryCatch(stats::lm(fml, data = data), error = function(e) NULL)
      if (is.null(m)) next
      
      v <- tryCatch(car::vif(m), error = function(e) NULL)
      if (is.null(v)) next
      
      max_v <- suppressWarnings(max(v, na.rm = TRUE))
      if (!is.finite(max_v)) next
      
      if (max_v < vif_th) {
        idx <- idx + 1L
        out[[idx]] <- data.frame(
          n_var   = length(vars),
          n_extra = length(add_vars),
          max_vif = max_v,
          vars    = paste(vars, collapse = ", "),
          formula_lm = paste(response, "~", paste(vars, collapse = " + ")),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(out) == 0) {
    return(data.frame(n_var=integer(), n_extra=integer(), max_vif=numeric(),
                      vars=character(), formula_lm=character()))
  }
  
  res <- do.call(rbind, out)
  res <- res[order(res$n_extra, res$max_vif), ]
  rownames(res) <- NULL
  res
}

fit_gam_comparison <- function(
    vif_sets,
    data,
    response,
    base_vars_required = c("VPD", "SM"),
    k = 6,
    # if allow Biome / Myc
    allow_biome = TRUE,
    allow_myc   = TRUE,
    # if allow interaction
    allow_interaction = TRUE,
    # set P 
    p_var = "P",
    # smooth 
    smooth_fun = function(x, k) sprintf("s(%s, k=%d)", x, k),
    family = Gamma(link="log"),
    method = "ML",
    verbose = TRUE
) {
  start_time <- Sys.time()
  cat("Start time:", format(start_time), "\n")
  
  total_iter <- nrow(vif_sets) *
    (if (allow_biome) 2 else 1) *
    (if (allow_myc) 2 else 1) *
    (if (allow_interaction) 2 else 1)
  
  pb <- txtProgressBar(min = 0, max = total_iter, style = 3)
  on.exit(close(pb), add = TRUE)
  progress <- 0
  
  stopifnot(is.data.frame(vif_sets))
  if (!("vars" %in% names(vif_sets))) stop("`vif_sets` must contain a column named 'vars'")
  if (!(response %in% names(data))) stop("The specified response variable is not found in `data`:  ", response)
  
  #  "VPD, SM, MAP" -> c("VPD","SM","MAP")
  parse_vars <- function(s) {
    trimws(unlist(strsplit(as.character(s), ",")))
  }
  
  # Biome/Myc 
  biome_opts <- if (allow_biome) c(FALSE, TRUE) else FALSE
  myc_opts   <- if (allow_myc)   c(FALSE, TRUE) else FALSE
  int_opts   <- if (allow_interaction) c(FALSE, TRUE) else FALSE
  
  out <- list()
  idx <- 0
  
  for (r in seq_len(nrow(vif_sets))) {
    
    vv <- parse_vars(vif_sets$vars[r])
    
    for (req in base_vars_required) {
      if (!req %in% vv) vv <- c(req, vv)
    }
    vv <- unique(vv)
    
    vv <- vv[vv %in% names(data)]
    if (length(vv) == 0) next
    
    cont_vars <- setdiff(vv, c("Biome", "Myc"))
    
    main_terms <- c()
    for (v in cont_vars) {
      if (v == p_var) {
        main_terms <- c(main_terms, v)           
      } else {
        main_terms <- c(main_terms, smooth_fun(v, k))
      }
    }
    
    make_formula <- function(add_biome, add_myc, add_int) {
      terms <- main_terms
      
      if (add_int) {
        if (!any(grepl("^s\\(VPD", terms))) terms <- c(terms, smooth_fun("VPD", k))
        if (!any(grepl("^s\\(SM",  terms))) terms <- c(terms, smooth_fun("SM",  k))
        terms <- c(terms, sprintf("ti(VPD, SM, k=c(%d,%d))", k, k))
      }
      
      if (add_biome && ("Biome" %in% names(data))) terms <- c(terms, "Biome")
      if (add_myc   && ("Myc"   %in% names(data))) terms <- c(terms, "Myc")
      
      terms <- unique(terms)
      
      as.formula(paste(response, "~", paste(terms, collapse = " + ")))
    }
    
    for (add_biome in biome_opts) {
      for (add_myc in myc_opts) {
        for (add_int in int_opts) {
          progress <- progress + 1
          setTxtProgressBar(pb, progress)
          
          fml <- make_formula(add_biome, add_myc, add_int)
          fml_txt <- paste(deparse(fml), collapse = "")
          fit <- tryCatch(
            mgcv::gam(fml, data = data, family = family, method = method),
            error = function(e) e
          )
          
          idx <- idx + 1
          
          if (inherits(fit, "error")) {
            out[[idx]] <- data.frame(
              set_id = r,
              vars = paste(vv, collapse = ", "),
              add_biome = add_biome,
              add_myc = add_myc,
              add_int = add_int,
              formula_gam = fml_txt,
              AIC = NA_real_,
              edf = NA_real_,
              dev_expl = NA_real_,
              ok = FALSE,
              err = fit$message,
              stringsAsFactors = FALSE
            )
            if (verbose) message("failed: ", fml_txt, " | ", fit$message)
          } else {
            sm <- summary(fit)
            out[[idx]] <- data.frame(
              set_id = r,
              vars = paste(vv, collapse = ", "),
              add_biome = add_biome,
              add_myc = add_myc,
              add_int = add_int,
              formula_gam = fml_txt,
              AIC = AIC(fit),
              edf = sum(fit$edf),
              dev_expl = sm$dev.expl,
              ok = TRUE,
              err = NA_character_,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
  
  res <- do.call(rbind, out)
  res <- res[order(res$AIC), ]
  rownames(res) <- NULL
  end_time <- Sys.time()
  cat("\nEnd time:", format(end_time), "\n")
  cat("Total time:", round(difftime(end_time, start_time, units = "mins"), 2), "mins\n")
  return(res)
}


#####################

preds <- c("MAP","TMP","PET","AET","MI","VPD","SM","PAR","CNr","P")
response <- "eco2_re"
response <- "eco2_ab"

vif_sets <- vif_subsets(
  response   = response,
  predictors = preds,
  data       = data20,
  require    = c("VPD","SM"),
  vif_th     = 5,
  min_extra  = 0
)

vif_sets

gam_aic_table <- fit_gam_comparison(
  vif_sets = vif_sets,
  data = data20,          
  response = response,
  k = 6,
  allow_biome = TRUE,
  allow_myc   = TRUE,
  allow_interaction = TRUE,
  p_var = "P",
  method = "ML"           
)

top10_table <- gam_aic_table %>%
  filter(ok) %>%
  mutate(
    delta_AIC = AIC - min(AIC),
    rel_lik = exp(-0.5 * delta_AIC),
    wi = rel_lik / sum(rel_lik),
    K = round(edf, 0)
  ) %>%
  arrange(delta_AIC) %>%
  slice(1:10) %>%
  mutate(
    Model = formula_gam,
    AIC = round(AIC, 1),
    delta_AIC = round(delta_AIC, 2),
    wi = round(wi, 2)
  ) %>%
  select(Model, K, AIC, delta_AIC, wi) %>%
  rename(`ΔAIC` = delta_AIC, `w_i` = wi)


write.csv(gam_aic_table,"results_mr1/MODIS_gpp_re_gam_aic_table.csv")
write.csv(top10_table,"results_mr1/MODIS_gpp_re_gam_aic_table_t10.csv")
write.csv(gam_aic_table,"results_mr1/MODIS_gpp_ab_gam_aic_table.csv")
write.csv(top10_table,"results_mr1/MODIS_gpp_ab_gam_aic_table_t10.csv")

write.csv(gam_aic_table,"results_mr1/GLASS_gpp_re_gam_aic_table.csv")
write.csv(top10_table,"results_mr1/GLASS_gpp_re_gam_aic_table_t10.csv")
write.csv(gam_aic_table,"results_mr1/GLASS_gpp_ab_gam_aic_table.csv")
write.csv(top10_table,"results_mr1/GLASS_gpp_ab_gam_aic_table_t10.csv")

write.csv(gam_aic_table,"results_mr1/FLUXCOM_gpp_re_gam_aic_table.csv")
write.csv(top10_table,"results_mr1/FLUXCOM_gpp_re_gam_aic_table_t10.csv")
write.csv(gam_aic_table,"results_mr1/FLUXCOM_gpp_ab_gam_aic_table.csv")
write.csv(top10_table,"results_mr1/FLUXCOM_gpp_ab_gam_aic_table_t10.csv")


#gam model
f_gam_re <- gam(
  eco2_re ~
    s(VPD, k=6) +
    s(SM,  k=6) +
    s(AET, k=6) +
    s(MI, k=6) +
    s(PAR, k=6) +
    s(CNr, k=6) +
    ti(VPD, SM, k=c(6,6))+
    P+
    Biome + Myc,
  data = data20,
  family = Gamma(link="log"),
  method = "REML"
)


f_gam_ab <- gam(
  eco2_ab ~
    s(VPD, k=6) +
    s(SM,  k=6) +
    s(AET, k=6) +
    s(MAP, k=6) +
    s(PET, k=6) +
    ti(VPD, SM, k=c(6,6))+
    P+
    Biome + Myc,
  data = data20,
  family = Gamma(link="log"),
  method = "REML"
)

g_gam_re <- gam(
  eco2_re ~
    s(VPD, k=6) +
    s(SM,  k=6) +
    s(AET, k=6) +
    s(MI, k=6) +
    s(PAR, k=6) +
    s(CNr, k=6) +
    ti(VPD, SM, k=c(6,6))+
    P+
    Biome + Myc,
  data = data20,
  family = Gamma(link="log"),
  method = "REML"
)

g_gam_ab <- gam(
  eco2_ab ~
    s(VPD, k=6) +
    s(SM,  k=6) +
    s(AET, k=6) +
    s(MI, k=6) +
    s(PAR, k=6) +
    s(CNr, k=6) +
    ti(VPD, SM, k=c(6,6))+
    P+
    Biome + Myc,
  data = data20,
  family = Gamma(link="log"),
  method = "REML"
)


m_gam_re <- gam(
  eco2_re ~
    s(VPD, k=6) +
    s(SM,  k=6) +
    s(AET, k=6) +
    s(MI, k=6) +
    s(PAR, k=6) +
    s(CNr, k=6) +
    ti(VPD, SM, k=c(6,6))+
    P+
    Biome + Myc,
  data = data20,
  family = Gamma(link="log"),
  method = "REML"
)

m_gam_ab <- gam(
  eco2_ab ~
    s(VPD, k=6) +
    s(SM,  k=6) +
    s(AET, k=6) +
    s(MI, k=6) +
    s(PAR, k=6) +
    s(CNr, k=6) +
    ti(VPD, SM, k=c(6,6))+
    Biome + Myc,
  data = data20,
  family = Gamma(link="log"),
  method = "REML"
)


#############
mods <- list(
  f_gam_ab = f_gam_ab,
  f_gam_re = f_gam_re,
  g_gam_ab = g_gam_ab,
  g_gam_re = g_gam_re,
  m_gam_ab = m_gam_ab,
  m_gam_re = m_gam_re
)

# 
fit_tab <- lapply(names(mods), function(nm) {
  sm <- summary(mods[[nm]])
  data.frame(
    model = nm,
    adj_r2 = sm$r.sq,
    dev_expl = sm$dev.expl
  )
})
fit_tab <- do.call(rbind, fit_tab)
fit_tab

# 
smooth_tab <- lapply(names(mods), function(nm) {
  sm <- summary(mods[[nm]])
  out <- as.data.frame(sm$s.table)
  out$term <- rownames(out)
  out$model <- nm
  rownames(out) <- NULL
  out
})
smooth_tab <- do.call(rbind, smooth_tab)
smooth_tab

#model diagnostics
gam.check(f_gam_ab)
gam.check(f_gam_re)
gam.check(g_gam_ab)
gam.check(g_gam_re)
gam.check(m_gam_ab)
gam.check(m_gam_re)

#s(VPD)
library(gratia)

sm_list <- lapply(names(mods), function(nm) {
  df <- smooth_estimates(mods[[nm]], smooth = "s(VPD)")
  df$model <- nm
  df
})
sm_df <- do.call(rbind, sm_list)
sm_df$type <- ifelse(grepl("_ab$", sm_df$model), "Absolute", "Relative")
sm_df$dataset <- NA
sm_df$dataset[grepl("^f_", sm_df$model)] <- "FLUXCOM GPP"
sm_df$dataset[grepl("^g_", sm_df$model)] <- "GLASS GPP"
sm_df$dataset[grepl("^m_", sm_df$model)] <- "MODIS GPP"

#Fig. S13
p <- ggplot(sm_df, aes(x = VPD, y = .estimate, color = dataset, fill = dataset)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = .estimate - 2 * .se,
                  ymax = .estimate + 2 * .se),
              alpha = 0.15, colour = NA) +
  facet_wrap(~type, scales = "free_y") +
  scale_color_manual(values = c(
    "FLUXCOM GPP" = "#845EC2",
    "GLASS GPP"   = "#4e8397",
    "MODIS GPP"   = "#d5cabd"
  )) +
  scale_fill_manual(values = c(
    "FLUXCOM GPP" = "#845EC2",
    "GLASS GPP"   = "#4e8397",
    "MODIS GPP"   = "#d5cabd"
  )) +
  labs(
    x = "VPD",
    y = "Partial effect of s(VPD)",
    color = "Dataset",
    fill = "Dataset"
  ) +
  theme_classic(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )



