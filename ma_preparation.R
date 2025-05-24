# load required packages
library(readxl)
library(epitools)
library(metaHelper)
library(metafor)
library(meta)
library(dmetar)
library(metasens)
library(robvis)
library(ggplot2)

# set working directory and read in the results table
setwd("C:\\Users\\Lenovo\\Documents\\Oxyfoxy")

# define variables
# define how many years of schooling are equivaent to a low vs high group comparison
year_group_diff <- 5
# define the factor for changing per level effect sizes to a two-group comparison
level_factor <- 0.5
# define the value of rho for the aggregation of multiple results in the same study
rho_val <- 0.6
# define the baseline risk of PD
baseline_risk <- 0.0165 # this is based on Elbaz et al. 2002

# load the extracted data
data <- read_xlsx("extractions_condensed.xlsx")
# remove unnecessary parts of strings
data$Authors <- sub("[, ].*$", "", data$Authors)
# turn the Authors field for de Oliveira Souza et al. to de Oliveira Souza
data[data$Authors == "de", "Authors"] <- "de Oliveira Souza"
data[,c(42:47)] <- lapply(data[,42:47], function(x) sub("\\[.*$", "", x))
names(data) <- sub("- Summary.*$", "RoB", names(data))
# reduce to relevant cols
df <- data[,c("Citation Name", "Authors", "Publication Date", "Key Questions", 
              "representative_pop", "smoking", "all_pd_assessed", "effect", "LL", "UL",
              "effect type", "direction", "unit", "n_controls_low_ed", 
              "n_controls_high_ed", "n_cases_low_ed", "n_cases_high_ed",
              "beta", "SE", "p-value", "Study Participation RoB", "Study Attrition RoB", 
              "Education measurement RoB", "PD measurement RoB", "Study confounding RoB", 
              "Statistical Analysis and Reporting RoB", "total sample size", "location", 
              "country income category", "PD diagnostic criteria", "comparison", 
              "Study type", "symptom_severity_assessment")]

# save  this df 
write.csv(df, "extractions_uncalculated.csv", row.names = F)
# load the df and run the script from here
df <- read.csv("extractions_uncalculated.csv", )

# align the direction of effect
# select the relevant rows
sel <- which(df$direction == 0)
# For rows where direction equals 0, update effect, LL, and UL
# swap UL and LL
df[sel, c("effect", "LL", "UL")] <- 1 / df[sel, c("effect", "UL", "LL")]
# calculate the OR with CI using the default median-unbiased estimation
# from case and control numbers in the high and low education groups
OR <- oddsratio(matrix(as.numeric(df[df$Authors == "Zou", c("n_controls_low_ed", "n_cases_low_ed", "n_controls_high_ed", "n_cases_high_ed")]), nrow = 2, byrow = TRUE))
# insert the results into df
df[df$Authors == "Zou", c("effect", "LL", "UL")] <- OR$measure["Exposed2",]
# repeat the same process for Gupta et al.
OR <- oddsratio(matrix(as.numeric(df[df$Authors == "Gupta", c("n_controls_low_ed", "n_cases_low_ed", "n_controls_high_ed", "n_cases_high_ed")]), nrow = 2, byrow = TRUE))
df[df$Authors == "Gupta", c("effect", "LL", "UL")] <- OR$measure["Exposed2",]
# calculate the OR with CI for Li et al. 2024 using the provided beta values and SEs
df[df$Authors == "Li" & df$Publication.Date == 2024, c("effect")] <- exp(df[df$Authors == "Li" & df$Publication.Date == 2024, c("beta")])
df[df$Authors == "Li" & df$Publication.Date == 2024, c("LL")] <- exp(df[df$Authors == "Li" & df$Publication.Date == 2024, c("beta")] - 
                                                                       1.96 * df[df$Authors == "Li" & df$Publication.Date == 2024, c("SE")])
df[df$Authors == "Li" & df$Publication.Date == 2024, c("UL")] <- exp(df[df$Authors == "Li" & df$Publication.Date == 2024, c("beta")] + 
                                                                       1.96 * df[df$Authors == "Li" & df$Publication.Date == 2024, c("SE")])
# calculate the OR's CI from the OR and associated p-value provided by
# Li et al. 2025 and Tan 2007
# Calculate the corresponding z-value for a two-tailed test of Li et al. 2025
z_value <- abs(qnorm(df[df$Authors == "Li" & df$Publication.Date == 2025, c("p.value")] / 2))

# Estimate the standard error on the log(Odds Ratio) scale:
se_log_or <- abs(log(df[df$Authors == "Li" & df$Publication.Date == 2025, c("effect")])) / z_value

# Calculate the lower and upper bounds of the 95% confidence interval on the log scale:
log_lower <- log(df[df$Authors == "Li" & df$Publication.Date == 2025, c("effect")]) - 1.96 * se_log_or
log_upper <- log(df[df$Authors == "Li" & df$Publication.Date == 2025, c("effect")]) + 1.96 * se_log_or

# Exponentiate to move back to the odds ratio scale:
df[df$Authors == "Li" & df$Publication.Date == 2025, c("LL")] <- exp(log_lower)
df[df$Authors == "Li" & df$Publication.Date == 2025, c("UL")] <- exp(log_upper)

# repeat the same process for Tan 2007
z_value <- abs(qnorm(df[df$Authors == "Tan", c("p.value")] / 2))

# Estimate the standard error on the log(Odds Ratio) scale:
se_log_or <- abs(log(df[df$Authors == "Tan", c("effect")])) / z_value

# Calculate the lower and upper bounds of the 95% confidence interval on the log scale:
log_lower <- log(df[df$Authors == "Tan", c("effect")]) - 1.96 * se_log_or
log_upper <- log(df[df$Authors == "Tan", c("effect")]) + 1.96 * se_log_or

# Exponentiate to move back to the odds ratio scale:
df[df$Authors == "Tan", c("LL")] <- exp(log_lower)
df[df$Authors == "Tan", c("UL")] <- exp(log_upper)

# modify the OR of the studies that used years of education as a predictor
# we assume a five year difference between the high and low group
# this number can be changed at the top of the script
rows_to_update <- which(df$unit == "year of education")
df[rows_to_update, c("effect", "LL", "UL")] <- df[rows_to_update, c("effect", "LL", "UL")]^year_group_diff
rows_to_update <- which(df$unit == "3 years of education")
df[rows_to_update, c("effect", "LL", "UL")] <- df[rows_to_update, c("effect", "LL", "UL")]^(year_group_diff/3)
# Llibre-Guerra et al. gave the effect size per level having compared five levels of education
# we decided to approximate the effect size of a two group comparison by setting the effect to OR ^ groups * 0.5
# the factor 0.5 can be changed at the top of the script
df[df$Authors == "Llibre-Guerra", c("effect", "LL", "UL")] <- 
  df[df$Authors == "Llibre-Guerra", c("effect", "LL", "UL")]^(5 * level_factor)
# determine overall RoB
# education measurement, PD measurement, study confounding & 
# analysis and reporting are considered the most important domains
# the worst score in any of these will determine the overall RoB
# if multiple of these domains have moderate RoB, it may also be deemed high RoB overall
# all these studies were deemed not to be at high RoB, though
df$Overall <- apply(df[, c("Education.measurement.RoB", "PD.measurement.RoB", 
                      "Study.confounding.RoB", "Statistical.Analysis.and.Reporting.RoB")], 1, function(x) {
  # If any column is "high", set to "high"
  if ("high RoB" %in% x) {
    "high RoB"
    # If there is more than one "moderate" and no "high", return NA
  } else if ("moderate RoB" %in% x) {
    "moderate RoB"
    # If all are "low", return "low"
  } else if (all(x == "low RoB")) {
    "low RoB"
  } else {
    NA
  }
})

# concatenate the authors and oublication year cols to create an individual clustering col
df$report <- paste(df$Authors, df$Publication.Date, sep = " ")

# order alphabetically
# sort alphabetically
df <- df[order(df$report), ]

# create a df with only references for the primary/secondary research question
df_meta <- df[df$Key.Questions == "Is more education associated with an increased risk of idiopathic Parkinson’s disease (PD)?",]
df_sec <- df[df$Key.Questions != "Is more education associated with an increased risk of idiopathic Parkinson’s disease (PD)?",]

# save these dfs
write.csv(df_meta, "outcome1.csv", row.names = F)
write.csv(df_sec, "outcome2.csv", row.names = F)

# create a traffic light plot for RoB for each outcome and save it
trafmet <- rob_traffic_light(data = df_meta[,c("report", "Study.Participation.RoB", 
                                              "Study.Attrition.RoB", "Education.measurement.RoB", 
                                              "PD.measurement.RoB", "Study.confounding.RoB", 
                                              "Statistical.Analysis.and.Reporting.RoB", "Overall")], 
                        tool = "ROB1", psize = 12)
tiff(filename = paste0("traf_met_all.tiff"), res = 600, units = "in", width = 7, height = 12)
trafmet + theme(strip.text.y.left=element_text(angle = 0))
dev.off()

trafsec <- rob_traffic_light(data = df_sec[,c("report", "Study.Participation.RoB", 
                                               "Study.Attrition.RoB", "Education.measurement.RoB", 
                                               "PD.measurement.RoB", "Study.confounding.RoB", 
                                               "Statistical.Analysis.and.Reporting.RoB", "Overall")], 
                             tool = "ROB1")
tiff(filename = paste0("traf_sec.tiff"), res = 600, units = "in", width = 7, height = 12)
trafsec + theme(strip.text.y.left=element_text(angle = 0))
dev.off()

# create an unweighted summary RoB plot for our secondary outcome
tiff(filename = paste0("sec_rob_summary.tiff"), res = 600, units = "in", width = 12, height = 7)
rob_summary(data = df_sec[,c("report", "Study.Participation.RoB", 
                             "Study.Attrition.RoB", "Education.measurement.RoB", 
                             "PD.measurement.RoB", "Study.confounding.RoB", 
                             "Statistical.Analysis.and.Reporting.RoB", "Overall")],
            weighted = F, tool = "ROB1", colour = "colourblind", overall = T)
dev.off()

# calculate log ORs
df_meta$logOR <- log(df_meta$effect)
# calculate the log OR SE
# the function goes from OR limits to the SMD SE
# to get the logOR SE it has to be additionally multiplied by pi/sqrt(3)
df_meta$SElogOR <- SE.SMD_from_OR.CI(CI_low = df_meta$LL, CI_up = df_meta$UL,
                                  sig_level = 0.05, two_tailed = TRUE) * pi/sqrt(3)
# bring data in proper format for aggregation
df_meta <- escalc(data = df_meta, yi = logOR, sei = SElogOR)
# the rho val can be modified at the top of the script
df_meta_agg <- aggregate(df_meta, cluster = report, rho = rho_val)
# calculate SEs from aggregated variances
df_meta_agg$sei <- sqrt(df_meta_agg$vi)

# remove all the references that are based on the same studies as others and 
# have a lower RoB assessment
# studies to exclude: Li et al. 2025, Li et al. 2024 and Huang et al. 2023
df_meta_agg <- df_meta_agg[!(df_meta_agg$report %in% c("Li 2025", "Li 2024", "Huang 2023")),]

# create a traffic light plot for RoB for each outcome and save it only including studies used in the MA
trafmet <- rob_traffic_light(data = df_meta_agg[,c("report", "Study.Participation.RoB", 
                                               "Study.Attrition.RoB", "Education.measurement.RoB", 
                                               "PD.measurement.RoB", "Study.confounding.RoB", 
                                               "Statistical.Analysis.and.Reporting.RoB", "Overall")], 
                             tool = "ROB1", psize = 12)
tiff(filename = paste0("traf_met_no_linked.tiff"), res = 600, units = "in", width = 7, height = 12)
trafmet + theme(strip.text.y.left=element_text(angle = 0))
dev.off()

# turn "HIC" and "LMIC" into 0s and 1s
df_meta_agg$country.income.category <- ifelse(df_meta_agg$country.income.category == "HIC", 0, 1)

# save  this df 
write.csv(df_meta_agg, "extractions_calculated.csv", row.names = F)
# load the df 
main <- read.csv("extractions_calculated.csv", )

# prepare a function to perform calculatios done in every analysis
metaanalyse <- function(dataframe, name){
  # dataframe is the respective dataframe to use
  # name gives the name of the analysis and is used to save results
  # this function performs a random-effects MA as well as three subgroup analyses
  # it produces a forest plot, a drapery plot and a funnel plot (unadjusted and adjusted) and saves them
  # it performs Rückert's limit MA
  # it performs a p-curve analysis
  # it returns the meta-analysis, the three subgroup analyses, the limit MA and the p-curve objects
  # produces a weighted summary risk of bias plot
  
  # calculate the RE MA
  m.gen <- metagen(TE = yi,
                   seTE = sei,
                   studlab = report,
                   data = dataframe,
                   sm = "OR",
                   common = F,
                   random = T,
                   method.tau = "PM",
                   method.random.ci = "HK",
                   title = name, 
                   prediction = T)
  
  # create a summary RoB plot using the weights of the RE MA
  # prepare a dataframe
  rob_df <- dataframe[,c("report", "Study.Participation.RoB", 
                         "Study.Attrition.RoB", "Education.measurement.RoB", 
                         "PD.measurement.RoB", "Study.confounding.RoB", 
                         "Statistical.Analysis.and.Reporting.RoB", "Overall")]
  rob_df$weights <- m.gen$w.random / (sum(m.gen$w.random)/100)
  tiff(filename = paste0("rob_summary_", name, ".tiff"), res = 600, units = "in", width = 12, height = 7)
  rob_summary(data = rob_df, tool = "ROB1", colour = "colourblind")
  dev.off()
  # calculate subgroup analyses
  m.gen_sub_smoking <- update(m.gen, subgroup = smoking, tau.common = T)
  m.gen_sub_repr <- update(m.gen, subgroup = representative_pop, tau.common = T)
  m.gen_sub_assess <- update(m.gen, subgroup = all_pd_assessed, tau.common = T)
  m.gen_sub_income <- update(m.gen, subgroup = country.income.category, tau.common = T)
  # produce and save a forest plot for the main and subgroup analyses
  tiff(filename = paste0("forest_", name, ".tiff"), res = 600, units = "in", width = 12, height = 7)
  meta::forest(m.gen, sortvar = TE, prediction = TRUE, print.tau2 = FALSE, 
               leftlabs = c("Author", "SE"))
  dev.off()
  tiff(filename = paste0("forest_smoking_", name, ".tiff"), res = 600, units = "in", width = 12, height = 7)
  meta::forest(m.gen_sub_smoking, sortvar = TE, prediction = TRUE, print.tau2 = FALSE, 
               leftlabs = c("Author", "SE"))
  dev.off()
  tiff(filename = paste0("forest_repr_", name, ".tiff"), res = 600, units = "in", width = 12, height = 7)
  meta::forest(m.gen_sub_repr, sortvar = TE, prediction = TRUE, print.tau2 = FALSE, 
               leftlabs = c("Author", "SE"))
  dev.off()
  tiff(filename = paste0("forest_assess_", name, ".tiff"), res = 600, units = "in", width = 12, height = 7)
  meta::forest(m.gen_sub_assess, sortvar = TE, prediction = TRUE, print.tau2 = FALSE, 
               leftlabs = c("Author", "SE"))
  dev.off()
  # produce and save a drapery plot
  tiff(filename = paste0("drapery_", name, ".tiff"), res = 600, units = "in", width = 12, height = 7)
  drapery(m.gen, 
          labels = "studlab",
          type = "pval", 
          legend = FALSE)
  dev.off()
  # produce and save a funnel plot with contour levels of alpha = 0.1, 0.05 and 0.01
  tiff(filename = paste0("funnel_", name, ".tiff"), res = 600, units = "in", width = 12, height = 7)
  meta::funnel(m.gen, studlab = T, contour.levels = c(0.9, 0.95, 0.99))
  dev.off()
  # perform Rückert's limit meta-analysis method
  lmeta <- limitmeta(m.gen)
  # produce and save an adjusted funnel plot
  tiff(filename = paste0("lmeta_funnel_", name, ".tiff"), res = 600, units = "in", width = 12, height = 7)
  funnel.limitmeta(lmeta, shrunken = T)
  dev.off()
  # calculate p-curve and save associated plots
  tiff(filename = paste0("pcurve_", name, ".tiff"), res = 600, units = "in", width = 12, height = 12)
  par(mfrow=c(2,1)) # allow both produced plots to be saved
  # create backup variable if pcurve fails with too few studies
  pcur <- try(pcurve(m.gen, effect.estimation = T, 
                     N = dataframe$total.sample.size), silent = TRUE)
  if (inherits(pcur, "try-error")) {
    pcur <- ""
  }
  dev.off()
  # prepare a dataframe to store the results of the meta-analyses
  # sorry this is some dirty code but i have no idea hoe else to extractthe results
  result_df <- data.frame(matrix(ncol = 11, nrow = 14))
  colnames(result_df) <- c("model", "n", 	"OR", "t_z_value", 
                           "p_value", "pred_interval",	"tau2", "I2", "Q", "df",	"p_value_q")
  result_df$model <- c("m.gen", "m.gen_sub_smoking", "smoking0", "smoking1", "m.gen_sub_repr", "repr0", "repr1",
                       "m.gen_sub_assess", "assess0", "assess1", "m.gen_sub_income", "income0", "income1", 
                       "lmeta")
  result_df[result_df$model %in% c("m.gen", "m.gen_sub_smoking", "m.gen_sub_repr", 
                                   "m.gen_sub_assess", "m.gen_sub_income", 
                                   "lmeta"), "n"] <- m.gen$k
                       
  result_df[result_df$model %in% c("m.gen", "m.gen_sub_smoking", "m.gen_sub_repr", 
                                   "m.gen_sub_assess", "m.gen_sub_income"), "OR"] <- paste0(round(exp(m.gen$TE.random), digits = 2), " (",
                                                       round(exp(m.gen$lower.random), digits = 2), " - ",
                                                       round(exp(m.gen$upper.random), digits = 2), ")")
  result_df[result_df$model %in% c("m.gen", "m.gen_sub_smoking", "m.gen_sub_repr", 
                                   "m.gen_sub_assess", "m.gen_sub_income"), "pred_interval"] <- paste0(round(exp(m.gen$lower.predict), digits = 2), " - ",
                                                       round(exp(m.gen$upper.predict), digits = 2))
  result_df[result_df$model %in% c("m.gen", "m.gen_sub_smoking", "m.gen_sub_repr", 
                                   "m.gen_sub_assess", "m.gen_sub_income"), "t_z_value"] <- m.gen$statistic.random
  result_df[result_df$model %in% c("m.gen", "m.gen_sub_smoking", "m.gen_sub_repr", 
                                   "m.gen_sub_assess", "m.gen_sub_income"), "p_value"] <- m.gen$pval.random
  result_df[result_df$model == "m.gen", c("tau2", "I2", "Q", "df", "p_value_q")] <-
    c(m.gen$tau2, m.gen$I2, m.gen$Q, m.gen$df.Q, m.gen$pval.Q)
  meta_list <- list(m.gen_sub_assess, m.gen_sub_income, m.gen_sub_repr, m.gen_sub_smoking)
  names(meta_list) <- c("m.gen_sub_assess", "m.gen_sub_income", "m.gen_sub_repr", "m.gen_sub_smoking")
  for (i in seq_along(meta_list)){
    obj <- meta_list[[i]]  
    result_df[result_df$model == names(meta_list)[[i]], c("tau2", "I2", "Q", "df", "p_value_q")]  <- 
      c(obj[["tau2.resid"]], obj[["I2.resid"]], obj[["Q.b.random"]], obj[["df.Q.b.random"]], obj[["pval.Q.b.random"]])
  result_df[result_df$model == paste0(sub(".*_", "", names(meta_list)[[i]]), "0"), c("n", "OR", "t_z_value", "p_value", "pred_interval", "tau2", "I2", "Q", "df", "p_value_q")] <-
      c(obj[["k.TE.w"]][[1]], paste0(round(exp(obj[["TE.random.w"]][[1]]), digits = 2), " (",
                                     round(exp(obj[["lower.random.w"]][[1]]), digits = 2), " - ",
                                     round(exp(obj[["upper.random.w"]][[1]]), digits = 2), ")"),
        obj[["statistic.random.w"]][[1]], obj[["pval.random.w"]][[1]], paste0(round(exp(obj[["lower.predict.w"]][[1]]), digits = 2), " - ", round(exp(obj[["upper.predict.w"]][[1]]), digits = 2)),
        obj[["tau2.w"]][[1]], obj[["I2.w"]][[1]], obj[["Q.w"]][[1]], obj[["df.Q.w"]][[1]], obj[["pval.Q.w"]][[1]])
  result_df[result_df$model == paste0(sub(".*_", "", names(meta_list)[[i]]), "1"), c("n", "OR", "t_z_value", "p_value", "pred_interval", "tau2", "I2", "Q", "p_value_q")] <-
    c(obj[["k.TE.w"]][[2]], paste0(round(exp(obj[["TE.random.w"]][[2]]), digits = 2), " (",
                                   round(exp(obj[["lower.random.w"]][[2]]), digits = 2), " - ",
                                   round(exp(obj[["upper.random.w"]][[2]]), digits = 2), ")"),
      obj[["statistic.random.w"]][[2]], obj[["pval.random.w"]][[2]], paste0(round(exp(obj[["lower.predict.w"]][[2]]), digits = 2), " - ", round(exp(obj[["upper.predict.w"]][[2]]), digits = 2)),
      obj[["tau2.w"]][[2]], obj[["I2.w"]][[2]], obj[["Q.w"]][[2]], obj[["pval.Q.w"]][[2]])
  }
  result_df[result_df$model == "lmeta", c("OR", "t_z_value", "p_value", "tau2", "I2", "Q", "df", "p_value_q")] <-
    c(paste0(round(exp(lmeta$TE.adjust), digits = 2), " (", round(exp(lmeta$lower.adjust), digits = 2), " - ", round(exp(lmeta$lower.adjust), digits = 2), ")"),
      lmeta$statistic.adjust, lmeta$pval.adjust, m.gen$tau2, m.gen$I2, lmeta$Q.small, 1, pchisq(lmeta$Q.small, df = 1, lower.tail = FALSE))
  write.csv(result_df, file = paste0("results_df_", name, ".csv"), row.names = F)
  # create a df to store the results of the p-curve analysis if it has been created
  if (!is.character(pcur)) {
    pcurve_df <- data.frame(matrix(nrow = 2, ncol = 12))
    colnames(pcurve_df) <- c("n", "n_sig", "n_sub_0.025", "test", "pBinomial",  
                             "zFull",	"pFull", "zHalf",	"pHalf", "Power_estimate",	"Evidential_value", "Effect_estimate")
    pcurve_df[, c("n", "n_sig", "n_sub_0.025")] <- rep(c(pcur$kInput, pcur$kAnalyzed, pcur$kp0.25), each = 2)
    pcurve_df$test <- c("right-skewness", "flatness")
    pcurve_df[,c("pBinomial", "zFull",	"pFull", "zHalf",	"pHalf")] <- 
      pcur$pcurveResults
    pcurve_df$Power_estimate <- paste0(pcur$Power[[1]], " (", pcur$Power[[2]], " - ", pcur$Power[[3]], ")")
    pcurve_df$Evidential_value <- c(pcur$EvidencePresent, pcur$EvidenceAbsent)
    pcurve_df$Effect_estimate <- pcur$dEstimate
    # save the df as a csv file
    write.csv(pcurve_df, file = paste0("results_pcurve_", name, ".csv"), row.names = F)
  }
  # return the produced objects
  to_return <- list(m.gen, m.gen_sub_smoking, m.gen_sub_repr, m.gen_sub_assess, m.gen_sub_income, lmeta, pcur)
  names(to_return) <- c("ma", "smoking", "representative", "assessed", "income", "lmeta", "pcurve")
  return(to_return)
}

# run MA for main dataset
main_results <- metaanalyse(dataframe = main, name = "main")

# check for outliers in the main analysis
outliers <- find.outliers(main_results$ma)

# create a df without outliers
no_outliers <- main[!(main$report %in% outliers$out.study.random),]

# do a sensitivity analysis without outliers
no_outliers_results <- metaanalyse(dataframe = no_outliers, name = "no_outliers")

# create a df with only studies reporting effect sizes as ORs and do a sensitivity analysis
only_or <- main[main$effect.type == "OR",]
only_or_results <- metaanalyse(dataframe = only_or, name = "only_or")

# do the same excluding studies using mendelian randomization
no_mr <- main[!(main$report %in% c("Shi 2022", "Zhang 2022")),]
no_mr_results <- metaanalyse(dataframe = no_mr, name = "no_mr")

# do the same excluding studies with a high overall RoB
no_high_rob <- main[main$Overall != "high RoB",]
no_high_rob_results <- metaanalyse(dataframe = no_high_rob, name = "no_high_rob")

# lastly do a sensitivity analysis with only 1 Swedish register-based study
one_reg <- main[!(main$report %in% c("Fardell 2020", "Wirdefeldt 2005")),]
one_reg_results <- metaanalyse(dataframe = one_reg, name = "one_reg")

# calculate the risk of PD based on the OR from the meta-analysis and the baseline risk 
absolute_risk <- ((baseline_risk * exp(main_results$ma$TE.random)) / 
  (1 - baseline_risk + exp(main_results$ma$TE.random) * baseline_risk))*100
ll_absolute_risk <- ((baseline_risk * exp(main_results$ma$lower.random)) / 
                       (1 - baseline_risk + exp(main_results$ma$lower.random) * baseline_risk))*100
ul_absolute_risk <- ((baseline_risk * exp(main_results$ma$upper.random)) / 
                       (1 - baseline_risk + exp(main_results$ma$upper.random) * baseline_risk))*100

risk_diff <- absolute_risk - baseline_risk*100
ll_risk_diff <- ll_absolute_risk - baseline_risk*100
ul_risk_diff <- ul_absolute_risk - baseline_risk*100

risk_df <- data.frame(
  level = c("mean", "upper", "lower"),
  absolute = c(absolute_risk, ul_absolute_risk, ll_absolute_risk),
  difference = c(risk_diff, ul_risk_diff, ll_risk_diff),
  OR = c(exp(main_results$ma$TE.random), exp(main_results$ma$upper.random),
         exp(main_results$ma$lower.random))
)
write.csv(risk_df, "risk_differences.csv", row.names = F)
