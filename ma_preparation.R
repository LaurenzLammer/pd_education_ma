# load required packages
library(readxl)
library(epitools)
library(metaHelper)
library(metafor)
library(meta)
library(dmetar)
library(metasens)
library(robvis)

# set working directory and read in the results table
setwd("C:\\Users\\Lenovo\\Documents\\Oxyfoxy")

# define variables
# define how many years of schooling are equivaent to a low vs high group comparison
year_group_diff <- 5
# define the factor for changing per level effect sizes to a two-group comparison
level_factor <- 0.5
# define the value of rho for the aggregation of multiple results in the same study
rho_val <- 0.6

# load the extracted data
data <- read_xlsx("extractions_condensed.xlsx")
# remove unnecessary parts of strings
data$Authors <- sub("[, ].*$", "", data$Authors)
# turn the Authors field for de Oliveira Souza et al. to de Oliveira Souza
data[data$Authors == "de", "Authors"] <- "de Oliveira Souza"
data[,38:43] <- lapply(data[,38:43], function(x) sub("\\[.*$", "", x))
names(data) <- sub("- Summary.*$", "RoB", names(data))
# reduce to relevant cols
df <- data[,c("Citation Name", "Authors", "Publication Date", "Key Questions", 
              "representative_pop", "smoking", "all_pd_assessed", "effect", "LL", "UL",
              "effect type", "direction", "unit", "n_controls_low_ed", 
              "n_controls_high_ed", "n_cases_low_ed", "n_cases_high_ed",
              "beta", "SE", "p-value", "Study Participation RoB", "Study Attrition RoB", 
              "Education measurement RoB", "PD measurement RoB", "Study confounding RoB", 
              "Statistical Analysis and Reporting RoB", "total sample size")]

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
  # calculate p-curve and save associated plot
  tiff(filename = paste0("pcurve_", name, ".tiff"), res = 600, units = "in", width = 12, height = 7)
  pcur <- pcurve(m.gen)
  dev.off()
  # return the produced objects
  to_return <- list(m.gen, m.gen_sub_smoking, m.gen_sub_repr, m.gen_sub_assess, lmeta, pcur)
  names(to_return) <- c("ma", "smoking", "representative", "assessed", "lmeta", "pcurve")
  return(to_return)
}

# run MA for main dataset
main_results <- metaanalyse(dataframe = main, name = "main")

# check for outliers in the main analysis
outliers <- find.outliers(main_results$ma)

# create a df without outliers
no_outliers <- main[!(main$clustering %in% outliers$out.study.random),]

# do a sensitivity analysis without outliers
no_outliers_results <- metaanalyse(dataframe = no_outliers, name = "no_outliers")

# create a df with only studies reporting effect sizes as ORs and do a sensitivity analysis
only_or <- main[main$effect.type == "OR",]
only_or_results <- metaanalyse(dataframe = only_or, name = "only_or")

# do the same excluding studies using mendelian randomization
no_mr <- main[!(main$clustering %in% c("Shi 2022", "Zhang 2022")),]
no_mr_results <- metaanalyse(dataframe = no_mr, name = "no_mr")

# do the same excluding studies with a high overall RoB
no_high_rob <- main[main$Overall != "high RoB",]
no_high_rob_results <- metaanalyse(dataframe = no_high_rob, name = "no_high_rob")

# lastly do a sensitivity analysis with only 1 Swedish register-based study
one_reg <- main[!(main$clustering %in% c("Fardell 2020", "Wirdefeldt 2005")),]
one_reg_results <- metaanalyse(dataframe = one_reg, name = "one_reg")

