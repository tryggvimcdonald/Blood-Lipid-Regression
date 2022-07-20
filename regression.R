library(stats)
library(dplyr)
library(ggplot2)
library(ggpubr)

#Load in bd_pheno dataset with summed scores already added in
load("C:/Users/pierc/OneDrive/Desktop/bd_pheno_scores.RData")

#Select important columns
score_phenos <- select(bd_pheno, Sex, Fish_oil_baseline, Tot_Chol, LDL, HDL, TAGs, df_TC, df_TAG, df_HDL, df_LDL)

#Standardize blood lipid phenotypes
score_phenos <- mutate(score_phenos, Stan_TC = scale(Tot_Chol), Stan_LDL = scale(LDL), Stan_HDL = scale(HDL), Stan_TAG = scale(TAGs))

#Male and Female separation
score_phenos <- mutate(score_phenos, Male_TC = case_when(Sex == 1 ~ Tot_Chol), Male_LDL = case_when(Sex == 1 ~ LDL), Male_HDL = case_when(Sex == 1 ~ HDL), Male_TAG = case_when(Sex == 1 ~ TAGs), Female_TC = case_when(Sex == 0 ~ Tot_Chol), Female_LDL = case_when(Sex == 0 ~ LDL), Female_HDL = case_when(Sex == 0 ~ HDL), Female_TAG = case_when(Sex == 0 ~ TAGs))

# Fish Oil Separation
score_phenos <- mutate(score_phenos, Fish_TC = case_when(Fish_oil_baseline == 1 ~ Tot_Chol), Fish_LDL = case_when(Fish_oil_baseline == 1 ~ LDL), Fish_HDL = case_when(Fish_oil_baseline == 1 ~ HDL), Fish_TAG = case_when(Fish_oil_baseline == 1 ~ TAGs), No_Fish_TC = case_when(Fish_oil_baseline == 0 ~ Tot_Chol), No_Fish_LDL = case_when(Fish_oil_baseline == 0 ~ LDL), No_Fish_HDL = case_when(Fish_oil_baseline == 0 ~ HDL), No_Fish_TAG = case_when(Fish_oil_baseline == 0 ~ TAGs))

#Take sample for plotting
plot_samples <- slice_sample(score_phenos, n = 5000)
plot_samples$Sex <- as.numeric(plot_samples$Sex)
plot_samples <- mutate(plot_samples, Sex = factor(pull(plot_samples, Sex),
                                                  levels = c(0,1),
                                                  labels = c("Female","Male")))
plot_samples$Fish_oil_baseline <- as.factor(plot_samples$Fish_oil_baseline)


#Unstandardized Plots
TC_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_TC), !is.na(Tot_Chol)), x = "df_TC", y = "Tot_Chol", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "Total Cholesterol PRS", ylab = "Total Cholesterol Level", alpha = 0.3)
HDL_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_HDL), !is.na(HDL)), x = "df_HDL", y = "HDL", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "HDL PRS", ylab = "HDL Level", alpha = 0.3)
LDL_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_LDL), !is.na(LDL)), x = "df_LDL", y = "LDL", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "LDL PRS", ylab = "LDL Level", alpha = 0.3)
TAG_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_TAG), !is.na(TAGs)), x = "df_TAG", y = "TAGs", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "Triglycerides PRS", ylab = "Triglycerides Level", alpha = 0.3)

scatters <- ggarrange(TC_scatter, HDL_scatter, LDL_scatter, TAG_scatter, ncol = 2, nrow = 2)
scatters <- annotate_figure(scatters, top = "Phenotype by PRS Scatterplots")

ggexport(scatters, filename = paste0(getwd(),"/Blood-Lipid-Phe-Score-Scatters.png"), width = 1600, height = 900)


# Standardized plots
Stan_TC_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_TC), !is.na(Stan_TC)), x = "df_TC", y = "Stan_TC", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "Total Cholesterol PRS", ylab = "Standardized Total Cholesterol Level", alpha = 0.3)
Stan_HDL_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_HDL), !is.na(Stan_HDL)), x = "df_HDL", y = "Stan_HDL", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "HDL PRS", ylab = "Standardized HDL Level", alpha = 0.3)
Stan_LDL_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_LDL), !is.na(Stan_LDL)), x = "df_LDL", y = "Stan_LDL", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "LDL PRS", ylab = "Standardized LDL Level", alpha = 0.3)
Stan_TAG_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_TAG), !is.na(Stan_TAG)), x = "df_TAG", y = "Stan_TAG", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "Triglycerides PRS", ylab = "Standardized Triglycerides Level", alpha = 0.3)

stan_scatters <- ggarrange(Stan_TC_scatter, Stan_HDL_scatter, Stan_LDL_scatter, Stan_TAG_scatter, ncol = 2, nrow = 2)
stan_scatters <- annotate_figure(stan_scatters, top = "Standardized Phenotype by PRS Scatterplots")

ggexport(stan_scatters, filename = paste0(getwd(),"/Standardized-Blood-Lipid-Phe-Score-Scatters.png"), width = 1600, height = 900)

#Male and Female Plots
Male_TC_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_TC), !is.na(Male_TC)), x = "df_TC", y = "Male_TC", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "Total Cholesterol PRS", ylab = "Total Cholesterol Level", color = "Sex", palette = c("#55c1e6"), alpha = 0.3)
Male_HDL_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_HDL), !is.na(Male_HDL)), x = "df_HDL", y = "Male_HDL", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "HDL PRS", ylab = "HDL Level", color = "Sex", palette = c("#55c1e6"), alpha = 0.3)
Male_LDL_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_LDL), !is.na(Male_LDL)), x = "df_LDL", y = "Male_LDL", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "LDL PRS", ylab = "LDL Level", color = "Sex", palette = c("#55c1e6"), alpha = 0.3)
Male_TAG_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_TAG), !is.na(Male_TAG)), x = "df_TAG", y = "Male_TAG", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "Triglycerides PRS", ylab = "Triglycerides Level", color = "Sex", palette = c("#55c1e6"), alpha = 0.3)

Female_TC_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_TC), !is.na(Female_TC)), x = "df_TC", y = "Female_TC", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "Total Cholesterol PRS", ylab = "Total Cholesterol Level", color = "Sex", palette = c("#e65594"), alpha = 0.3) + theme(legend.title = element_text(color = "white", size = 1))
Female_HDL_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_HDL), !is.na(Female_HDL)), x = "df_HDL", y = "Female_HDL", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "HDL PRS", ylab = "HDL Level", color = "Sex", palette = c("#e65594"), alpha = 0.3)
Female_LDL_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_LDL), !is.na(Female_LDL)), x = "df_LDL", y = "Female_LDL", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "LDL PRS", ylab = "LDL Level", color = "Sex", palette = c("#e65594"), alpha = 0.3)
Female_TAG_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_TAG), !is.na(Female_TAG)), x = "df_TAG", y = "Female_TAG", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "Triglycerides PRS", ylab = "Triglycerides Level", color = "Sex", palette = c("#e65594"), alpha = 0.3)

Male_legend <- get_legend(Male_TC_scatter)
Female_legend <- get_legend(Female_TC_scatter)
Legends <- ggarrange(Male_legend, NULL, Female_legend, ncol = 3, widths = c(1, -0.95, 1))

sex_scatters <- ggarrange(Male_TC_scatter, Male_HDL_scatter, Male_LDL_scatter, Male_TAG_scatter, Female_TC_scatter, Female_HDL_scatter, Female_LDL_scatter, Female_TAG_scatter, ncol = 4, nrow = 2, legend = "bottom", common.legend = TRUE, legend.grob = element_grob(element))
sex_scatters <- annotate_figure(sex_scatters, top = "Phenotype by PRS Scatterplots, Stratified by Sex")

sex_scatter_legends <- ggarrange(sex_scatters, Legends, heights = c(0.95, 0.05), nrow = 2)

ggexport(sex_scatter_legends, filename = paste0(getwd(),"/Sex-Blood-Lipid-Phe-Score-Scatters.png"), width = 1600, height = 900)


# Fish Oil Plots
Fish_TC_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_TC), !is.na(Fish_TC)), x = "df_TC", y = "Fish_TC", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "Total Cholesterol PRS", ylab = "Total Cholesterol Level", color = "Fish_oil_baseline", palette = "#35b544", alpha = 0.3) + scale_color_manual(name = "Fish Oil Status", values = "#35b544", labels = "Yes")
Fish_HDL_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_HDL), !is.na(Fish_HDL)), x = "df_HDL", y = "Fish_HDL", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "HDL PRS", ylab = "HDL Level", color = "Fish_oil_baseline", palette = c("#35b544"), alpha = 0.3)
Fish_LDL_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_LDL), !is.na(Fish_LDL)), x = "df_LDL", y = "Fish_LDL", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "LDL PRS", ylab = "LDL Level", color = "Fish_oil_baseline", palette = c("#35b544"), alpha = 0.3)
Fish_TAG_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_TAG), !is.na(Fish_TAG)), x = "df_TAG", y = "Fish_TAG", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "Triglycerides PRS", ylab = "Triglycerides Level", color = "Fish_oil_baseline", palette = c("#35b544"), alpha = 0.3)

No_Fish_TC_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_TC), !is.na(No_Fish_TC)), x = "df_TC", y = "No_Fish_TC", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "Total Cholesterol PRS", ylab = "Total Cholesterol Level", color = "Fish_oil_baseline", palette = c("#FC7174"), alpha = 0.3) + theme(legend.title = element_text(color = "white", size = 0)) + scale_color_discrete(name = "Fish Oil Status", labels = "No")
No_Fish_HDL_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_HDL), !is.na(No_Fish_HDL)), x = "df_HDL", y = "No_Fish_HDL", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "HDL PRS", ylab = "HDL Level", color = "Fish_oil_baseline", palette = c("#FC7174"), alpha = 0.3)
No_Fish_LDL_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_LDL), !is.na(No_Fish_LDL)), x = "df_LDL", y = "No_Fish_LDL", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "LDL PRS", ylab = "LDL Level", color = "Fish_oil_baseline", palette = c("#FC7174"), alpha = 0.3)
No_Fish_TAG_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_TAG), !is.na(No_Fish_TAG)), x = "df_TAG", y = "No_Fish_TAG", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "Triglycerides PRS", ylab = "Triglycerides Level", color = "Fish_oil_baseline", palette = c("#FC7174"), alpha = 0.3)

Fish_legend <- get_legend(Fish_TC_scatter)
No_Fish_legend <- get_legend(No_Fish_TC_scatter)
Legends <- ggarrange(Fish_legend, NULL, No_Fish_legend, ncol = 3, widths = c(1, -0.94, 1))

fish_scatters <- ggarrange(Fish_TC_scatter, Fish_HDL_scatter, Fish_LDL_scatter, Fish_TAG_scatter, No_Fish_TC_scatter, No_Fish_HDL_scatter, No_Fish_LDL_scatter, No_Fish_TAG_scatter, nrow = 2, ncol = 4, legend = "none")
fish_scatters <- annotate_figure(fish_scatters, top = "Phenotype by PRS Scatterplots, Stratified by Fish Oil Supplementation Status")

fish_scatter_legends <- ggarrange(fish_scatters, Legends, heights = c(0.95,0.05), nrow = 2)

ggexport(fish_scatter_legends, filename = paste0(getwd(),"/Fish-Blood-Lipid-Phe-Score-Scatters.png"), width = 1600, height = 900)

# Sex-Fish Combo
Male_Fish_TC_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_TC), !is.na(Fish_TC), Sex == "Male"), x = "df_TC", y = "Fish_TC", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "Total Cholesterol PRS", ylab = "Total Cholesterol Level", color = c("chartreuse3"), alpha = 0.4, shape = "Sex") + scale_shape_manual(values = 15, guide = guide_legend(override.aes = list(color = "black")))
Male_Fish_HDL_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_HDL), !is.na(Fish_HDL), Sex == "Male"), x = "df_HDL", y = "Fish_HDL", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "HDL PRS", ylab = "HDL Level", color = "Fish_oil_baseline", palette = c("chartreuse3"), alpha = 0.4, shape = 15) + scale_color_manual(name = "Fish Oil Status", values = "chartreuse3", labels = "Yes")
Male_Fish_LDL_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_LDL), !is.na(Fish_LDL), Sex == "Male"), x = "df_LDL", y = "Fish_LDL", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "LDL PRS", ylab = "LDL Level", color = "Fish_oil_baseline", palette = c("chartreuse3"), alpha = 0.4, shape = "Sex") + scale_shape_manual(values = 15)
Male_Fish_TAG_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_TAG), !is.na(Fish_TAG), Sex == "Male"), x = "df_TAG", y = "Fish_TAG", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "Triglycerides PRS", ylab = "Triglycerides Level", color = "Fish_oil_baseline", palette = c("chartreuse3"), alpha = 0.4, shape = "Sex") + scale_shape_manual(values = 15)

Female_Fish_TC_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_TC), !is.na(Fish_TC), Sex == "Female"), x = "df_TC", y = "Fish_TC", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "Total Cholesterol PRS", ylab = "Total Cholesterol Level", color = "Fish_oil_baseline", palette = c("chartreuse3"), alpha = 0.4, shape = "Sex") + scale_shape_manual(values = 18)
Female_Fish_HDL_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_HDL), !is.na(Fish_HDL), Sex == "Female"), x = "df_HDL", y = "Fish_HDL", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "HDL PRS", ylab = "HDL Level", color = "Fish_oil_baseline", palette = c("chartreuse3"), alpha = 0.4, shape = "Sex") + scale_shape_manual(values = 18)
Female_Fish_LDL_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_LDL), !is.na(Fish_LDL), Sex == "Female"), x = "df_LDL", y = "Fish_LDL", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "LDL PRS", ylab = "LDL Level", color = "Fish_oil_baseline", palette = c("chartreuse3"), alpha = 0.4, shape = "Sex") + scale_shape_manual(values = 18)
Female_Fish_TAG_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_TAG), !is.na(Fish_TAG), Sex == "Female"), x = "df_TAG", y = "Fish_TAG", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "Triglycerides PRS", ylab = "Triglycerides Level", color = "Fish_oil_baseline", palette = c("chartreuse3"), alpha = 0.4, shape = "Sex") + scale_shape_manual(values = 18)

Male_No_Fish_TC_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_TC), !is.na(No_Fish_TC), Sex == "Male"), x = "df_TC", y = "No_Fish_TC", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "Total Cholesterol PRS", ylab = "Total Cholesterol Level", color = "Fish_oil_baseline", palette = c("#FC7174"), alpha = 0.4, shape = "Sex") + scale_shape_manual(values = 15)
Male_No_Fish_HDL_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_HDL), !is.na(No_Fish_HDL), Sex == "Male"), x = "df_HDL", y = "No_Fish_HDL", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "HDL PRS", ylab = "HDL Level", color = "Fish_oil_baseline", palette = c("#FC7174"), alpha = 0.4, shape = 15) + theme(legend.title = element_text(color = "white", size = 0)) + scale_color_manual(name = "Fish Oil Status", values = "#FC7174", labels = "No")
Male_No_Fish_LDL_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_LDL), !is.na(No_Fish_LDL), Sex == "Male"), x = "df_LDL", y = "No_Fish_LDL", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "LDL PRS", ylab = "LDL Level", color = "Fish_oil_baseline", palette = c("#FC7174"), alpha = 0.4, shape = "Sex") + scale_shape_manual(values = 15)
Male_No_Fish_TAG_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_TAG), !is.na(No_Fish_TAG), Sex == "Male"), x = "df_TAG", y = "No_Fish_TAG", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "Triglycerides PRS", ylab = "Triglycerides Level", color = "Fish_oil_baseline", palette = c("#FC7174"), alpha = 0.4, shape = "Sex") + scale_shape_manual(values = 15)

Female_No_Fish_TC_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_TC), !is.na(No_Fish_TC), Sex == "Female"), x = "df_TC", y = "No_Fish_TC", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "Total Cholesterol PRS", ylab = "Total Cholesterol Level", color = c("#FC7174"), alpha = 0.4, shape = "Sex") + scale_shape_manual(values = 18, guide = guide_legend(override.aes = list(color = "black"))) + theme(legend.title = element_text(color = "white", size = 0))
Female_No_Fish_HDL_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_HDL), !is.na(No_Fish_HDL), Sex == "Female"), x = "df_HDL", y = "No_Fish_HDL", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "HDL PRS", ylab = "HDL Level", color = "Fish_oil_baseline", palette = c("#FC7174"), alpha = 0.4, shape = "Sex") + scale_shape_manual(values = 18)
Female_No_Fish_LDL_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_LDL), !is.na(No_Fish_LDL), Sex == "Female"), x = "df_LDL", y = "No_Fish_LDL", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "LDL PRS", ylab = "LDL Level", color = "Fish_oil_baseline", palette = c("#FC7174"), alpha = 0.4, shape = "Sex") + scale_shape_manual(values = 18)
Female_No_Fish_TAG_scatter <- ggscatter(data = filter(plot_samples, !is.na(df_TAG), !is.na(No_Fish_TAG), Sex == "Female"), x = "df_TAG", y = "No_Fish_TAG", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", add.params = list(color = "red"), xlab = "Triglycerides PRS", ylab = "Triglycerides Level", color = "Fish_oil_baseline", palette = c("#FC7174"), alpha = 0.4, shape = "Sex") + scale_shape_manual(values = 18)

sex_fish_scatters <- ggarrange(Male_Fish_TC_scatter, Male_Fish_HDL_scatter, Male_Fish_LDL_scatter, Male_Fish_TAG_scatter, Female_Fish_TC_scatter, Female_Fish_HDL_scatter, Female_Fish_LDL_scatter, Female_Fish_TAG_scatter, Male_No_Fish_TC_scatter, Male_No_Fish_HDL_scatter, Male_No_Fish_LDL_scatter, Male_No_Fish_TAG_scatter, Female_No_Fish_TC_scatter, Female_No_Fish_HDL_scatter, Female_No_Fish_LDL_scatter, Female_No_Fish_TAG_scatter, nrow = 4, ncol = 4, legend = "none")
sex_fish_scatters <- annotate_figure(sex_fish_scatters, top = "Phenotype by PRS Scatterplots, Stratified by Fish Oil Supplementation Status and Sex")


Male_legend <- get_legend(Male_Fish_TC_scatter)
Female_legend <- get_legend(Female_No_Fish_TC_scatter)
Fish_legend <- get_legend(Male_Fish_HDL_scatter)
No_Fish_legend <- get_legend(Male_No_Fish_HDL_scatter)

Legends <- ggarrange(Male_legend, NULL, Female_legend, NULL, Fish_legend, NULL, No_Fish_legend, ncol = 7, widths = c(1,-0.9,1,-0.5,1,-0.9,1))

sex_fish_scatters_legends <- ggarrange(sex_fish_scatters, Legends, heights = c(0.95,0.05), nrow = 2)

ggexport(sex_fish_scatters_legends, filename = paste0(getwd(),"/Sex-Fish-Blood-Lipid-Phe-Score-Scatters.png"), width = 1600, height = 900)


# Normality Tests
norm_tests <- shapiro_test(plot_samples, Tot_Chol, LDL, HDL, TAGs, df_TC, df_LDL, df_HDL, df_TAG)
norm_tests <- add_significance(norm_tests, p.col = "p", output.col = "signif")
norm_tests

#Correlation Tests
TC_cor <- cor_test(score_phenos, Tot_Chol, df_TC, method = "spearman")
TC_cor <- select(TC_cor,!method)

LDL_cor <- cor_test(score_phenos, LDL, df_LDL, method = "spearman")
LDL_cor <- select(LDL_cor,!method)

HDL_cor <- cor_test(score_phenos, HDL, df_HDL, method = "spearman")
HDL_cor <- select(HDL_cor,!method)

TAG_cor <- cor_test(score_phenos, TAGs, df_TAG, method = "spearman")
TAG_cor <- select(TAG_cor,!method)

cor_tests <- bind_rows(TC_cor, LDL_cor, HDL_cor, TAG_cor)

cat(format(cor_tests)[c(-1L,-3L)], sep = "\n")




# Level 1 Linear Models
TC_1_lm <- lm(Tot_Chol ~ df_TC, data = bd_pheno, na.action = "na.omit")
HDL_1_lm <- lm(HDL ~ df_HDL, data = bd_pheno, na.action = "na.omit")
LDL_1_lm <- lm(LDL ~ df_LDL, data = bd_pheno, na.action = "na.omit")
TAG_1_lm <- lm(TAGs ~ df_TAG, data = bd_pheno, na.action = "na.omit")

TC_1_lm <- tidy(TC_1_lm)
HDL_1_lm <- tidy(HDL_1_lm)
LDL_1_lm <- tidy(LDL_1_lm)
TAG_1_lm <- tidy(TAG_1_lm)

lm1 <- bind_rows(TC_1_lm, HDL_1_lm, LDL_1_lm, TAG_1_lm)





# Level 2 Linear Models
TC_2_lm <- lm(Tot_Chol ~ df_TC + Sex + Age + I(Age^2) + Array + 
                   Assessment_centres_10003 + Assessment_centres_11001 + Assessment_centres_11002 + 
                   Assessment_centres_11003 + Assessment_centres_11004 + Assessment_centres_11005 + 
                   Assessment_centres_11006 + Assessment_centres_11007 + Assessment_centres_11008 + 
                   Assessment_centres_11009 + Assessment_centres_11010 + Assessment_centres_11011 + 
                   Assessment_centres_11011 + Assessment_centres_11012 + Assessment_centres_11013 + 
                   Assessment_centres_11014 + Assessment_centres_11016 + Assessment_centres_11017 + 
                   Assessment_centres_11018 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + 
                   PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + 
                   PCA19 + PCA20, data = bd_pheno, na.action = "na.omit")

HDL_2_lm <- lm(HDL ~ df_HDL + Sex + Age + I(Age^2) + Array + 
                   Assessment_centres_10003 + Assessment_centres_11001 + Assessment_centres_11002 + 
                   Assessment_centres_11003 + Assessment_centres_11004 + Assessment_centres_11005 + 
                   Assessment_centres_11006 + Assessment_centres_11007 + Assessment_centres_11008 + 
                   Assessment_centres_11009 + Assessment_centres_11010 + Assessment_centres_11011 + 
                   Assessment_centres_11011 + Assessment_centres_11012 + Assessment_centres_11013 + 
                   Assessment_centres_11014 + Assessment_centres_11016 + Assessment_centres_11017 + 
                   Assessment_centres_11018 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + 
                   PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + 
                   PCA19 + PCA20, data = bd_pheno, na.action = "na.omit")

LDL_2_lm <- lm(LDL ~ df_LDL + Sex + Age + I(Age^2) + Array + 
                   Assessment_centres_10003 + Assessment_centres_11001 + Assessment_centres_11002 + 
                   Assessment_centres_11003 + Assessment_centres_11004 + Assessment_centres_11005 + 
                   Assessment_centres_11006 + Assessment_centres_11007 + Assessment_centres_11008 + 
                   Assessment_centres_11009 + Assessment_centres_11010 + Assessment_centres_11011 + 
                   Assessment_centres_11011 + Assessment_centres_11012 + Assessment_centres_11013 + 
                   Assessment_centres_11014 + Assessment_centres_11016 + Assessment_centres_11017 + 
                   Assessment_centres_11018 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + 
                   PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + 
                   PCA19 + PCA20, data = bd_pheno, na.action = "na.omit")

TAG_2_lm <- lm(TAGs ~ df_TAG + Sex + Age + I(Age^2) + Array + 
                   Assessment_centres_10003 + Assessment_centres_11001 + Assessment_centres_11002 + 
                   Assessment_centres_11003 + Assessment_centres_11004 + Assessment_centres_11005 + 
                   Assessment_centres_11006 + Assessment_centres_11007 + Assessment_centres_11008 + 
                   Assessment_centres_11009 + Assessment_centres_11010 + Assessment_centres_11011 + 
                   Assessment_centres_11011 + Assessment_centres_11012 + Assessment_centres_11013 + 
                   Assessment_centres_11014 + Assessment_centres_11016 + Assessment_centres_11017 + 
                   Assessment_centres_11018 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + 
                   PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + 
                   PCA19 + PCA20, data = bd_pheno, na.action = "na.omit")

TC_2_lm <- tidy(TC_2_lm)
HDL_2_lm <- tidy(HDL_2_lm)
LDL_2_lm <- tidy(LDL_2_lm)
TAG_2_lm <- tidy(TAG_2_lm)


# Level 3 Linear Models
TC_3_lm <- lm(Tot_Chol ~ df_TC + Sex + Age + I(Age^2) + Array + 
                Assessment_centres_10003 + Assessment_centres_11001 + Assessment_centres_11002 + 
                Assessment_centres_11003 + Assessment_centres_11004 + Assessment_centres_11005 + 
                Assessment_centres_11006 + Assessment_centres_11007 + Assessment_centres_11008 + 
                Assessment_centres_11009 + Assessment_centres_11010 + Assessment_centres_11011 + 
                Assessment_centres_11011 + Assessment_centres_11012 + Assessment_centres_11013 + 
                Assessment_centres_11014 + Assessment_centres_11016 + Assessment_centres_11017 + 
                Assessment_centres_11018 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + 
                PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + 
                PCA19 + PCA20 + BMI + Townsend, data = bd_pheno, na.action = "na.omit")

HDL_3_lm <- lm(HDL ~ df_HDL + Sex + Age + I(Age^2) + Array + 
                     Assessment_centres_10003 + Assessment_centres_11001 + Assessment_centres_11002 + 
                     Assessment_centres_11003 + Assessment_centres_11004 + Assessment_centres_11005 + 
                     Assessment_centres_11006 + Assessment_centres_11007 + Assessment_centres_11008 + 
                     Assessment_centres_11009 + Assessment_centres_11010 + Assessment_centres_11011 + 
                     Assessment_centres_11011 + Assessment_centres_11012 + Assessment_centres_11013 + 
                     Assessment_centres_11014 + Assessment_centres_11016 + Assessment_centres_11017 + 
                     Assessment_centres_11018 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + 
                     PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + 
                     PCA19 + PCA20 + BMI + Townsend, data = bd_pheno, na.action = "na.omit")

LDL_3_lm <- lm(LDL ~ df_LDL + Sex + Age + I(Age^2) + Array + 
                     Assessment_centres_10003 + Assessment_centres_11001 + Assessment_centres_11002 + 
                     Assessment_centres_11003 + Assessment_centres_11004 + Assessment_centres_11005 + 
                     Assessment_centres_11006 + Assessment_centres_11007 + Assessment_centres_11008 + 
                     Assessment_centres_11009 + Assessment_centres_11010 + Assessment_centres_11011 + 
                     Assessment_centres_11011 + Assessment_centres_11012 + Assessment_centres_11013 + 
                     Assessment_centres_11014 + Assessment_centres_11016 + Assessment_centres_11017 + 
                     Assessment_centres_11018 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + 
                     PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + 
                     PCA19 + PCA20 + BMI + Townsend, data = bd_pheno, na.action = "na.omit")

TAG_3_lm <- lm(TAGs ~ df_TAG + Sex + Age + I(Age^2) + Array + 
                     Assessment_centres_10003 + Assessment_centres_11001 + Assessment_centres_11002 + 
                     Assessment_centres_11003 + Assessment_centres_11004 + Assessment_centres_11005 + 
                     Assessment_centres_11006 + Assessment_centres_11007 + Assessment_centres_11008 + 
                     Assessment_centres_11009 + Assessment_centres_11010 + Assessment_centres_11011 + 
                     Assessment_centres_11011 + Assessment_centres_11012 + Assessment_centres_11013 + 
                     Assessment_centres_11014 + Assessment_centres_11016 + Assessment_centres_11017 + 
                     Assessment_centres_11018 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + 
                     PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + 
                     PCA19 + PCA20 + BMI + Townsend, data = bd_pheno, na.action = "na.omit")

TC_3_lm <- tidy(TC_3_lm)
HDL_3_lm <- tidy(HDL_3_lm)
LDL_3_lm <- tidy(LDL_3_lm)
TAG_3_lm <- tidy(TAG_3_lm)


# Sex Interaction
TC_Male_lm <- lm(Tot_Chol ~ df_TC + Sex + Age + I(Age^2) + Array + 
                     Assessment_centres_10003 + Assessment_centres_11001 + Assessment_centres_11002 + 
                     Assessment_centres_11003 + Assessment_centres_11004 + Assessment_centres_11005 + 
                     Assessment_centres_11006 + Assessment_centres_11007 + Assessment_centres_11008 + 
                     Assessment_centres_11009 + Assessment_centres_11010 + Assessment_centres_11011 + 
                     Assessment_centres_11011 + Assessment_centres_11012 + Assessment_centres_11013 + 
                     Assessment_centres_11014 + Assessment_centres_11016 + Assessment_centres_11017 + 
                     Assessment_centres_11018 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + 
                     PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + 
                     PCA19 + PCA20 + BMI + Townsend, data = filter(bd_pheno, Fish_oil_baseline == 1), na.action = "na.omit")

HDL_Male_lm <- lm(HDL ~ df_HDL + Sex + Age + I(Age^2) + Array + 
                      Assessment_centres_10003 + Assessment_centres_11001 + Assessment_centres_11002 + 
                      Assessment_centres_11003 + Assessment_centres_11004 + Assessment_centres_11005 + 
                      Assessment_centres_11006 + Assessment_centres_11007 + Assessment_centres_11008 + 
                      Assessment_centres_11009 + Assessment_centres_11010 + Assessment_centres_11011 + 
                      Assessment_centres_11011 + Assessment_centres_11012 + Assessment_centres_11013 + 
                      Assessment_centres_11014 + Assessment_centres_11016 + Assessment_centres_11017 + 
                      Assessment_centres_11018 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + 
                      PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + 
                      PCA19 + PCA20 + BMI + Townsend, data = filter(bd_pheno, Fish_oil_baseline == 1), na.action = "na.omit")

LDL_Male_lm <- lm(LDL ~ df_LDL + Sex + Age + I(Age^2) + Array + 
                      Assessment_centres_10003 + Assessment_centres_11001 + Assessment_centres_11002 + 
                      Assessment_centres_11003 + Assessment_centres_11004 + Assessment_centres_11005 + 
                      Assessment_centres_11006 + Assessment_centres_11007 + Assessment_centres_11008 + 
                      Assessment_centres_11009 + Assessment_centres_11010 + Assessment_centres_11011 + 
                      Assessment_centres_11011 + Assessment_centres_11012 + Assessment_centres_11013 + 
                      Assessment_centres_11014 + Assessment_centres_11016 + Assessment_centres_11017 + 
                      Assessment_centres_11018 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + 
                      PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + 
                      PCA19 + PCA20 + BMI + Townsend, data = filter(bd_pheno, Fish_oil_baseline == 1), na.action = "na.omit")

TAG_Male_lm <- lm(TAGs ~ df_TAG + Sex + Age + I(Age^2) + Array + 
                      Assessment_centres_10003 + Assessment_centres_11001 + Assessment_centres_11002 + 
                      Assessment_centres_11003 + Assessment_centres_11004 + Assessment_centres_11005 + 
                      Assessment_centres_11006 + Assessment_centres_11007 + Assessment_centres_11008 + 
                      Assessment_centres_11009 + Assessment_centres_11010 + Assessment_centres_11011 + 
                      Assessment_centres_11011 + Assessment_centres_11012 + Assessment_centres_11013 + 
                      Assessment_centres_11014 + Assessment_centres_11016 + Assessment_centres_11017 + 
                      Assessment_centres_11018 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + 
                      PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + 
                      PCA19 + PCA20 + BMI + Townsend, data = filter(bd_pheno, Fish_oil_baseline == 1), na.action = "na.omit")


TC_Female_lm <- lm(Tot_Chol ~ df_TC + Sex + Age + I(Age^2) + Array + 
                        Assessment_centres_10003 + Assessment_centres_11001 + Assessment_centres_11002 + 
                        Assessment_centres_11003 + Assessment_centres_11004 + Assessment_centres_11005 + 
                        Assessment_centres_11006 + Assessment_centres_11007 + Assessment_centres_11008 + 
                        Assessment_centres_11009 + Assessment_centres_11010 + Assessment_centres_11011 + 
                        Assessment_centres_11011 + Assessment_centres_11012 + Assessment_centres_11013 + 
                        Assessment_centres_11014 + Assessment_centres_11016 + Assessment_centres_11017 + 
                        Assessment_centres_11018 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + 
                        PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + 
                        PCA19 + PCA20 + BMI + Townsend, data = filter(bd_pheno, Fish_oil_baseline == 0), na.action = "na.omit")

HDL_Female_lm <- lm(HDL ~ df_HDL + Sex + Age + I(Age^2) + Array + 
                         Assessment_centres_10003 + Assessment_centres_11001 + Assessment_centres_11002 + 
                         Assessment_centres_11003 + Assessment_centres_11004 + Assessment_centres_11005 + 
                         Assessment_centres_11006 + Assessment_centres_11007 + Assessment_centres_11008 + 
                         Assessment_centres_11009 + Assessment_centres_11010 + Assessment_centres_11011 + 
                         Assessment_centres_11011 + Assessment_centres_11012 + Assessment_centres_11013 + 
                         Assessment_centres_11014 + Assessment_centres_11016 + Assessment_centres_11017 + 
                         Assessment_centres_11018 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + 
                         PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + 
                         PCA19 + PCA20 + BMI + Townsend, data = filter(bd_pheno, Fish_oil_baseline == 0), na.action = "na.omit")

LDL_Female_lm <- lm(LDL ~ df_LDL + Sex + Age + I(Age^2) + Array + 
                         Assessment_centres_10003 + Assessment_centres_11001 + Assessment_centres_11002 + 
                         Assessment_centres_11003 + Assessment_centres_11004 + Assessment_centres_11005 + 
                         Assessment_centres_11006 + Assessment_centres_11007 + Assessment_centres_11008 + 
                         Assessment_centres_11009 + Assessment_centres_11010 + Assessment_centres_11011 + 
                         Assessment_centres_11011 + Assessment_centres_11012 + Assessment_centres_11013 + 
                         Assessment_centres_11014 + Assessment_centres_11016 + Assessment_centres_11017 + 
                         Assessment_centres_11018 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + 
                         PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + 
                         PCA19 + PCA20 + BMI + Townsend, data = filter(bd_pheno, Fish_oil_baseline == 0), na.action = "na.omit")

TAG_Female_lm <- lm(TAGs ~ df_TAG + Sex + Age + I(Age^2) + Array + 
                         Assessment_centres_10003 + Assessment_centres_11001 + Assessment_centres_11002 + 
                         Assessment_centres_11003 + Assessment_centres_11004 + Assessment_centres_11005 + 
                         Assessment_centres_11006 + Assessment_centres_11007 + Assessment_centres_11008 + 
                         Assessment_centres_11009 + Assessment_centres_11010 + Assessment_centres_11011 + 
                         Assessment_centres_11011 + Assessment_centres_11012 + Assessment_centres_11013 + 
                         Assessment_centres_11014 + Assessment_centres_11016 + Assessment_centres_11017 + 
                         Assessment_centres_11018 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + 
                         PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + 
                         PCA19 + PCA20 + BMI + Townsend, data = filter(bd_pheno, Fish_oil_baseline == 0), na.action = "na.omit")

TC_Male_lm <- tidy(TC_Male_lm)
HDL_Male_lm <- tidy(HDL_Male_lm)
LDL_Male_lm <- tidy(LDL_Male_lm)
TAG_Male_lm <- tidy(TAG_Male_lm)

TC_Female_lm <- tidy(TC_Female_lm)
HDL_Female_lm <- tidy(HDL_Female_lm)
LDL_Female_lm <- tidy(LDL_Female_lm)
TAG_Female_lm <- tidy(TAG_Female_lm)


# Fish Interaction
TC_Fish_lm <- lm(Tot_Chol ~ df_TC + Sex + Age + I(Age^2) + Array + 
                        Assessment_centres_10003 + Assessment_centres_11001 + Assessment_centres_11002 + 
                        Assessment_centres_11003 + Assessment_centres_11004 + Assessment_centres_11005 + 
                        Assessment_centres_11006 + Assessment_centres_11007 + Assessment_centres_11008 + 
                        Assessment_centres_11009 + Assessment_centres_11010 + Assessment_centres_11011 + 
                        Assessment_centres_11011 + Assessment_centres_11012 + Assessment_centres_11013 + 
                        Assessment_centres_11014 + Assessment_centres_11016 + Assessment_centres_11017 + 
                        Assessment_centres_11018 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + 
                        PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + 
                        PCA19 + PCA20 + BMI + Townsend, data = filter(bd_pheno, Fish_oil_baseline == 1), na.action = "na.omit")

HDL_Fish_lm <- lm(HDL ~ df_HDL + Sex + Age + I(Age^2) + Array + 
                         Assessment_centres_10003 + Assessment_centres_11001 + Assessment_centres_11002 + 
                         Assessment_centres_11003 + Assessment_centres_11004 + Assessment_centres_11005 + 
                         Assessment_centres_11006 + Assessment_centres_11007 + Assessment_centres_11008 + 
                         Assessment_centres_11009 + Assessment_centres_11010 + Assessment_centres_11011 + 
                         Assessment_centres_11011 + Assessment_centres_11012 + Assessment_centres_11013 + 
                         Assessment_centres_11014 + Assessment_centres_11016 + Assessment_centres_11017 + 
                         Assessment_centres_11018 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + 
                         PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + 
                         PCA19 + PCA20 + BMI + Townsend, data = filter(bd_pheno, Fish_oil_baseline == 1), na.action = "na.omit")

LDL_Fish_lm <- lm(LDL ~ df_LDL + Sex + Age + I(Age^2) + Array + 
                         Assessment_centres_10003 + Assessment_centres_11001 + Assessment_centres_11002 + 
                         Assessment_centres_11003 + Assessment_centres_11004 + Assessment_centres_11005 + 
                         Assessment_centres_11006 + Assessment_centres_11007 + Assessment_centres_11008 + 
                         Assessment_centres_11009 + Assessment_centres_11010 + Assessment_centres_11011 + 
                         Assessment_centres_11011 + Assessment_centres_11012 + Assessment_centres_11013 + 
                         Assessment_centres_11014 + Assessment_centres_11016 + Assessment_centres_11017 + 
                         Assessment_centres_11018 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + 
                         PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + 
                         PCA19 + PCA20 + BMI + Townsend, data = filter(bd_pheno, Fish_oil_baseline == 1), na.action = "na.omit")

TAG_Fish_lm <- lm(TAGs ~ df_TAG + Sex + Age + I(Age^2) + Array + 
                         Assessment_centres_10003 + Assessment_centres_11001 + Assessment_centres_11002 + 
                         Assessment_centres_11003 + Assessment_centres_11004 + Assessment_centres_11005 + 
                         Assessment_centres_11006 + Assessment_centres_11007 + Assessment_centres_11008 + 
                         Assessment_centres_11009 + Assessment_centres_11010 + Assessment_centres_11011 + 
                         Assessment_centres_11011 + Assessment_centres_11012 + Assessment_centres_11013 + 
                         Assessment_centres_11014 + Assessment_centres_11016 + Assessment_centres_11017 + 
                         Assessment_centres_11018 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + 
                         PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + 
                         PCA19 + PCA20 + BMI + Townsend, data = filter(bd_pheno, Fish_oil_baseline == 1), na.action = "na.omit")


TC_No_Fish_lm <- lm(Tot_Chol ~ df_TC + Sex + Age + I(Age^2) + Array + 
                        Assessment_centres_10003 + Assessment_centres_11001 + Assessment_centres_11002 + 
                        Assessment_centres_11003 + Assessment_centres_11004 + Assessment_centres_11005 + 
                        Assessment_centres_11006 + Assessment_centres_11007 + Assessment_centres_11008 + 
                        Assessment_centres_11009 + Assessment_centres_11010 + Assessment_centres_11011 + 
                        Assessment_centres_11011 + Assessment_centres_11012 + Assessment_centres_11013 + 
                        Assessment_centres_11014 + Assessment_centres_11016 + Assessment_centres_11017 + 
                        Assessment_centres_11018 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + 
                        PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + 
                        PCA19 + PCA20 + BMI + Townsend, data = filter(bd_pheno, Fish_oil_baseline == 0), na.action = "na.omit")

HDL_No_Fish_lm <- lm(HDL ~ df_HDL + Sex + Age + I(Age^2) + Array + 
                         Assessment_centres_10003 + Assessment_centres_11001 + Assessment_centres_11002 + 
                         Assessment_centres_11003 + Assessment_centres_11004 + Assessment_centres_11005 + 
                         Assessment_centres_11006 + Assessment_centres_11007 + Assessment_centres_11008 + 
                         Assessment_centres_11009 + Assessment_centres_11010 + Assessment_centres_11011 + 
                         Assessment_centres_11011 + Assessment_centres_11012 + Assessment_centres_11013 + 
                         Assessment_centres_11014 + Assessment_centres_11016 + Assessment_centres_11017 + 
                         Assessment_centres_11018 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + 
                         PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + 
                         PCA19 + PCA20 + BMI + Townsend, data = filter(bd_pheno, Fish_oil_baseline == 0), na.action = "na.omit")

LDL_No_Fish_lm <- lm(LDL ~ df_LDL + Sex + Age + I(Age^2) + Array + 
                         Assessment_centres_10003 + Assessment_centres_11001 + Assessment_centres_11002 + 
                         Assessment_centres_11003 + Assessment_centres_11004 + Assessment_centres_11005 + 
                         Assessment_centres_11006 + Assessment_centres_11007 + Assessment_centres_11008 + 
                         Assessment_centres_11009 + Assessment_centres_11010 + Assessment_centres_11011 + 
                         Assessment_centres_11011 + Assessment_centres_11012 + Assessment_centres_11013 + 
                         Assessment_centres_11014 + Assessment_centres_11016 + Assessment_centres_11017 + 
                         Assessment_centres_11018 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + 
                         PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + 
                         PCA19 + PCA20 + BMI + Townsend, data = filter(bd_pheno, Fish_oil_baseline == 0), na.action = "na.omit")

TAG_No_Fish_lm <- lm(TAGs ~ df_TAG + Sex + Age + I(Age^2) + Array + 
                         Assessment_centres_10003 + Assessment_centres_11001 + Assessment_centres_11002 + 
                         Assessment_centres_11003 + Assessment_centres_11004 + Assessment_centres_11005 + 
                         Assessment_centres_11006 + Assessment_centres_11007 + Assessment_centres_11008 + 
                         Assessment_centres_11009 + Assessment_centres_11010 + Assessment_centres_11011 + 
                         Assessment_centres_11011 + Assessment_centres_11012 + Assessment_centres_11013 + 
                         Assessment_centres_11014 + Assessment_centres_11016 + Assessment_centres_11017 + 
                         Assessment_centres_11018 + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + 
                         PCA9 + PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA16 + PCA17 + PCA18 + 
                         PCA19 + PCA20 + BMI + Townsend, data = filter(bd_pheno, Fish_oil_baseline == 0), na.action = "na.omit")

TC_Fish_lm <- tidy(TC_Fish_lm)
HDL_Fish_lm <- tidy(HDL_Fish_lm)
LDL_Fish_lm <- tidy(LDL_Fish_lm)
TAG_Fish_lm <- tidy(TAG_Fish_lm)

TC_No_Fish_lm <- tidy(TC_No_Fish_lm)
HDL_No_Fish_lm <- tidy(HDL_No_Fish_lm)
LDL_No_Fish_lm <- tidy(LDL_No_Fish_lm)
TAG_No_Fish_lm <- tidy(TAG_No_Fish_lm)

