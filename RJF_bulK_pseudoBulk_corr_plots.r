library("ggpubr")
library("dplyr")
library("ggplot2")
library("ggExtra")
library('cowplot')
require("gridExtra")


setwd('/Users/reza/Documents/methylation_project/IR_TCells/2nd_Run/03_Bulk_vs._pseudo_Bulk/Bulk_vs._pseudoBulk_corr')

# filenames <- list.files(path=getwd())

# for (i in filenames){
#   df_name <- strsplit(i,'_')[[1]][2]
#   assign(df_name,read.csv(i, header=FALSE, sep="\t", col.names = c('pseudo', 'bulk')))
#   #names(df_name) <- c('pseudo', 'bulk')
#   
#   #df_nonZero <- filter(df_name, bulk > 0, pseudo > 0)
#   
#   ggscatter(CFC, x = "bulk", y = "pseudo", 
#             add = "reg.line", conf.int = TRUE, 
#             cor.coef = TRUE, cor.method = "pearson",
#             xlab = "Bulk sample", ylab = "Pseudo Bulk",
#             size= .05, color = "pink", shape = 21,
#             add.params = list(color = "black", fill = "lightgray", size= .1),
#             cor.coeff.args = list(method = "pearson", label.sep = "\n"))
# 
# }



CFC_df <- read.csv("merged_CFC_pseudoBulk_vs._bulk_methylation_level_CpG_bin_3k_Homo_sapiens_sorted.txt", header=FALSE, sep="\t")
names(CFC_df) <- c('pseudo', 'bulk')

CFC_df_nonZero <- filter(CFC_df, bulk > 0, pseudo > 0)

CFC_p <- ggscatter(CFC_df_nonZero[sample(nrow(CFC_df_nonZero), 10000),],
                   x = "bulk", y = "pseudo", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "Bulk sample", ylab = "Pseudo Bulk", title = "CFC",
                   size= .05, color = "darkblue", shape = 21,
                   add.params = list(color = "black", fill = "lightgray", size= .1),
                   cor.coeff.args = list(method = "pearson", label.sep = "\n"))

ggMarginal(CFC_p, type = "density")
###################################

CFL_df <- read.csv("merged_CFL_pseudoBulk_vs._bulk_methylation_level_CpG_bin_3k_Homo_sapiens_sorted.txt", header=FALSE, sep="\t")
names(CFL_df) <- c('pseudo', 'bulk')

CFL_df_nonZero <- filter(CFL_df, bulk > 0, pseudo > 0)

CFL_p <- ggscatter(CFL_df_nonZero[sample(nrow(CFL_df_nonZero), 10000),],
                   x = "bulk", y = "pseudo", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "Bulk sample", ylab = "Pseudo Bulk", title = "CFL",
                   size= .05, color = "darkblue", shape = 21,
                   add.params = list(color = "black", fill = "lightgray", size= .1),
                   cor.coeff.args = list(method = "pearson", label.sep = "\n"))

ggMarginal(CFL_p, type = "density")
###################################

OFC_df <- read.csv("merged_OFC_pseudoBulk_vs._bulk_methylation_level_CpG_bin_3k_Homo_sapiens_sorted.txt", header=FALSE, sep="\t")
names(OFC_df) <- c('pseudo', 'bulk')

OFC_df_nonZero <- filter(OFC_df, bulk > 0, pseudo > 0)

OFC_p <- ggscatter(OFC_df_nonZero[sample(nrow(OFC_df_nonZero), 10000),],
                   x = "bulk", y = "pseudo", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "Bulk sample", ylab = "Pseudo Bulk", title = "OFC",
                   size= .05, color = "darkblue", shape = 21,
                   add.params = list(color = "black", fill = "lightgray", size= .1),
                   cor.coeff.args = list(method = "pearson", label.sep = "\n"))

ggMarginal(OFC_p, type = "density")
###################################

OFL_df <- read.csv("merged_OFL_pseudoBulk_vs._bulk_methylation_level_CpG_bin_3k_Homo_sapiens_sorted.txt", header=FALSE, sep="\t")
names(OFL_df) <- c('pseudo', 'bulk')

OFL_df_nonZero <- filter(OFL_df, bulk > 0, pseudo > 0)

OFL_p <- ggscatter(OFL_df_nonZero[sample(nrow(OFL_df_nonZero), 10000),],
                   x = "bulk", y = "pseudo", 
                   add = "reg.line", conf.int = TRUE,
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "Bulk sample", ylab = "Pseudo Bulk", title = "OFL",
                   size= .05, color = "darkblue", shape = 21,
                   add.params = list(color = "black", fill = "lightgray", size= .1),
                   cor.coeff.args = list(method = "pearson", label.sep = "\n"))

ggMarginal(OFL_p, type = "density")


###################################
YFC_df <- read.csv("merged_YFC_pseudoBulk_vs._bulk_methylation_level_CpG_bin_3k_Homo_sapiens_sorted.txt", header=FALSE, sep="\t")
names(YFC_df) <- c('pseudo', 'bulk')

YFC_df_nonZero <- filter(YFC_df, bulk > 0, pseudo > 0)

YFC_p <- ggscatter(YFC_df_nonZero[sample(nrow(YFC_df_nonZero), 10000),],
                   x = "bulk", y = "pseudo", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "Bulk sample", ylab = "Pseudo Bulk", title = "YFC",
                   size= .05, color = "darkblue", shape = 21,
                   add.params = list(color = "black", fill = "lightgray", size= .1),
                   cor.coeff.args = list(method = "pearson", label.sep = "\n"))

ggMarginal(YFC_p, type = "density")
###################################
YFL_df <- read.csv("merged_YFL_pseudoBulk_vs._bulk_methylation_level_CpG_bin_3k_Homo_sapiens_sorted.txt", header=FALSE, sep="\t")
names(YFL_df) <- c('pseudo', 'bulk')

YFL_df_nonZero <- filter(YFL_df, bulk > 0, pseudo > 0)

YFL_p <- ggscatter(YFL_df_nonZero[sample(nrow(YFL_df_nonZero), 10000),],
                   x = "bulk", y = "pseudo", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "Bulk sample", ylab = "Pseudo Bulk", title = "YFL",
                   size= .05, color = "darkblue", shape = 21,
                   add.params = list(color = "black", fill = "lightgray", size= .1),
                   cor.coeff.args = list(method = "pearson", label.sep = "\n"))

ggMarginal(YFL_p, type = "density")

### Creating the plot and saving

Figure_corr <- plot_grid(CFC_p, CFL_p, OFC_p, OFL_p, YFC_p, YFL_p, labels = "AUTO", ncol = 2, nrow= 3, align = 'h')

png("/Users/reza/Documents/methylation_project/IR_TCells/2nd_Run/Plots/Bulk_pseudoBulk_corr.png", width = 30, height = 30, units = 'cm', res = 300)
grid.arrange(Figure_corr) # Make plot
dev.off()

### Noise-variance plot
library(ggrepel)
setwd('/Users/reza/Documents/methylation_project/IR_TCells/2nd_Run/04_DMW_variance_noise')
variance_noise_df <- read.csv("IR-TCells_noise_variance_DMWs.txt", sep='\t')

# Creat a new column follwiong Xiao's suggestion on 12/09/2020 (mm/dd/yyyy):
variance_noise_df$diff= variance_noise_df$variance - variance_noise_df$noise


ggplot(variance_noise_df, aes(x=noise, y=variance, color=batch, label=Pair)) + 
  geom_point()+ geom_text(vjust = 0, nudge_y = 0.)+ ggtitle("IR_TCells variance-noise plot") +theme_bw()

p <- ggplot(variance_noise_df, aes(x=noise, y=variance, label=Pair)) + geom_point(color = "red")+ theme_bw()
p + geom_text_repel(min.segment.length = 0, seed = 42, box.padding = 0.5) + labs(title = "IR_TCells variance-noise plot")

ggscatter(variance_noise_df, x = "noise", y = "variance",
          add = "reg.line",
          conf.int = FALSE,
          color = "batch", palette = "jco",
          label = "Pair", repel = TRUE,
          shape = "batch")+ stat_cor(aes(color = batch))

# (Noise)
my_comparisons <- list( c("YFC", "YFL"), c("OFC", "OFL"), c("CFC", "CFL"), c("YFC", "OFC"), c("YFC", "CFC"), c("YFL", "OFL"),
                         c("YFL", "CFL"))

noise_plot <- ggboxplot(variance_noise_df, x = "batch", y = "noise",
                        color = "batch", palette = c("#E7B800", "#e7b90090", "#0073C2FF", "#0074c2ba", "#868686FF", "#86868658"),
                        add = "jitter")+stat_compare_means(comparisons = my_comparisons)+
                        scale_x_discrete(limits=c("YFC", "YFL", "OFC", "OFL", "CFC", "CFL"))+
                        xlab("Treatment")


# ggboxplot(variance_noise_df, x = "batch", y = "noise",
#           color = "batch", palette = "jco",
#           add = "jitter")+stat_compare_means(comparisons = my_comparisons)

# (Variance)

variance_plot <- ggboxplot(variance_noise_df, x = "batch", y = "variance",
                        color = "batch", palette = c("#E7B800", "#e7b90090", "#0073C2FF", "#0074c2ba", "#868686FF", "#86868658"),
                        add = "jitter")+stat_compare_means(comparisons = my_comparisons)+
                        scale_x_discrete(limits=c("YFC", "YFL", "OFC", "OFL", "CFC", "CFL"))+
                        xlab("Treatment")

# ggboxplot(variance_noise_df, x = "batch", y = "variance",
#           color = "batch", palette = "jco",
#           add = "jitter")+stat_compare_means(comparisons = my_comparisons)

# ggboxplot(variance_noise_df, x = "batch", y = "variance",
#           color = "batch", palette = "jco")+
#   stat_compare_means(method = "anova")+      # Add global p-value
#   stat_compare_means(label = "p.signif", method = "t.test",
#                      ref.group = ".all.") + stat_compare_means(comparisons = my_comparisons)

### DMW plot

DMW_plot <- ggboxplot(variance_noise_df, x = "batch", y = "DMW",
                        color = "batch", palette = c("#E7B800", "#e7b90090", "#0073C2FF", "#0074c2ba", "#868686FF", "#86868658"),
                        add = "jitter")+stat_compare_means(comparisons = my_comparisons)+
                        scale_x_discrete(limits=c("YFC", "YFL", "OFC", "OFL", "CFC", "CFL"))+
                        xlab("Treatment")

# ggboxplot(variance_noise_df, x = "batch", y = "DMW",
#           color = "batch", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#           add = "jitter")+ stat_compare_means(comparisons = my_comparisons)


# ggboxplot(variance_noise_df2, x = "batch", y = "DMW",
#           color = "batch", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#           add = "jitter")+ stat_compare_means(comparisons = my_comparisons)

diff_plot <- ggboxplot(variance_noise_df, x = "batch", y = "diff",
                      color = "batch", palette = c("#E7B800", "#e7b90090", "#0073C2FF", "#0074c2ba", "#868686FF", "#86868658"),
                      # palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                      add = "jitter")+ stat_compare_means(comparisons = my_comparisons)+
                      scale_x_discrete(limits=c("YFC", "YFL", "OFC", "OFL", "CFC", "CFL"))+
                      xlab("Treatment")+
                      ylab("(variance - noise)")


Figure_var_noise <- plot_grid(variance_plot, noise_plot, diff_plot, DMW_plot, labels = "AUTO", ncol = 4, align = 'h')

png("/Users/reza/Documents/methylation_project/IR_TCells/2nd_Run/Plots/IR-TCells_variance_noise_DMWs_modified.png", width = 40, height = 15, units = 'cm', res = 300)
grid.arrange(Figure_var_noise) # Make plot
dev.off()


### Methylation percentage violine plot
library("ggpubr")
meth_prc_df <- read.csv('/Users/reza/Documents/methylation_project/IR_TCells/2nd_Run/Plots/IR_TCells_BATCH_methylation_percentage_stat_report.txt', sep = '\t')

meth_prc <- ggviolin(meth_prc_df, x = "batch", y = "CpG_meth", fill = "batch",
                    palette = c("#E7B800", "#e7b90090", "#0073C2FF", "#0074c2ba", "#868686FF", "#86868658"),
                    add = "boxplot", add.params = list(fill = "white"))+
                    stat_compare_means(comparisons = my_comparisons)+
                    scale_x_discrete(limits=c("YFC", "YFL", "OFC", "OFL", "CFC", "CFL"))+
                    xlab("Treatment")+
                    ylab("%CpG methylation")

nonCpG_meth_prc <- ggviolin(meth_prc_df, x = "batch", y = "nonCpG_meth", fill = "batch",
                    palette = c("#E7B800", "#e7b90090", "#0073C2FF", "#0074c2ba", "#868686FF", "#86868658"),
                    add = "boxplot", add.params = list(fill = "white"))+
                    stat_compare_means(comparisons = my_comparisons)+
                    scale_x_discrete(limits=c("YFC", "YFL", "OFC", "OFL", "CFC", "CFL"))+
                    xlab("Treatment")+
                    ylab("%non.CpG methylation")

Figure_meth_prc <- plot_grid(meth_prc, labels = "AUTO", ncol = 1, align = 'h')

png("/Users/reza/Documents/methylation_project/IR_TCells/2nd_Run/Plots/IR-TCells_methylation_prc.png", width = 15, height = 15, units = 'cm', res = 300)
grid.arrange(Figure_meth_prc) # Make plot
dev.off()


# Bar plot of mean +/-se
ggbarplot(variance_noise_df, x = "batch", y = "DMW", add = "mean_se")+
  stat_compare_means()

# Line plot of mean +/-se
ggline(variance_noise_df, x = "batch", y = "DMW", add = "mean_se")+
  stat_compare_means()
