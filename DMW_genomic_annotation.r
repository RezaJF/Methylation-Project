library(AnnotationHub)
library(annotatr)
library(regioneR)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(svMisc)

setwd("/Users/reza/Documents/methylation_project/IR_TCells/2nd_Run/04_DMW_variance_noise/DMR")

################## STEP 01: PREPARING THE INPUT FILES FOR ANNOTATION ##################

#################
### CFC Batch ###
#################

# Get a list of all text files with the "dmr-CFCB" prefix
# CFC_file_list <- list.files(pattern = "dmr-CFCB.*\\.txt")

# Read the header from the first text file
# CFC_header <- readLines(CFC_file_list[1], n = 1, sep="\t")

# Read the rest of the files, skipping the header
# CFC_data_list <- lapply(CFC_file_list[-1], read.table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Combine the data from all files into a single dataframe
# CFC_data <- do.call(rbind, CFC_data_list)

# Add the header as the column names of the dataframe
# colnames(CFC_data) <- c("chr", "binid", "binstart",	"binend", "meth.x", "unmeth.x",	"meth.y", "unmeth.y",
#                         "m1", "n1", "m2", "n2", "se1", "se2", "sepool", "s", "cid", "ciu", "p",	"statistic", "diff")

# CFC_data %>% distinct() %>% dplyr::select("chr", "binstart", "binend", "p", "diff", "m1", "m2") %>%
#         filter( p <= .05/nrow(CFC_data)) %>% 
#         arrange(chr,binstart) -> CFC_DMW_df

# CFC_DMW_df$chr <- sub("^", "chr", CFC_DMW_df$chr )

# write.table(CFC_DMW_df, file = "./genomic_annotation/CFC_batch_DMW.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names =  FALSE)


pattern_list = list("dmr-YFLB.*\\.txt", "dmr-YFCB.*\\.txt","dmr-OFCB.*\\.txt", "dmr-OFLB.*\\.txt", "dmr-CFCB.*\\.txt", "dmr-CFLB.*\\.txt")
file_names = list("YFL_batch_DMW.txt", "YFC_batch_DMW.txt", "OFC_batch_DMW.txt", "OFL_batch_DMW.txt", "CFC_batch_DMW.txt", "CFL_batch_DMW.txt")

for (i in 1:length(pattern_list)) {
  file_list <- list.files(pattern = paste(pattern_list[i]))
  data_list <- lapply(file_list[-1], read.table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  data <- do.call(rbind, data_list)

  data %>% distinct() %>% dplyr::select("chr", "binstart", "binend", "p", "diff", "m1", "m2") %>%
        filter( p <= .05/nrow(data)) %>% 
        arrange(chr,binstart) -> DMW_df
  
  DMW_df$chr <- sub("^", "chr", DMW_df$chr )

  write.table(DMW_df, file = paste("./genomic_annotation/",file_names[i], sep=""), sep = "\t", row.names = FALSE, quote = FALSE, col.names =  FALSE)
  progress(i)
}

################## STEP 02: ANNOTATING DIFFERENTIALLY METHYLATED WINDOWS ##################

# Create a named vector for the AnnotationHub accession codes with desired names
# These annotation tracks are specific to T-Cells:

H3K4me3_codes = c('Dnd41' = 'AH29695')
H3K27ac_codes = c('Dnd41' = 'AH31682')
H3K27me3_codes = c('Dnd41' = 'AH31683')
H3K4me1_codes = c('Dnd41' = 'AH29693')
H3K36me3_codes = c('Dnd41' = 'AH29700')
H3K9me3_codes = c('Dnd41' = 'AH30737')
H3K9ac_codes = c('Dnd41' = 'AH29696')
H3K4me2_codes = c('Dnd41' = 'AH29694')
H4K20me1_codes = c('Dnd41' = 'AH29702')
H3K79me2_codes = c('Dnd41' = 'AH29701')

# Fetch ah_codes from AnnotationHub and create annotations annotatr understands
histone_marks= c('H3K4me3', 'H3K27ac', 'H3K27me3', 'H3K4me1', 'H3K36me3',
                 'H3K9me3', 'H3K9ac', 'H3K4me2', 'H4K20me1', 'H3K79me2')

build_ah_annots(genome = 'hg19', ah_codes = H3k4me3_codes, annotation_class = "H3K4me3")
build_ah_annots(genome = 'hg19', ah_codes = H3K27ac_codes, annotation_class = "H3K27ac")
build_ah_annots(genome = 'hg19', ah_codes = H3K27me3_codes, annotation_class = "H3K27me3")
build_ah_annots(genome = 'hg19', ah_codes = H3K4me1_codes, annotation_class = "H3K4me1")
build_ah_annots(genome = 'hg19', ah_codes = H3K36me3_codes, annotation_class = "H3K36me3")
build_ah_annots(genome = 'hg19', ah_codes = H3K9me3_codes, annotation_class = "H3K9me3")
build_ah_annots(genome = 'hg19', ah_codes = H3K9ac_codes, annotation_class = "H3K9ac")
build_ah_annots(genome = 'hg19', ah_codes = H3K4me2_codes, annotation_class = "H3K4me2")
build_ah_annots(genome = 'hg19', ah_codes = H4K20me1_codes, annotation_class = "H4K20me1")
build_ah_annots(genome = 'hg19', ah_codes = H3K79me2_codes, annotation_class = "H3K79me2")


ah_names = c('hg19_H3K4me3_Dnd41', 'hg19_H3K27ac_Dnd41', 'hg19_H3K27me3_Dnd41', 'hg19_H3K4me1_Dnd41',
             'hg19_H3K36me3_Dnd41', 'hg19_H3K9me3_Dnd41', 'hg19_H3K9ac_Dnd41', 'hg19_H3K4me2_Dnd41',
             'hg19_H4K20me1_Dnd41', 'hg19_H3K79me2_Dnd41')


dm_regions_YFC = read_regions(con = "./genomic_annotation/YFC_batch_DMW.txt", genome = 'hg19', extraCols = c(p = 'numeric', diff= 'numeric', m1= 'numeric', m2= 'numeric'), format = 'bed')
dm_regions_YFL = read_regions(con = "./genomic_annotation/YFL_batch_DMW.txt", genome = 'hg19', extraCols = c(p = 'numeric', diff= 'numeric', m1= 'numeric', m2= 'numeric'), format = 'bed')

dm_regions_OFC = read_regions(con = "./genomic_annotation/OFC_batch_DMW.txt", genome = 'hg19', extraCols = c(p = 'numeric', diff= 'numeric', m1= 'numeric', m2= 'numeric'), format = 'bed')
dm_regions_OFL = read_regions(con = "./genomic_annotation/OFL_batch_DMW.txt", genome = 'hg19', extraCols = c(p = 'numeric', diff= 'numeric', m1= 'numeric', m2= 'numeric'), format = 'bed')

dm_regions_CFC = read_regions(con = "./genomic_annotation/CFC_batch_DMW.txt", genome = 'hg19', extraCols = c(p = 'numeric', diff= 'numeric', m1= 'numeric', m2= 'numeric'), format = 'bed')
dm_regions_CFL = read_regions(con = "./genomic_annotation/CFL_batch_DMW.txt", genome = 'hg19', extraCols = c(p = 'numeric', diff= 'numeric', m1= 'numeric', m2= 'numeric'), format = 'bed')




# Select annotations for intersection with regions
# Note inclusion of custom annotation, and use of shortcuts
annots = c('hg19_cpg_shores', 'hg19_cpg_islands' , 'hg19_cpg_shelves', 
           'hg19_basicgenes', 'hg19_genes_intergenic', 'hg19_genes_3UTRs', 'hg19_genes_exons', 'hg19_genes_introns',
           'hg19_genes_intronexonboundaries', 'hg19_enhancers_fantom', 'hg19_genes_intergenic', 
           'hg19_genes_5UTRs', 'hg19_genes_promoters',
           'hg19_chromatin_Nhlf-ActivePromoter', 'hg19_chromatin_Nhlf-PoisedPromoter',
           'hg19_chromatin_Nhlf-WeakEnhancer', 'hg19_chromatin_Nhlf-Heterochrom/lo',
           'hg19_Nhlf-chromatin',
           'hg19_H3K4me3_Dnd41', 'hg19_H3K27ac_Dnd41', 'hg19_H3K27me3_Dnd41', 'hg19_H3K4me1_Dnd41',
           'hg19_H3K36me3_Dnd41', 'hg19_H3K9me3_Dnd41', 'hg19_H3K9ac_Dnd41', 'hg19_H3K4me2_Dnd41',
           'hg19_H4K20me1_Dnd41', 'hg19_H3K79me2_Dnd41')

# Build the annotations (a single GRanges object)
annotations = build_annotations(genome = 'hg19', annotations = annots)

# Intersect the regions we read in with the annotations
## 1. YFC
dm_annotated_YFC = annotate_regions(
  regions = dm_regions_YFC,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

## 2. YFL
dm_annotated_YFL = annotate_regions(
  regions = dm_regions_YFL,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

## 3. OFC
dm_annotated_OFC = annotate_regions(
  regions = dm_regions_OFC,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

## 4. OFL
dm_annotated_OFL = annotate_regions(
  regions = dm_regions_OFL,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

## 5. CFC
dm_annotated_CFC = annotate_regions(
  regions = dm_regions_CFC,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

## 6. CFL
dm_annotated_CFL = annotate_regions(
  regions = dm_regions_CFL,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

  # The `annotate_regions()` function returns a `GRanges`, but it may be more convenient to manipulate a coerced `data.frame`.
# Coerce to a data.frame

df_dm_annotated_YFC = data.frame(dm_annotated_YFC)
df_dm_annotated_YFL = data.frame(dm_annotated_YFL)

df_dm_annotated_OFC = data.frame(dm_annotated_OFC)
df_dm_annotated_OFL = data.frame(dm_annotated_OFL)

df_dm_annotated_CFC = data.frame(dm_annotated_CFC)
df_dm_annotated_CFL = data.frame(dm_annotated_CFL)


# # Take the mean of the diff column across all regions occurring in an annotation.
## 1. YFC
 dm_numsum_YFC = summarize_numerical(
  annotated_regions = dm_annotated_YFC,
  by = c('annot.type', 'annot.id'),
  over = c('diff'),
  quiet = TRUE)

## 2. YFL
 dm_numsum_YFL = summarize_numerical(
  annotated_regions = dm_annotated_YFL,
  by = c('annot.type', 'annot.id'),
  over = c('diff'),
  quiet = TRUE)

## 3. OFC
dm_numsum_OFC = summarize_numerical(
  annotated_regions = dm_annotated_OFC,
  by = c('annot.type', 'annot.id'),
  over = c('diff'),
  quiet = TRUE)

## 4. OFL
dm_numsum_OFL = summarize_numerical(
  annotated_regions = dm_annotated_OFL,
  by = c('annot.type', 'annot.id'),
  over = c('diff'),
  quiet = TRUE)

## 5. CFC
dm_numsum_CFC = summarize_numerical(
  annotated_regions = dm_annotated_CFC,
  by = c('annot.type', 'annot.id'),
  over = c('diff'),
  quiet = TRUE)

## 6. CFL
dm_numsum_CFL = summarize_numerical(
  annotated_regions = dm_annotated_CFL,
  by = c('annot.type', 'annot.id'),
  over = c('diff'),
  quiet = TRUE)

dm_numsum_YFC$batch <- "YFC"
dm_numsum_YFL$batch <- "YFL"

dm_numsum_OFC$batch <- "OFC"
dm_numsum_OFL$batch <- "OFL"

dm_numsum_CFC$batch <- "CFC"
dm_numsum_CFL$batch <- "CFL"

DMW_annotation <- bind_rows(dm_numsum_YFC, dm_numsum_YFL,
                            dm_numsum_OFC, dm_numsum_OFL,
                            dm_numsum_CFC, dm_numsum_CFL)


ggplot(DMW_annotation, aes(annot.type, sd)) +
  geom_boxplot(aes(fill= batch), alpha=0.7, show.legend = TRUE, outlier.shape = NA) +
  theme_classic()+
  #scale_fill_viridis_d()+
  xlab("annotation category") + ylab("variance")+
  rotate_x_text(45)+
  coord_cartesian(ylim =  c(0, .30))+
  geom_signif(comparisons = list(c("YFC", "YFL"), c("OFC", "OFL"), c("CFC", "CFL")), 
              map_signif_level=TRUE)+
              scale_fill_manual(values = c("#E7B800", "#e7b90090", "#0073C2FF", "#0074c2ba", "#868686FF", "#86868658"),
              breaks=c("CFC", "CFL","OFC", "OFL","YFC", "YFL")) -> Figure_annotation

require("gridExtra")
library('cowplot')

Figure <- plot_grid(Figure_annotation, labels = "AUTO", ncol = 1, align = 'h')

png("/Users/reza/Documents/methylation_project/IR_TCells/2nd_Run/Plots/DMW_genomic_annotation.png", width = 75, height = 35, units = 'cm', res = 300)
grid.arrange(Figure) # Make plot
dev.off()



 