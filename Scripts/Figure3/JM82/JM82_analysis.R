library(barcodetrackR)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(SummarizedExperiment)

# read in the modified heatmap function 
source("/Volumes/dirgroup/TSCBB/LAB/Annie/5_Tissue_NK_Study/Figure4/ZJ31/barcode_ggheatmap_annie_v2.1.R")

# figure number 
figure_num <- 3

# input information about monkey and samples
monkey_ID <- "JM82"

######### INITIAL SETUP (Loading in files + locating samples) #################

data_dir <- paste("/Volumes/dirgroup/TSCBB/LAB/Annie/6_Tissue_NK_Manuscript/raw_data", sep = "")
save_dir <- paste("/Volumes/dirgroup/TSCBB/LAB/Annie/6_Tissue_NK_Manuscript/", "Figure", figure_num, sep = "")

# read in sample list 
sample_list <- read.delim(file = file.path(save_dir,paste(monkey_ID, "/", monkey_ID, "_B_samples.txt", sep="")),
                          sep = "\t", header=F)
colnames(sample_list) <- "sample_name"
sample_list <- sample_list$sample_name

# read in metadata file 
metadata_df <- read.delim(file = file.path(data_dir,paste("counts_data/", monkey_ID, "_sample_metadata.txt", sep = "")),
                          sep = "\t", header = T)

# read in counts file 
count_data <- read.delim(file = file.path(data_dir,paste("counts_data/", monkey_ID, "_sample_counts.txt", sep = "")),
                         sep = "\t", header = T, row.names = 1)
count_data <- count_data[,metadata_df$SAMPLENAME]

# Make sure everything is numeric
count_data <- data.frame(apply(count_data, 2,            
                               function(x) as.numeric(as.character(x))),
                         row.names = rownames(count_data))

# find the indices corresponding to each sample in the metadata & countdata  
sample_index <- rep(NA, length(sample_list))
for (i in 1:length(sample_list)){
  sample_index[i] <- which(metadata_df$SAMPLENAME == sample_list[i])
}

# creating plot labels 
plot_labels_B <- paste( paste(metadata_df[sample_index, ]$Timepoint, "m", sep = ""), metadata_df[sample_index, ]$Experiment, metadata_df[sample_index, ]$Cell.Type)

metadata_df[["cell_source"]] <- paste(metadata_df$Experiment, metadata_df$Cell.Type)
metadata_df[["source"]] <- metadata_df$Experiment
metadata_df$Experiment <- paste(metadata_df$source, metadata_df$Cell.Type)

# loading in the data into BarcodetrackR
my_SE <- barcodetrackR::create_SE(your_data = count_data, meta_data = metadata_df, threshold = 0)


#################HEATMAPS + CORRELATION + DIVERSITY PLOTS ########################

# Top 10 clones heatmap 
barcodetrackR::barcode_ggheatmap(my_SE[,sample_index], n_clones = 10, label_size = 10, 
                                 plot_labels = plot_labels_B) # regular heatmap 

# Modified Top 10 clones heatmap (Fold change = 5)
barcode_ggheatmap_annie_v2.1(my_SE[,sample_index], n_clones = 10, label_size = 25, 
                             plot_labels = plot_labels_B, 
                             tissue_specific_clones = 1, foldchange_threshold = 5)

# Correlation plot (negatives turned off)
cor_plot(my_SE[,sample_index], method_corr = "pearson", 
         plot_type = "color", assay = "proportions", 
         plot_labels = plot_labels_B, label_size = 25, no_negatives = T)

# Diversity plots 
clonal_diversity(my_SE[,sample_index],
                 plot_over = "Timepoint",
                 group_by = "Experiment",
                 index_type = "shannon", text_size = 30)


############################# PCA PLOT #######################################

coldata_SE <- SummarizedExperiment::colData(my_SE)

normalized_values <- SummarizedExperiment::assays(my_SE)[["normalized"]]

norm_val_samples <- t(normalized_values[, sample_index])
all_data_pca <- prcomp(norm_val_samples)
all_data_pca_df <- as.data.frame(all_data_pca$x)
all_data_pca_df$source <- my_SE[,sample_index]$Experiment
all_data_pca_df$cell_type <- my_SE[,sample_index]$Cell.Type

timepoints <- factor(metadata_df[sample_index, ]$Timepoint)
location <- factor(my_SE[,sample_index]$Experiment)
cell.type <- factor(my_SE[,sample_index]$Cell.Type)

p_all <- ggplot(all_data_pca_df,aes(x=PC1,y=PC2, shape = timepoints, color = location)) + 
  geom_point(alpha = 0.6, size = 12) + #stat_ellipse(aes(x=PC1, y=PC2,color=location),type = "norm") + 
  labs(shape = "timepoints", color = "Tissue", x = "PC1", y = "PC2") + 
  scale_color_manual(values = c("PB Gr" = "grey", "PB B" = "#F8766D", "BAL B" = "#B79F00", "Liver B" = "#00BA38", 
                                "Spleen B" = "#00BFC4", "Jejunum B" = "#619CFF", "Colon B" = "#F564E3")) + 
  scale_shape_manual(values = c("45.5" = 16, "46" = 17, "47" = 15, "48" = 18, "49" = 8, "49.5" = 7)) +   
  theme(
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 25),
    legend.text=element_text(size=25),
    legend.title = element_text(size=25), 
    legend.position="right", 
    legend.box = "vertical"
  )  + guides(color = "none")
print(p_all)