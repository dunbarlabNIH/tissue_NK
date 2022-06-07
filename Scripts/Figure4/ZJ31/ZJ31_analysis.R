library(barcodetrackR)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(SummarizedExperiment)

# read in the modified heatmap function 
source("/Volumes/dirgroup/TSCBB/LAB/Annie/5_Tissue_NK_Study/Figure4/ZJ31/barcode_ggheatmap_annie_v2.1.R")

# figure number 
figure_num <- 4

# input information about monkey and samples
monkey_ID <- "ZJ31"

######### INITIAL SETUP (Loading in files + locating samples) #################

data_dir <- paste("/Volumes/dirgroup/TSCBB/LAB/Annie/6_Tissue_NK_Manuscript/raw_data", sep = "")
save_dir <- paste("/Volumes/dirgroup/TSCBB/LAB/Annie/6_Tissue_NK_Manuscript/", "Figure", figure_num, sep = "")

# read in sample list 
sample_list <- read.delim(file = file.path(save_dir,paste(monkey_ID, "/", monkey_ID, "_NK_samples.txt", sep="")),
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
plot_labels_NK <- paste( paste(metadata_df[sample_index, ]$months, "m", sep = ""), metadata_df[sample_index, ]$source, metadata_df[sample_index, ]$cell_type)

metadata_df[["cell_type_source"]] <- paste(metadata_df$source, metadata_df$cell_type)
metadata_df[["Experiment"]] <- paste(metadata_df$source, metadata_df$cell_type)

# loading in the data into BarcodetrackR
my_SE <- barcodetrackR::create_SE(your_data = count_data, meta_data = metadata_df, threshold = 0)


#################HEATMAPS + CORRELATION + DIVERSITY PLOTS ########################

# Top 10 clones heatmap 
barcodetrackR::barcode_ggheatmap(my_SE[,sample_index], n_clones = 10, label_size = 10, 
                                 plot_labels = plot_labels_NK) # regular heatmap 

# Modified Top 10 clones heatmap (Fold change = 5)
barcode_ggheatmap_annie_v2.1(my_SE[,sample_index], n_clones = 10, label_size = 25, 
                             plot_labels = plot_labels_NK, 
                             tissue_specific_clones = 1, foldchange_threshold = 5)

# Correlation plot (negatives turned off)
cor_plot(my_SE[,sample_index], method_corr = "pearson", 
         plot_type = "color", assay = "proportions", 
         plot_labels = plot_labels_NK, label_size = 25, no_negatives = T)

# Diversity plots 
clonal_diversity(my_SE[,sample_index],
                 plot_over = "months",
                 group_by = "Experiment",
                 index_type = "shannon", text_size = 30)


############################# PCA PLOT #######################################


coldata_SE <- SummarizedExperiment::colData(my_SE)

normalized_values <- SummarizedExperiment::assays(my_SE)[["normalized"]]

norm_val_samples <- t(normalized_values[, sample_index])
all_data_pca <- prcomp(norm_val_samples)
all_data_pca_df <- as.data.frame(all_data_pca$x)
all_data_pca_df$source <- my_SE[,sample_index]$Experiment
all_data_pca_df$cell_type <- my_SE[,sample_index]$cell_type

timepoints <- factor(metadata_df[sample_index, ]$months)
location <- factor(my_SE[,sample_index]$Experiment)
cell.type <- factor(my_SE[,sample_index]$cell_type)


p_all <- ggplot(all_data_pca_df,aes(x=PC1,y=PC2, shape = timepoints, color = location)) + 
  geom_point(size = 8, alpha = 0.75) + #stat_ellipse(aes(x=PC1, y=PC2,color=location),type = "norm") + 
  labs(shape = "Timepoints",color = "Tissue", x = "PC1", y = "PC2") + 
  scale_shape_manual(values = c("51" = 16, "63" = 17, "73.5" = 15, "75" = 18, "76" = 8)) + 
  scale_color_manual(values = c("PB Gr" = "grey", "PB CD16p" = "#F8766D", "PB CD56p" = "pink", "PB DN" = "orange", "BAL NK" = "#B79F00", "Liver NK" = "#00BA38", 
                                "Spleen NK" = "#00BFC4", "Jejunum NK" = "#619CFF", "Colon NK" = "#F564E3")) + 
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
  )  + guides(size = "none")
print(p_all)

######################### CUM CONTRIBUTION OF AGGREGATE TOP X CLONES ######################

sample_list <- sample_list[-c(1,2,3)]
sample_index <- sample_index[-c(1,2,3)]

top_x_clone_num <- 10

PB_CD16p_samples <- sample_list[1:4]
PB_CD56p_samples <- sample_list[5:8]
PB_DN_samples <- sample_list[9:10]
BAL_samples <- sample_list[11]
Liver_samples <- sample_list[12:13]
Spleen_samples <- sample_list[14]
Jejunum_samples <- sample_list[15:16]


PB_samples <- list("PB CD16p" = PB_CD16p_samples, "PB CD56p" = PB_CD56p_samples, 
                   "PB DN" = PB_DN_samples, "BAL NK" = BAL_samples, 
                   "Liver NK" = Liver_samples, "Spleen NK" = Spleen_samples, 
                   "Jejunum NK" = Jejunum_samples)

top_clones_list <- list()

for (i in 1:length(PB_samples)){
  
  top_clones_temp <- c()
  
  for (j in 1:length(PB_samples[[i]])){
    
    top_clones_temp_temp <- get_top_clones(my_SE, PB_samples[[i]][j], top_x_clone_num)
    top_clones_temp <- c(top_clones_temp, top_clones_temp_temp)
    
  }
  
  top_clones_df <- data.frame(unique(top_clones_temp))
  colnames(top_clones_df) <- names(PB_samples)[i]
  
  top_clones_list <- append(top_clones_list, top_clones_df)
  
}


sample_proportions <- SummarizedExperiment::assays(my_SE[,sample_index])[["proportions"]]

clones_biased_list <- top_clones_list
for (i in 1:length(clones_biased_list)){
  tissue <- names(clones_biased_list)[i]
  sample_prop_sum <- colSums(sample_proportions[clones_biased_list[[i]], ]) * 100
  sample_prop_sum <- data.frame(sample_prop_sum, source = my_SE$source[sample_index], 
                                cell.type = my_SE$cell_type[sample_index], 
                                cell.type.source = my_SE$cell_type_source[sample_index], 
                                class = factor(c(rep("PB", 10), rep("Tissue", 6)), levels = c("PB", "Tissue", "Sample")))
  print(tissue)
  sample_prop_sum$class[sample_prop_sum$cell.type.source == tissue] <- "Sample"
  sample_prop_sum$source <- factor(sample_prop_sum$source, levels = c("PB", "BAL", "Liver", "Spleen", "Jejunum", "Colon"))
  sample_prop_sum$cell.type.source <- factor(sample_prop_sum$cell.type.source, 
                                             levels = c("PB CD16p", "PB CD56p", 
                                                        "PB DN", "BAL NK", 
                                                        "Liver NK", "Spleen NK", 
                                                        "Jejunum NK", "Colon NK"))
  print(sample_prop_sum)

  p <- ggplot(sample_prop_sum, aes(y = sample_prop_sum, x = cell.type.source, fill = class)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    labs(y = "Cumulative Contribution (%)", x  = "cell type & source", title = paste(tissue, "top 10 clones",  sep = " "), 
         subtitle = paste("Top", top_x_clone_num, "clones across all samples of a PB subtype")) + 
    theme_classic() + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + ylim(0, 100) + 
    scale_fill_manual(values = c("PB" = "#F8766D", "Tissue" = "#619CFF", "Sample" = "grey"))
  print(p)
}


