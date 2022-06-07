library(barcodetrackR)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(SummarizedExperiment)

# read in the modified heatmap function 
source("/Volumes/dirgroup/TSCBB/LAB/Annie/5_Tissue_NK_Study/Figure4/ZJ31/barcode_ggheatmap_annie_v2.1.R")

# figure number 
figure_num <- 5

# input information about monkey and samples
monkey_ID <- "ZJ31"

######### INITIAL SETUP (Loading in files + locating samples) #################

data_dir <- paste("/Volumes/dirgroup/TSCBB/LAB/Annie/6_Tissue_NK_Manuscript/raw_data", sep = "")
save_dir <- paste("/Volumes/dirgroup/TSCBB/LAB/Annie/6_Tissue_NK_Manuscript/", "Figure", figure_num, sep = "")

# read in sample list 
sample_list <- read.delim(file = file.path(save_dir,paste(monkey_ID, "/", monkey_ID, "_nk_subset_tissue.txt", sep="")),
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


#################HEATMAPS + CORRELATION ########################

# Top 10 clones heatmap 
barcodetrackR::barcode_ggheatmap(my_SE[,sample_index], n_clones = 10, label_size = 10, 
                                 plot_labels = plot_labels_NK) # regular heatmap 


# Correlation plot (negatives turned off)
cor_plot(my_SE[,sample_index], method_corr = "pearson", 
         plot_type = "color", assay = "proportions", 
         plot_labels = plot_labels_NK, label_size = 25, no_negatives = T)


