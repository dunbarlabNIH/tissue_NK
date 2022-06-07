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
sample_list <- read.delim(file = file.path(save_dir,paste(monkey_ID, "/", monkey_ID, "_nk_subset_tissue_reordered.txt", sep="")),
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
metadata_df[["cell_type_source_v2"]] <- paste(metadata_df$source_v2, metadata_df$cell_type)

metadata_df[["Experiment"]] <- metadata_df$source
metadata_df$Experiment <- paste(metadata_df$source, metadata_df$cell_type)

# loading in the data into BarcodetrackR
my_SE <- barcodetrackR::create_SE(your_data = count_data, meta_data = metadata_df, threshold = 0)


#################### TESTING TO MAKE SURE BARCODETRACKR WORKS ####################

barcodetrackR::barcode_ggheatmap(my_SE[,sample_index], n_clones = 10, label_size = 10, 
                                 plot_labels = plot_labels_NK)

######################### CUM CONTRIBUTION OF AGGREGATE TOP X CLONES ######################

top_x_clone_num <- 10
PB_CD16p_samples <- sample_list[1]
PB_CD56p_samples <- sample_list[2]
PB_DN_samples <- sample_list[3]

Liver_CD16p_samples <- sample_list[4:5]
Liver_CD56p_samples <- sample_list[6:7]
Liver_DN_samples <- sample_list[8:9]

Spleen_CD16p_samples <- sample_list[10]
Spleen_CD56p_samples <- sample_list[11]
Spleen_DN_samples <- sample_list[12]

LN_CD16p_samples <- sample_list[13]
LN_CD56p_samples <- sample_list[14]
LN_DN_samples <- sample_list[15]


PB_samples <- list("PB CD16p" = PB_CD16p_samples, "PB CD56p" = PB_CD56p_samples, 
                   "PB DN" = PB_DN_samples, 
                   "Liver CD16p" = Liver_CD16p_samples, "Liver CD56p" = Liver_CD56p_samples, 
                   "Liver DN" = Liver_DN_samples, 
                   "Spleen CD16p" = Spleen_CD16p_samples, "Spleen CD56p" = Spleen_CD56p_samples,
                   "Spleen DN" = Spleen_DN_samples, 
                   "LN CD16p" = LN_CD16p_samples, "LN CD56p" = LN_CD56p_samples,
                   "LN DN" = LN_DN_samples)

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
  sample_prop_sum <- data.frame(sample_prop_sum, source = my_SE$source_v2[sample_index], 
                                cell.type = my_SE$cell_type[sample_index], 
                                cell.type.source = my_SE$cell_type_source[sample_index], 
                                class = factor(c(rep("PB", 3), rep("Tissue", 12)), levels = c("PB", "Tissue", "Sample")))
  print(tissue)
  
  Liver_cd16_samples <- which(sample_prop_sum$source == "Liver" & sample_prop_sum$cell.type == "CD16p")
  Liver_cd56_samples <- which(sample_prop_sum$source == "Liver" & sample_prop_sum$cell.type == "CD56p")
  Liver_DN_samples <- which(sample_prop_sum$source == "Liver" & sample_prop_sum$cell.type == "DN")
  print(sample_prop_sum)
  sample_prop_sum_test <- sample_prop_sum
  sample_prop_sum_test[Liver_cd16_samples[1], ]$sample_prop_sum <- mean(sample_prop_sum[Liver_cd16_samples, ]$sample_prop_sum )
  sample_prop_sum_test[Liver_cd16_samples[1], ]$cell.type.source <- "Liver CD16p"
  sample_prop_sum_test[Liver_cd56_samples[1], ]$sample_prop_sum <- mean(sample_prop_sum[Liver_cd56_samples, ]$sample_prop_sum )
  sample_prop_sum_test[Liver_cd56_samples[1], ]$cell.type.source <- "Liver CD56p"
  sample_prop_sum_test[Liver_DN_samples[1], ]$sample_prop_sum <- mean(sample_prop_sum[Liver_DN_samples, ]$sample_prop_sum )
  sample_prop_sum_test[Liver_DN_samples[1], ]$cell.type.source <- "Liver DN"
  
  sample_prop_sum_test <- sample_prop_sum_test[-c(Liver_cd16_samples[2], Liver_cd56_samples[2], Liver_DN_samples[2]), ]
  
  #sample_prop_sum_test$source[sample_prop_sum_test$cell.type.source == tissue] <- tissue
  sample_prop_sum_test[dim(sample_prop_sum_test)[1] + 1, ] <- sample_prop_sum_test[sample_prop_sum_test$cell.type.source == tissue, ]
  sample_prop_sum_test$source[dim(sample_prop_sum_test)[1]] <- tissue
  sample_prop_sum_test$cell.type[dim(sample_prop_sum_test)[1]] <- "Sample"
  rownames(sample_prop_sum_test)[dim(sample_prop_sum_test)[1]] <- "Sample_being_compared"
  #sample_prop_sum_test$cell.type[sample_prop_sum_test$cell.type.source == tissue] <- "Sample"
  sample_prop_sum_test$cell.type <- factor(sample_prop_sum_test$cell.type, levels = c("Sample", "CD16p", "CD56p", "DN"))
  sample_prop_sum_test$source <- factor(sample_prop_sum_test$source, levels = c(tissue, " ", "PB", "Liver", "Spleen", "LN"))
  sample_prop_sum_test$cell.type.source <- factor(sample_prop_sum_test$cell.type.source,
                                                  levels = c("PB CD16p", "PB CD56p",
                                                             "PB DN", "Liver CD16p",
                                                             "Liver CD56p", "Liver DN",
                                                             "Spleen CD16p", "Spleen CD56p",
                                                             "Spleen DN", "LN CD16p",
                                                             "LN CD56p", "LN DN"))
  sample_prop_sum_test[dim(sample_prop_sum_test)[1] + 1, ]$sample_prop_sum <- 0
  sample_prop_sum_test[dim(sample_prop_sum_test)[1], ]$source <- tissue
  sample_prop_sum_test[dim(sample_prop_sum_test)[1], ]$cell.type <- "CD16p"
  sample_prop_sum_test[dim(sample_prop_sum_test)[1] + 1, ]$sample_prop_sum <- 0
  sample_prop_sum_test[dim(sample_prop_sum_test)[1], ]$source <- tissue
  sample_prop_sum_test[dim(sample_prop_sum_test)[1], ]$cell.type <- "CD56p"
  print(sample_prop_sum_test)
  # #png(file = paste(save_dir, "/", monkey_ID, "/", monkey_ID, "_figure_", figure_num, "_cum_contribution_barplot_", tissue, "_top10_clones_v2.png", sep = ""),
  # #    width=1300, height=900, res = 200)
  # p <- ggplot(sample_prop_sum, aes(y = sample_prop_sum, x = cell.type.source, fill = class)) +
  #   geom_bar(stat = "identity", position = "dodge") +
  #   labs(y = "Cumulative Contribution (%)", x  = "cell type & source", title = paste(tissue, "top 10 clones",  sep = " "),
  #        subtitle = paste("Top", top_x_clone_num, "clones across all samples of a PB subtype")) +
  #   theme_classic() + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + ylim(0, 100) +
  #   scale_fill_manual(values = c("PB" = "#F8766D", "Tissue" = "#619CFF", "Sample" = "purple"))
  # print(p)
  
  p_test <- ggplot(sample_prop_sum_test, aes(y = sample_prop_sum, x = source, fill = cell.type)) + geom_bar(stat = "identity", position = "dodge") + 
    labs(y = "Cumulative Contribution (%)", x  = "Location", 
         title = paste(tissue, "top 10 clones",  sep = " "), 
         subtitle = paste("Top", top_x_clone_num, "clones across all samples of", tissue)) +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + ylim(0, 100) + 
    theme(axis.title = element_text(size = 15), axis.text.x = element_text(size = 15, angle = 90), axis.text = element_text(size = 15), 
          legend.text = element_text(size = 15), legend.title = element_text(size = 15), title = element_text(size = 15))
  print(p_test)
  

  p <- ggplot(sample_prop_sum_test, aes(y = sample_prop_sum, x = source, fill = cell.type)) + geom_bar(stat = "identity", position = "dodge") + 
    labs(y = "Cumulative Contribution (%)", x  = "Location") +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + ylim(0, 100) + 
    theme(axis.title = element_text(size = 25), axis.text.x = element_text(size = 20, angle = 90), axis.text = element_text(size = 20), 
          legend.text = element_text(size = 20), legend.title = element_text(size = 20),  axis.title.y = element_blank())
  print(p)
}

