## The purpose of this script is to section out the data relating to 
## particular samples (located in the sample_list.txt file)
## The script will output the counts and metadata files corresponding to the 
## list of samples provided 

library(barcodetrackR)
library(ggplot2)
library(reshape2)

# input information about monkey and samples
monkey_ID <- "ZJ31"
process_date <- 20220512

# load sample information and data 
# file paths
read_dir <- paste("/Volumes/dirgroup/TSCBB/LAB/Barcoded Monkeys/", monkey_ID, "/Processed data/", process_date, sep = "")

# read in metadata file 
metadata_df <- read.delim(file = file.path(read_dir,paste(monkey_ID, "_combined_",process_date,"_metadata.txt", sep="")),
                          sep = "\t")

# read in counts file 
count_data <- read.delim(paste(read_dir, "/", monkey_ID, "_combined_", process_date, "_counts.txt", sep = ''), 
                         header = T, row.names = 1, sep = "\t")
count_data <- count_data[,metadata_df$SAMPLENAME]


# read in sample list 
save_dir <- paste("/Volumes/dirgroup/TSCBB/LAB/Annie/6_Tissue_NK_Manuscript/raw_data", sep = "")

sample_list <- read.delim(file = file.path(save_dir,paste("sample_lists/", monkey_ID, "_samples.txt", sep="")),
                          sep = "\t", header=F)

colnames(sample_list) <- "sample_name"
sample_list <- sample_list$sample_name

sample_index <- rep(NA, length(sample_list))
for (i in 1:length(sample_list)){
  sample_index[i] <- which(metadata_df$SAMPLENAME == sample_list[i])
}

sample_index_unique <- unique(sample_index)

df_sectioned <- count_data[, sample_index_unique] 
metadata_df_sectioned <- metadata_df[sample_index_unique, ]  
my_SE <- barcodetrackR::create_SE(your_data = df_sectioned, meta_data = metadata_df_sectioned, threshold = 0)   
barcode_ggheatmap(my_SE[, c(1:10)], label_size = 2)

write.table(df_sectioned, file = file.path(save_dir,paste("counts_data/", monkey_ID, "_sample_counts.txt", sep = "")), 
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)  

write.table(metadata_df_sectioned, file = file.path(save_dir,paste("counts_data/", monkey_ID, "_sample_metadata.txt", sep="")), 
            sep = "\t", row.names = F, quote = FALSE)
