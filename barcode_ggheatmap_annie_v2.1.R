#' barcode_ggheatmap
#'
#' Creates a heatmap displaying the log abundance of the top 'n' clones from each sample in the Summarized Experiment object, using ggplot2. Clones are on the y-axis and samples are on the x-axis. The ordering and clustering of clones on the y-axis as well as all aesthetics of the plot can be controlled through the arguments described below.
#'
#' @param your_SE A Summarized Experiment object.
#' @param plot_labels Vector of x axis labels. Defaults to colnames(your_SE).
#' @param n_clones The top 'n' clones to plot.
#' @param cellnote_assay Character. One of "stars", "counts", or "proportions." To have no cellnote, set cellnote_size to 0.
#' @param your_title The title for the plot.
#' @param grid Logical. Include a grid or not in the heatmap.
#' @param label_size The size of the column labels.
#' @param dendro Logical. Whether or not to show row dendrogram when hierarchical clustering.
#' @param cellnote_size The numerical size of the cell note labels. To have no cellnote, set cellnote_size to 0.
#' @param distance_method Character. Use summary(proxy::pr_DB) to see all possible options for distance metrics in clustering.
#' @param minkowski_power The power of the Minkowski distance (if minkowski is the distance method used).
#' @param hclust_linkage Character. One of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param row_order Character; "hierarchical" to perform hierarchical clustering on the output and order in that manner, "emergence" to organize rows  by order of presence in data (from left to right), or a character vector of rows within the summarized experiment to plot.
#' @param clusters How many clusters to cut hierarchical tree into for display when row_order is "hierarchical".
#' @param percent_scale A numeric vector through which to spread the color scale (values inclusive from 0 to 1). Must be same length as color_scale.
#' @param color_scale A character vector which indicates the colors of the color scale. Must be same length as percent_scale.
#' @param return_table Logical. Whether or not to return table of barcode sequences with their log abundance in the 'value' column and cellnote for each sample instead of displaying a plot.
#'
#' @return Displays a heatmap in the current plot window. Or if return_table is set to TRUE, returns a dataframe of the barcode sequences, log abundances, and cellnotes for each sample.
#'
#' @importFrom rlang %||%
#' @importFrom stats hclust
#' @importFrom stats cutree
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' data(wu_subset)
#' barcode_ggheatmap(
#'     your_SE = wu_subset, n_clones = 10,
#'     grid = TRUE, label_size = 6
#' )
barcode_ggheatmap_annie_v2.1 <- function(your_SE,
                                       plot_labels = NULL,
                                       n_clones = 10,
                                       cellnote_assay = "stars",
                                       your_title = NULL,
                                       grid = TRUE,
                                       label_size = 12,
                                       dendro = FALSE,
                                       cellnote_size = 4,
                                       distance_method = "Euclidean",
                                       minkowski_power = 2,
                                       hclust_linkage = "complete",
                                       row_order = "hierarchical",
                                       clusters = 0,
                                       clones_of_interest = NULL,
                                       tissue_specific_clones = NULL, 
                                       foldchange_threshold = 10,
                                       percent_scale = c(0, 0.000025, 0.001, 0.01, 0.1, 1),
                                       color_scale = c("#4575B4", "#4575B4", "lightblue", "#fefeb9", "#D73027", "red4"),
                                       return_table = FALSE,
                                       return_clones = FALSE) {
  
  # eliminate new lines from title
  if (!is.null(your_title)) {
    if (length(grep("\n", your_title)) > 0) {
      stop("your_title should not include newline characters")
    }
  }
  #print(foldchange_threshold)
  # get labels for heatmap
  #print(colnames(your_SE))
  if (is.null(plot_labels) == T) {
    plot_labels <- colnames(your_SE)
  }else{
    plot_labels <- plot_labels
  }
  
  if (length(plot_labels) != ncol(your_SE)) {
    stop("plot_labels must be same length as number of columns being plotted")
  }
  
  # error checking
  if (length(percent_scale) != length(color_scale)) {
    stop("percent_scale and color_scale must be vectors of the same length.")
  }
  
  
  # subset the rows of the summarized experiment and get the ordering of barcodes within the heatmap for plotting
  if (row_order == "hierarchical" | row_order == "emergence") {
    
    # subsets those barcodes that have at least one top N clone
    top_clones_choices <- apply(SummarizedExperiment::assays(your_SE)$ranks, 1, function(x) {
      any(x <= n_clones, na.rm = TRUE)
    })
    your_SE <- your_SE[top_clones_choices, ]
    
    # creates data frame with '*' for those cells w/ top clones
    cellnote_matrix <- SummarizedExperiment::assays(your_SE)$ranks
    cellnote_matrix[cellnote_matrix > n_clones] <- NA
    cellnote_matrix[cellnote_matrix <= n_clones] <- "*"
    SummarizedExperiment::assays(your_SE)$stars <- as.data.frame(cellnote_matrix)
    
    # this does the heavy duty plotting set-up. It sets the order of the data on the heatmap and the dendrogram/cluster cuts
    if (row_order == "hierarchical") {
      clustering_data <- SummarizedExperiment::assays(your_SE)[["logs"]]
      clustering_data.dist <- proxy::dist(clustering_data, method = distance_method, p = minkowski_power)
      hclustering <- hclust(clustering_data.dist, method = hclust_linkage)
      barcode_order <- rownames(your_SE)[hclustering$order]
      
      if (dendro) {
        dendro_data <- ggdendro::dendro_data(hclustering, type = "rectangle")
      }
      
      if (clusters > 0) {
        clustercuts_data <- data.frame(
          clusters = cutree(hclustering, clusters),
          assignment = factor(hclustering$labels, levels = hclustering$labels[(hclustering$order)])
        )
      }
    } else if (row_order == "emergence") {
      barcode_order <- rownames(your_SE)[do.call(order, SummarizedExperiment::assays(your_SE)$proportions)]
    }
  } else {
    message("using supplied row_order")
    your_SE <- your_SE[row_order, ]
    barcode_order <- row_order
  }
  
  # set column names as plot_labels
  colnames(your_SE) <- plot_labels
  
  # create scale for plotting
  log_used <- S4Vectors::metadata(your_SE)$log_base
  scale_factor_used <- S4Vectors::metadata(your_SE)$scale_factor
  log_scale <- log(percent_scale * scale_factor_used + 1, base = log_used)
  
  # organizing data for plotting
  plotting_data <- tibble::rownames_to_column(SummarizedExperiment::assays(your_SE)[["logs"]], var = "sequence")
  plotting_data <- tidyr::pivot_longer(plotting_data, cols = -sequence, names_to = "sample_name", values_to = "value")
  plotting_data$sample_name <- factor(plotting_data$sample_name, levels = colnames(your_SE))
  plotting_data$sequence <- factor(plotting_data$sequence, levels = barcode_order)
  
  # organizing labels for plotting overlay
  plotting_cellnote <- tibble::rownames_to_column(SummarizedExperiment::assays(your_SE)[[cellnote_assay]], var = "sequence")
  plotting_cellnote <- tidyr::pivot_longer(plotting_cellnote, cols = -sequence, names_to = "sample_name", values_to = "label")
  plotting_data$cellnote <- plotting_cellnote$label
  if (is.numeric(plotting_data$cellnote)) {
    if (cellnote_assay == "proportions") {
      plotting_data$cellnote <- paste0(round(plotting_data$cellnote * 100, digits = 2), "%")
    } else {
      plotting_data$cellnote <- round(plotting_data$cellnote, digits = 2)
    }
  }
  
  
  if (return_table) {
    return(plotting_data)
  }
  
  if (grid) grid_color <- "black" else grid_color <- NA
  
  # make a plot_label that is invisible -> use it in the dendrogram and cluster bars to make sure they are the same height as the heatmap
  invisible_label <- plot_labels[which(max(nchar(as.character(plot_labels))) == nchar(as.character(plot_labels)))[1]]
  
  
  g1_heatmap <- ggplot2::ggplot(plotting_data, ggplot2::aes(x = .data$sample_name, y = .data$sequence)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$value), color = grid_color) +
    ggplot2::geom_text(ggplot2::aes(label = .data$cellnote), vjust = 0.75, size = cellnote_size, color = "black", na.rm = TRUE) +
    ggplot2::scale_fill_gradientn(
      paste0("Percentage\nContribution"),
      colors = color_scale,
      values = scales::rescale(log_scale, to = c(0, 1)),
      breaks = log_scale,
      limits = c(min(log_scale), max(log_scale)),
      labels = paste0(percent_scale * 100, "%"),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_discrete(labels = NULL, breaks = NULL, expand = c(0, 0)) +
    ggplot2::scale_x_discrete(expand = c(0, 0), labels = plot_labels) +
    ggplot2::ylab(NULL) +
    ggplot2::xlab(NULL) +
    ggplot2::ggtitle(your_title) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = label_size),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = label_size),
      legend.title = ggplot2::element_text(size = label_size),
      legend.key.width = ggplot2::unit(0.2, "cm"),
      legend.text = ggplot2::element_text(size = label_size),
      axis.ticks = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(c(5.5, 5.5, 5.5, 5.5), "pt")
    )
  
  if (row_order != "emergence") {
    if (dendro) {
      g1_heatmap <- g1_heatmap + ggplot2::theme(plot.margin = ggplot2::unit(c(5.5, 5.5, 5.5, 1), "pt"))
      g2_dendrogram <- ggplot2::ggplot(ggdendro::segment(dendro_data)) +
        ggplot2::geom_segment(ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend)) +
        ggplot2::scale_x_discrete(expand = c(.5 / nrow(your_SE), 0.01)) +
        ggplot2::scale_y_reverse(expand = c(0.01, 0), labels = invisible_label, breaks = 1) +
        ggplot2::coord_flip() +
        ggplot2::ylab(NULL) +
        ggplot2::xlab(NULL) +
        ggplot2::theme(
          plot.margin = ggplot2::unit(c(5.5, 0.1, 5.5, 5.5), "pt"),
          plot.title = ggplot2::element_text(size = label_size),
          axis.text.x = ggplot2::element_text(colour = "white", angle = 90, hjust = 1, vjust = 0.5, size = label_size),
          panel.background = ggplot2::element_rect(fill = "white", colour = "white"),
          axis.ticks = ggplot2::element_blank()
        )
      if (!is.null(your_title)) {
        g2_dendrogram <- g2_dendrogram + ggplot2::ggtitle("")
      }
    }
    if (clusters > 0) {
      g1_heatmap <- g1_heatmap + ggplot2::theme(plot.margin = ggplot2::unit(c(5.5, 5.5, 5.5, 1), "pt"))
      g3_clusters <- ggplot2::ggplot(clustercuts_data, ggplot2::aes(x = 1, y = .data$assignment, fill = factor(clusters))) +
        ggplot2::geom_tile() +
        ggplot2::scale_x_continuous(expand = c(0, 0), labels = invisible_label, breaks = 1) +
        ggplot2::scale_y_discrete(expand = c(0, 0)) +
        ggplot2::theme(
          plot.margin = ggplot2::unit(c(5.5, 1, 5.5, 5.5), "pt"),
          plot.title = ggplot2::element_text(size = label_size),
          axis.title = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(colour = "white", angle = 90, hjust = 1, vjust = 0.5, size = label_size),
          legend.position = "none"
        )
      if (dendro) {
        g3_clusters <- g3_clusters + ggplot2::theme(plot.margin = ggplot2::unit(c(5.5, 1, 5.5, 1), "pt"))
      }
      if (!is.null(your_title)) {
        g3_clusters <- g3_clusters + ggplot2::ggtitle("")
      }
    }
  }
  if (!is.null(clones_of_interest)){
    
    num_subsets <- length(clones_of_interest) - 1
    
    for (i in 1:num_subsets){
      
      if (i == 1){
        clones_subset <- i*data.frame(as.integer(barcode_order %in% clones_of_interest[[i]]))
      } else if (i >= 2){
        clones_subset_2 <- i*data.frame(as.integer(barcode_order %in% clones_of_interest[[i]]))
        clones_subset <- clones_subset + clones_subset_2
      }
    }
    clones_of_interest_info <- clones_subset
    #clones_of_interest_info <- data.frame(as.integer(barcode_order %in% clones_of_interest))
    rownames(clones_of_interest_info) <- barcode_order
    #clones_of_interest_info <- clones_of_interest_info[hclustering$labels]
    clones_of_interest_info <- cbind(barcode_order, clones_of_interest_info)
    colnames(clones_of_interest_info) <- c("barcodes", "value")
    
    clones_of_interest_order <- data.frame(assignment = factor(hclustering$labels, levels = hclustering$labels[(hclustering$order)]))
    
    #clones_of_interest_info <- clones_of_interest_info[barcode_order, ]
    clones_of_interest_info <- cbind(clones_of_interest_info, assignment = factor(barcode_order, levels = hclustering$labels[(hclustering$order)]))
    #print(clones_of_interest_info)
    g1_heatmap <- g1_heatmap + ggplot2::theme(plot.margin = ggplot2::unit(c(5.5, 5.5, 5.5, 1), "pt"))
    cols <- c("0" = "grey79", "1" = "#782B9D", "2" = "#de7065ff", "3" = "#339933") #c("0" = "grey79", "1" = "mediumpurple1", "2" = "#7FC97F", "3" = "#FDC086")
    g4_clusters <- ggplot2::ggplot(clones_of_interest_info, ggplot2::aes(x = 1, y = .data$assignment, fill = factor(.data$value))) +
      ggplot2::geom_tile() +
      ggplot2::scale_x_continuous(expand = c(0, 0), labels = invisible_label, breaks = 1) +
      ggplot2::scale_y_discrete(expand = c(0, 0)) +
      ggplot2::theme(
        plot.margin = ggplot2::unit(c(5.5, 1, 5.5, 5.5), "pt"),
        plot.title = ggplot2::element_text(size = label_size),
        axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(colour = "white", angle = 90, hjust = 1, vjust = 0.5, size = label_size),
        legend.position = "bottomright"
      ) + scale_fill_manual(values= cols)#c("grey79", "purple1", "green4", "navyblue"))
    #print("hi")
    cowplot::plot_grid(g4_clusters, g1_heatmap, rel_widths = c(.2, 4), ncol = 2)
    if (dendro) {
      g4_clusters <- g4_clusters + ggplot2::theme(plot.margin = ggplot2::unit(c(5.5, 1, 5.5, 1), "pt"))
    }
    if (!is.null(tissue_specific_clones)) {
      g4_clusters <- g4_clusters + ggplot2::theme(plot.margin = ggplot2::unit(c(5.5, 1, 5.5, 5.5), "pt"))
    }
    if (!is.null(your_title)) {
      g4_clusters <- g4_clusters + ggplot2::ggtitle("")
    }
  }
  
  if (!is.null(tissue_specific_clones)){
    
    #possibleError <- tryCatch(expr = unique(your_SE$source) , error = function(e) e)
    #print(possibleError)
    #print(is.null(possibleError))
    #if(!inherits(possibleError, "error")){
    #print(unique(your_SE$Experiment))
    if (!is.null(unique(your_SE$Experiment))){
      uniq_tissue_list <- unique(your_SE$Experiment)
      tissue_list_all <- your_SE$Experiment
      print(uniq_tissue_list)
    }else if (is.null(unique(your_SE$Experiment))){
      uniq_tissue_list <- unique(your_SE$source)
      tissue_list_all <- your_SE$source
      #print("hello")
    }
    #print(uniq_tissue_list)
    tissue_index_all <- list()
    for (i in 1:length(uniq_tissue_list)){
      tissue_temp <- uniq_tissue_list[i]
      #print(tissue_temp)
      tissue_index <- data.frame(which(tissue_list_all == tissue_temp))
      colnames(tissue_index) <- tissue_temp
      tissue_index_all <- append(tissue_index_all, tissue_index)
    }
    
    #print(tissue_index_all)
    
    #print(tissue_index_all)
    plotting_data_numeric <- plotting_data
    #plotting_data_numeric <- tibble::rownames_to_column(SummarizedExperiment::assays(your_SE)[["logs"]], var = "sequence")
    #print(plotting_data_numeric)
    #print(plotting_data_numeric$value)
    plotting_data_numeric$value <- as.numeric(plotting_data_numeric$value)
    heatmap_data <- reshape2::dcast(plotting_data_numeric, sequence~sample_name)
    rownames(heatmap_data) <- heatmap_data$sequence
    heatmap_data <- heatmap_data[, -1]
    heatmap_data <- heatmap_data[barcode_order, ]
    
    #print(tissue_index_all)
    tissue_mean_clone_all <- data.frame(matrix(NA, nrow = length(barcode_order), ncol = 0))
    for (i in 1:length(tissue_index_all)){
      #print(tissue_index_all[[i]])
      heatmap_tissue_temp <- heatmap_data[, tissue_index_all[[i]]]
      #print(dim(heatmap_tissue_temp))
      if (is.null(dim(heatmap_tissue_temp))){
        tissue_mean_clone <- data.frame(heatmap_tissue_temp)
        colnames(tissue_mean_clone) <- uniq_tissue_list[i]
        rownames(tissue_mean_clone) <- barcode_order
        tissue_mean_clone_all <- cbind(tissue_mean_clone_all, tissue_mean_clone)
      } else{
        tissue_mean_clone <- data.frame(rowMeans(heatmap_tissue_temp))
        colnames(tissue_mean_clone) <- uniq_tissue_list[i]
        rownames(tissue_mean_clone) <- barcode_order
        tissue_mean_clone_all <- cbind(tissue_mean_clone_all, tissue_mean_clone)
      }
    }
    
    #print(tissue_mean_clone_all)
    
    clones_biased_list <- list()
    
    for (i in 1:length(tissue_mean_clone_all)){
      
      tissue_name <- uniq_tissue_list[i]
      print(tissue_name)
      tissue_rowname_index <- which(colnames(tissue_mean_clone_all) == tissue_name)
      #print(tissue_rowname_index)
      heatmap_tissue <- tissue_mean_clone_all[, tissue_rowname_index]
      heatmap_temp_other <- tissue_mean_clone_all[, -tissue_rowname_index]
      
      for (j in 1:length(heatmap_temp_other)){
        #print(j)
        if (j == 1) {
          fold_change <- exp(heatmap_tissue - heatmap_temp_other[, j])
          print(fold_change)
          clones_biased <- which(fold_change >= foldchange_threshold)
          #print(clones_biased)
        } else{
          fold_change <- exp(heatmap_tissue - heatmap_temp_other[, j])
          #print(fold_change)
          clones_biased_temp <- which(fold_change >= foldchange_threshold)
          #print(clones_biased_temp)
          #print(clones_biased)
          clones_biased <- c(clones_biased, clones_biased_temp)
          #print(clones_biased)
        }
        
      }
      # heatmap_temp_other_mean <- rowMeans(heatmap_temp_other)
      # fold_change <- exp(heatmap_tissue - heatmap_temp_other_mean)
      # clones_biased <- which(fold_change > 10)
      clones_biased_table <- table(clones_biased)
      #print(clones_biased_table)
      clones_biased_final <- as.numeric(names(which(clones_biased_table >= length(heatmap_temp_other))))
      #print(clones_biased_final)
      clones_biased_df <- data.frame(barcode_order[clones_biased_final])
      colnames(clones_biased_df) <- tissue_name
      clones_biased_list <- append(clones_biased_list, clones_biased_df)
    }
    
    num_subsets <- length(clones_biased_list)
    #tissue_dict <- c("PB Gr"=1, "PB B"=2, "BAL B"=3, "Liver B"=4, "Spleen B"=5, "Jejunum B"=6, "Colon B"=7)
    print(names(clones_biased_list))
    print("PB CD16p" %in% names(clones_biased_list))
    if ("PB B" %in% names(clones_biased_list) == T){
      tissue_dict <- c("PB Gr"=1, "PB B"=2, "BAL B"=3, "Liver B"=4, "Spleen B"=5, "Jejunum B"=6, "Colon B"=7)
    } else if ("PB CD16p" %in% names(clones_biased_list) == T){
      tissue_dict <- c("PB CD16p"=1, "PB CD56p"=2, "PB DN"=3, "BAL NK"=4, "Liver NK"=5, "Spleen NK"=6, "Jejunum NK"=7, "Colon NK"=8, "PB Gr" = 9)
    } else if ("PB NKG2p_CD16p" %in% names(clones_biased_list) == T){
      tissue_dict <- c("PB NKG2p_CD16p"=1, "PB NKG2p_CD56p"=2, "PB NKG2p_DN"=3, "BAL NKG2p_NK"=4, "Liver NKG2p_NK"=5, "Spleen NKG2p_NK"=6, "Jejunum NKG2p_NK"=7, "Colon NKG2p_NK"=8, "PB Gr" = 9)
    } else if ("PB" %in% names(clones_biased_list) == T){
      tissue_dict <- c("PB"=1, "BAL"=2, "Liver"=3, "Spleen"=4, "Jejunum"=5, "Colon"=6)
    }
    print(tissue_dict)
    #print(names(clones_biased_list))
    for (i in 1:num_subsets){
      
      if (i == 1){
        clones_subset <- tissue_dict[uniq_tissue_list[i]] * tissue_dict[uniq_tissue_list[i]] * data.frame(as.integer(barcode_order %in% clones_biased_list[[i]]))
        #print(clones_subset)
      } else if (i >= 2){
        clones_subset_2 <- tissue_dict[uniq_tissue_list[i]] * tissue_dict[uniq_tissue_list[i]] * data.frame(as.integer(barcode_order %in% clones_biased_list[[i]]))
        clones_subset <- clones_subset + clones_subset_2
        #print(clones_subset)
      }
    }
    
    print(clones_subset)
    clones_of_interest_info <- clones_subset
    #clones_of_interest_info <- data.frame(as.integer(barcode_order %in% clones_biased))
    #print(sum(as.integer(barcode_order %in% names(clones_biased))))
    rownames(clones_of_interest_info) <- barcode_order
    clones_of_interest_info <- cbind(barcode_order, clones_of_interest_info)
    
    colnames(clones_of_interest_info) <- c("barcodes", "value")
    
    clones_of_interest_order <- data.frame(assignment = factor(hclustering$labels, levels = hclustering$labels[(hclustering$order)]))
    
    #clones_of_interest_info <- clones_of_interest_info[barcode_order, ]
    clones_of_interest_info <- cbind(clones_of_interest_info, assignment = factor(barcode_order, levels = hclustering$labels[(hclustering$order)]))
    #print(clones_of_interest_info)
    g1_heatmap <- g1_heatmap + ggplot2::theme(plot.margin = ggplot2::unit(c(5.5, 5.5, 5.5, 1), "pt"))
    g5 <- ggplot2::ggplot(clones_of_interest_info, ggplot2::aes(x = 1, y = .data$assignment, fill = factor(.data$value))) +
      ggplot2::geom_tile() +
      ggplot2::scale_x_continuous(expand = c(0, 0), labels = invisible_label, breaks = 1) +
      ggplot2::scale_y_discrete(expand = c(0, 0)) +
      ggplot2::theme(
        plot.margin = ggplot2::unit(c(5.5, 1, 5.5, 5.5), "pt"), #c(5.5, 1, 5.5, 5.5)
        plot.title = ggplot2::element_text(size = label_size),
        axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(colour = "white", angle = 90, hjust = 1, vjust = 0.5, size = label_size),
        legend.position = "none"
      ) + ggplot2::scale_fill_manual(values= c("0" = "grey79","4" =  "navyblue", "9" = "#8DA0CB", "16" = "#E78AC3", "25" = "#FC8D62", "36" = "#FFD92F", "49" = "#5DC863", "64" = "#E5C494", "1" = "#21908CFF", "81" = "black"))
    #print("hi")
    # colors (0 = grey, 1 = lavendar, 4 = navy blue, 9 = light blue, 16 = light pink, 25 = orange, 36 = yellow; 49 = dark green, 64 = beige)
    #cowplot::plot_grid(g5, g1_heatmap, rel_widths = c(.2, 4), ncol = 2)
    if (dendro) {
      g5 <- g5 + ggplot2::theme(plot.margin = ggplot2::unit(c(5.5, 1, 5.5, 1), "pt"))
    }
    if (!is.null(clones_of_interest)) {
      g5 <- g5 + ggplot2::theme(plot.margin = ggplot2::unit(c(5.5, 3, 5.5, 1), "pt"))
    }
    if (!is.null(your_title)) {
      g5 <- g5 + ggplot2::ggtitle("")
    }
    
    
  }
  if (return_clones){
    return(clones_biased_list)
  }
  # now finally plot using cowplot
  if (row_order == "emergence") {
    g1_heatmap
  } else if (clusters > 0 & dendro ) {
    cowplot::plot_grid(g2_dendrogram, g3_clusters, g1_heatmap, rel_widths = c(1, .2, 4), ncol = 3)
  } else if (clusters == 0 & dendro) {
    cowplot::plot_grid(g2_dendrogram, g1_heatmap, rel_widths = c(1, 4), ncol = 2)
  } else if (clusters > 0 & !dendro) {
    cowplot::plot_grid(g3_clusters, g1_heatmap, rel_widths = c(.2, 4), ncol = 2)
  } else if (clusters == 0 & !dendro & is.null(clones_of_interest) & is.null(tissue_specific_clones)) {
    g1_heatmap
  } else if (clusters == 0 & !dendro & !is.null(clones_of_interest) & is.null(tissue_specific_clones)) {
    cowplot::plot_grid(g4_clusters, g1_heatmap, rel_widths = c(.2, 4), ncol = 2)
  } else if (clusters == 0 & !dendro & is.null(clones_of_interest) & !is.null(tissue_specific_clones)){
    cowplot::plot_grid(g5, g1_heatmap, rel_widths = c(.2, 4), ncol = 2)
  } else if (clusters == 0 & !dendro & !is.null(clones_of_interest) & !is.null(tissue_specific_clones)){
    cowplot::plot_grid(g4_clusters, g5, g1_heatmap, rel_widths = c(.2, .2, 4), ncol = 3)
  }
  
}