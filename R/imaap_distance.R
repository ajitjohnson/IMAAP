#' @title Average distance between cells of interest
#' @description The fuction plots the distribution for every marker along with the gaussian model fit.
#' @param data Dataframe of data (natural scale) containing samples as rows and markers as columns.
#' @param peaks_low Mean expression value of the left most distribution in the gaussian fit
#' @param peaks_high Mean expression value of the left most distribution in the gaussian fit
#' @param SD Standard Deviation. Default is set to 3. This determines to what extend of the data distribution is considered as positive to marker expression. If your data distribution is clearly not bi-modal, it is suggested to be conservative in the selection of an SD value.
#' @return Plots
#' @export


setwd("/Users/aj/Dropbox (Partners HealthCare)/Data/Vignesh_Lymphoma_tma/cell_type_calling/ALCL_ALK")
load('ALCL.RData')

centroids <- read.csv('alcl_centroids.csv', row.names = 1, header = T)
cell_annotation <- cell_collapsed
cell_annotation$Level.4[cell_annotation$Level.4 %in% c("CD45 Neg Tumor", "Other Tumor")] <- "Tumor cells"
cell_annotation$Level.4[cell_annotation$Level.4 %in% c("Dendritic cells")] <- "CD21+ cells"
cell_annotation$Level.4[cell_annotation$Level.4 %in% c("other")] <- "Unknown cells"

cell_ann_cent <- merge(centroids, cell_annotation, by= "row.names")
cell_ann_cent <- merge(cell_ann_cent, section, by.x= "Row.names", by.y= "universal_ID")
# 40 and 42 are good cores
cell_ann_40 <- cell_ann_cent[cell_ann_cent$core_ID == 40,]
cell_ann_42 <- cell_ann_cent[cell_ann_cent$core_ID == 42,]
row.names(cell_ann_40) <- cell_ann_40[,1]
row.names(cell_ann_42) <- cell_ann_42[,1]



imaap_distance <- function(centroids, cell_annotation) {
  require(plyr)
  require(dplyr)
  require(ggplot2)
  require(gridExtra)
  require(magrittr)
  # subset the centroid dataframe to include only the x and y corrdinates with the cell_id
  centroids_minimal <- as.matrix(centroids[,c('x','y')])
  # function to calculate the euclidian distance
  euc_dist <- function(m) {mtm <- Matrix::tcrossprod(m); sq <- rowSums(m*m);  sqrt(outer(sq,sq,"+") - 2*mtm)}
  # calculate the distance between every pair of cells
  distance_matrix <- euc_dist(centroids_minimal)
  # convert the matrix into a table
  d <- reshape2::melt(distance_matrix)[reshape2::melt(upper.tri(distance_matrix))$value,]
  names(d) <- c("cell-A","cell-B","dist")
  # rename the cell_id's to actual cell types
  d_renamed <- d
  d_renamed$`cell-A` <- mapvalues(d_renamed$`cell-A`, from= row.names(centroids),to=centroids[,cell_annotation])
  d_renamed$`cell-B` <- mapvalues(d_renamed$`cell-B`, from= row.names(centroids),to=centroids[,cell_annotation])
  # Find the combinations of cell types
  unique_combinations <- combn(unique(centroids[,cell_annotation]),2)
  # Intialise empty data frame to hold all the combinations
  u_combinations <- data.frame()
  for (i in 1: (length(unique_combinations)/2)){
    temp <- t(data.frame(unique_combinations[,i]))
    row.names(temp) <- i
    u_combinations <- rbind(u_combinations, temp)
  }
  names(u_combinations) <- c("cell-A","cell-B")
  # Subset all the pairs of cells
  x <- apply(u_combinations, 1, function(x) filter(d_renamed, as.character(d_renamed$`cell-A`) == as.character(x[[1]]) & as.character(d_renamed$`cell-B`) == as.character(x[[2]])))
  # Calculate the mean and median distances for all pairs
  u_combinations$mean_distance <- unlist(lapply(x, function(y) mean(y[,3])))
  u_combinations$median_distance <- unlist(lapply(x, function(y) median(y[,3])))
  # Plotting function
  print("Generating distance distribution plots between all cell types- check Distance distribution plots folder")
  my_plots <- list()
  for (i in 1: length(x)){
    P1 <- ggplot(x[[i]], aes(x = dist)) +  geom_density(fill = "#4271AE")+
      scale_x_continuous(name = "Distance")+ ggtitle(paste(as.character(x[[i]][1,1]),"\nand",as.character(x[[i]][1,2])))+
      geom_vline(xintercept = median(x[[i]][,3]), size = 1, colour = "#000000", linetype = "dashed") + theme_minimal() +
      theme(plot.title = element_text(color="black", size=10))
    my_plots[[i]] <- P1
  }
  dir.create("Distance distribution plots", showWarnings = FALSE)
  pdf(paste("./Distance distribution plots/median distance plot",".pdf" ), width = 9.5, height = 12) # Open a new pdf file
  grid.arrange(grobs = my_plots, ncol = 5) # Write the grid.arrange in the file
  dev.off() #
  return(distance = u_combinations)
  }


distance <- imaap_distance (centroids= cell_ann_40, cell_annotation = "Level.4")

