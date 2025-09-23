library(dplyr)
library(ggplot2)
library(gridExtra)
library(data.table)
library(CELESTA)
library(tidyverse)
library(Rmixmod)
library(spdep)
library(zeallot)


setwd('/Users/jabrand2/Desktop/ff2/')
set.seed(41)

# core 1
f1 <- paste0(getwd(), '/example_data/anal_cancer/summaries/c-7_expression.csv')
marker_exprs1 <- read.table(f1, sep = ",", header = T)
colnames(marker_exprs1) <- gsub('*_coord', '', colnames(marker_exprs1))

marker_exprs1 <- marker_exprs1[,grepl('^X|^Y|*_Cell_Mean', colnames(marker_exprs1))]
marker_exprs1$Y <- marker_exprs1$Y + max(abs(marker_exprs1$Y))
plot(marker_exprs1$X, marker_exprs1$Y, pch = '.')

### The pre-saved imaging data is taken from reg009 of the published CODEX data Schurch et al. Cell,2020

# design matrix
design_matrix <- read.csv(paste0(getwd(), '/celesta/design_matrix.csv'), row.names = 1)

colnames(marker_exprs1) <- gsub("\\.", "", colnames(marker_exprs1))
colnames(design_matrix) <- gsub("\\.", "", colnames(design_matrix)) 

design_matrix <- design_matrix %>%
                    rownames_to_column(var = "celltype")

# If the protein marker is known to be expressed for that cell type, then it 
# is denoted by “1”. If the protein marker is known to not express for a cell type,
# then it is denoted by “0”. If the protein marker is irrelevant or uncertain to
# express for a cell type, then it is denoted by “NA”.

CelestaObj <- CreateCelestaObject(project_title = "core_C-7", 
                                  prior_marker_info = design_matrix,
                                  imaging_data = marker_exprs1)


### Filter out questionable cells. 
### A cell with every marker having expression probability higher than 0.9 are filtered out. 
### And A cell with every marker having expression probability lower than 0.4 are filtered out. 
### User can define the thresholds based on inspecting their data. 
### **This step is optional.** We suggest starting without running this step to see whether there are many doublets/triplets.
# CelestaObj <- FilterCells(CelestaObj, high_marker_threshold = 0.9, low_marker_threshold = 0.4)

### Assign cell types. 
### max_iteration is used to define the maximum iterations allowed in the EM algorithm per round. 
### cell_change_threshold is a user-defined ending condition for the EM algorithm. 
### For example, 0.01 means that when fewer than 1% of the total number of cells do not change identity, the algorithm will stop.

# 11 celltypes
design_matrix[,1]

PlotExpProb(coords=CelestaObj@coords,
            marker_exp_prob = CelestaObj@marker_exp_prob,
            prior_marker_info = design_matrix,
            save_plot = FALSE)

CelestaObj@prior_info$celltype

# CelestaObj <- AssignCells(CelestaObj, 
#                           max_iteration = 10, 
#                           cell_change_threshold = 0.01, 
#                           high_expression_threshold_anchor = c(0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8),
#                           low_expression_threshold_anchor = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
#                           high_expression_threshold_index = c(0.8, 0.8, 0.7, 0.8, 0.9, 1, 0.9, 0.9, 0.9, 0.9, 0.9),
#                           low_expression_threshold_index = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))

# pure defaults
CelestaObj <- AssignCells(CelestaObj)

### Plot cells with CELESTA assigned cell types.
PlotCellsAnyCombination(cell_type_assignment_to_plot = CelestaObj@final_cell_type_assignment[,(CelestaObj@total_rounds+1)],
                        coords = CelestaObj@coords,
                        prior_info = design_matrix,
                        cell_number_to_use=c(1,2,3),
                        cell_type_colors=c("yellow","red","blue"),
                        test_size=1)

### To include unknown cells
PlotCellsAnyCombination(cell_type_assignment_to_plot = CelestaObj@final_cell_type_assignment[,(CelestaObj@total_rounds+1)],
                        coords = CelestaObj@coords,
                        prior_info = design_matrix,
                        cell_type_colors = hcl.colors(palette = 'teal rose', n = 11),
                        cell_number_to_use = 0:11)


# read in labeled data
c7_annos <- read.csv(paste0(getwd(),'/example_data/anal_cancer/temp_annotations/C7_annotations.csv'))
colnames(c7_annos) <- c('id', 'celltype') 

c7_predicted <- data.frame(id = read.csv(f1)$Object.ID, 
                           predicted = CelestaObj@final_cell_type_assignment[,5])

merged_df <- merge(c7_annos, c7_predicted, by = 'id', type = 'left')

table(merged_df$celltype, merged_df$predicted)

plot_df <- data.frame(x = marker_exprs1$X, 
                      y = marker_exprs1$Y, 
                      id =  read.csv(f1)$Object.ID, 
                      celesta_predicted = CelestaObj@final_cell_type_assignment[,5])

plot_df$celesta_predicted <- gsub('epithelia', 'Epi', plot_df$celesta_predicted)
plot_df$celesta_predicted <- gsub('endothelia', 'Endothelia', plot_df$celesta_predicted)
plot_df$celesta_predicted <- gsub('myeloid', 'Myeloid', plot_df$celesta_predicted)
plot_df$celesta_predicted <- gsub('stroma', 'Stroma', plot_df$celesta_predicted)
plot_df$celesta_predicted <- gsub('Bcells', 'xx_Bcell', plot_df$celesta_predicted)
plot_df$celesta_predicted <- gsub('T_cells', 'xx_Tcell', plot_df$celesta_predicted)
plot_df$celesta_predicted <- gsub('NKcells', 'xx_NK', plot_df$celesta_predicted)
plot_df$celesta_predicted <- gsub('immune', 'xx_Immune', plot_df$celesta_predicted)
plot_df$celesta_predicted <- gsub('Unknown', 'xx_Unknown', plot_df$celesta_predicted)

clrs <- tab20_colors <- c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c",
                          "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5",
                          "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f",
                          "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5")

plot_df  %>%
  ggplot(aes(x = x, y = y)) + 
  geom_point(aes(color = celesta_predicted)) + 
  scale_color_manual(values = clrs) +
  theme_classic()

grid.arrange(
merge(plot_df, c7_annos) %>%
  ggplot(aes(x = x, y = y)) +
  geom_point(aes(color = celesta_predicted)) +
  scale_color_manual(values = clrs) +
  theme_classic(), 

merge(plot_df, c7_annos) %>%
  ggplot(aes(x = x, y = y)) +
  geom_point(aes(color = celltype)) +
  scale_color_manual(values = clrs) +
  theme_classic(), 
ncol = 2)



# write.csv(c7_predicted, '/Users/jabrand2/Desktop/fluoroforest/celesta_c7_predicted.csv')

# core 2

f2 <- paste0(getwd(), '/example_data/anal_cancer/summaries/n-12_expression.csv')
marker_exprs2 <- read.table(f2, sep = ",", header = T)
colnames(marker_exprs2) <- gsub('*_coord', '', colnames(marker_exprs2))

marker_exprs2 <- marker_exprs2[,grepl('^X|^Y|*_Cell_Mean', colnames(marker_exprs2))]
marker_exprs2$Y <- marker_exprs2$Y + max(abs(marker_exprs2$Y))


design_matrix <- read.csv(paste0(getwd(), '/celesta/design_matrix.csv'), row.names = 1)

colnames(marker_exprs2) <- gsub("\\.", "", colnames(marker_exprs2))
colnames(design_matrix) <- gsub("\\.", "", colnames(design_matrix)) 

design_matrix <- design_matrix %>%
  rownames_to_column(var = "celltype")

# If the protein marker is known to be expressed for that cell type, then it 
# is denoted by “1”. If the protein marker is known to not express for a cell type,
# then it is denoted by “0”. If the protein marker is irrelevant or uncertain to
# express for a cell type, then it is denoted by “NA”.

CelestaObj <- CreateCelestaObject(project_title = "core_N-12", 
                                  prior_marker_info = design_matrix,
                                  imaging_data = marker_exprs2)


PlotExpProb(coords=CelestaObj@coords,
            marker_exp_prob = CelestaObj@marker_exp_prob,
            prior_marker_info = design_matrix,
            save_plot = FALSE)

CelestaObj@prior_info$celltype

CelestaObj <- AssignCells(CelestaObj)


PlotCellsAnyCombination(cell_type_assignment_to_plot = CelestaObj@final_cell_type_assignment[,(CelestaObj@total_rounds+1)],
                        coords = CelestaObj@coords,
                        prior_info = design_matrix,
                        cell_number_to_use=c(1,2,3),
                        cell_type_colors=c("yellow","red","blue"),
                        test_size=1)

### To include unknown cells
PlotCellsAnyCombination(cell_type_assignment_to_plot = CelestaObj@final_cell_type_assignment[,(CelestaObj@total_rounds+1)],
                        coords = CelestaObj@coords,
                        prior_info = design_matrix,
                        cell_type_colors = hcl.colors(palette = 'teal rose', n = 11),
                        cell_number_to_use = 0:11)


# read in labeled data
n12_annos <- read.csv(paste0(getwd(),'/example_data/anal_cancer/temp_annotations/N12_annotations.csv'))
colnames(n12_annos) <- c('id', 'celltype') 

n12_predicted <- data.frame(id = read.csv(f2)$Object.ID, 
                           predicted = CelestaObj@final_cell_type_assignment[,5])

merged_df <- merge(n12_annos, n12_predicted, by = 'id', type = 'left')

table(merged_df$celltype, merged_df$predicted)

plot_df <- data.frame(x = marker_exprs2$X, 
                      y = marker_exprs2$Y, 
                      id =  read.csv(f2)$Object.ID, 
                      celesta_predicted = CelestaObj@final_cell_type_assignment[,5])


plot_df$celesta_predicted <- gsub('epithelia', 'Epi', plot_df$celesta_predicted)
plot_df$celesta_predicted <- gsub('endothelia', 'Endothelia', plot_df$celesta_predicted)
plot_df$celesta_predicted <- gsub('myeloid', 'Myeloid', plot_df$celesta_predicted)
plot_df$celesta_predicted <- gsub('stroma', 'Stroma', plot_df$celesta_predicted)
plot_df$celesta_predicted <- gsub('T_cells', 'xx_Tcell', plot_df$celesta_predicted)
plot_df$celesta_predicted <- gsub('NKcells', 'xx_NK', plot_df$celesta_predicted)
plot_df$celesta_predicted <- gsub('immune', 'xx_Immune', plot_df$celesta_predicted)
plot_df$celesta_predicted <- gsub('Unknown', 'xx_Unknown', plot_df$celesta_predicted)

clrs <- tab20_colors <- c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c",
                          "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5",
                          "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f",
                          "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5")

plot_df  %>%
  ggplot(aes(x = x, y = y)) + 
  geom_point(aes(color = celesta_predicted)) + 
  theme_classic()

grid.arrange(
  merge(plot_df, n12_annos) %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = celesta_predicted)) +
    scale_color_manual(values = clrs) +
    theme_classic(), 
  
  merge(plot_df, n12_annos) %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = celltype)) +
    scale_color_manual(values = clrs) +
    theme_classic(), 
  ncol = 2)

# write.csv(n12_predicted, '/Users/jabrand2/Desktop/fluoroforest/celesta_n12_predicted.csv')
