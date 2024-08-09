set.seed(2023)

easypackages::libraries("Seurat","tidyverse","SingleR", "ggridges", "viridis", "ggh4x", "celldex", "scCustomize", "SeuratWrappers", "ggExtra", "textTinyR", "patchwork", "pheatmap", "ggrepel", "tidyseurat", "ggpubr", "viridis", "writexl", "readxl", "presto", "scCustomize", "qs", "sctransform","glmGamPoi", "clusterProfiler", "enrichplot", "ExperimentHub", "scuttle", 'scmap',"cowplot",
                        'clusterProfiler',
                        'org.Mm.eg.db',
                        'tidyseurat',
                        'enrichplot',
                        'xlsx', 'clustree', "DESeq2")
obj.v5 = read_rds("../../Documents/machiels_lab_viral/intermediate_data/seurat_obj_central_reclustered.rds")
AM = read_rds("../../Documents/machiels_lab_viral/intermediate_data/seurat_obj_central_am.rds")
obj.v5 |>
  join_features(c("Pparg", "Fabp4", "Itgax", "Mafb", "Ccr2"), assay = "RNA", slot = "scale.data") |> 
  group_by(.feature, day_mock, condition, .cell) |>
  summarise(mean_scaled_dsb = mean(.abundance_RNA)) |>
  group_by(.feature) |> 
  mutate(zscore_mean_scaled_dsb = (mean_scaled_dsb - mean(mean_scaled_dsb)) / sd(mean_scaled_dsb)) |> 
  na.omit() |> 
  mutate(zscore_mean_scaled_dsb = ifelse(zscore_mean_scaled_dsb > 2, 2,
                                         ifelse(zscore_mean_scaled_dsb < -2, -2, zscore_mean_scaled_dsb))) |> 
  ggplot(aes(.feature, day_mock, fill = zscore_mean_scaled_dsb)) +
  geom_tile() +
  scale_fill_viridis() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ggtitle("heatmap_cluster5") +
  ggh4x::facet_nested_wrap(vars(condition, .cell), nrow = 1)


AM$all = paste(AM$day, AM$sample_type, AM$condition, AM$seurat_clusters, sep="_")
ggplot(AM@meta.data, aes(all))+geom_bar(stat="count") +
  theme(axis.text.x = element_text( size = 15, angle = 90), axis.text.y = element_text(angle = 90), axis.title = element_text(size = 18))
ggsave("cellcounts_am_allmetadata.png", bg = "white", width = 50, height = 20, limitsize = FALSE, path = "../../Documents/machiels_lab_viral/Case_study_annotation//")

AM %>%
  as_tibble() %>%z
  # Convert obj.v5 to a tibble
  separate(orig.ident, sep = "__", into = c("day", "sample_type"), remove = FALSE) %>%
  mutate(day = str_replace_all(day, c("viral.experiment.1" = "d60", "viral.experiment.2" = "d8"))) %>%
  mutate(day = ifelse(condition == "Mock", condition, day)) %>%
  mutate(day = factor(day, levels = c("Mock", "d8", "d60")),
         sample_type = as_factor(sample_type),
         condition = as_factor(condition),
         sampletag_Ms4a3 = as_factor(sampletag_Ms4a3),
         seurat_clusters = as_factor(seurat_clusters)
  ) %>%
  dplyr::group_by(day, sample_type, condition, sampletag_Ms4a3, seurat_clusters) %>%
  dplyr::count(day, sample_type, condition, sampletag_Ms4a3, seurat_clusters, .drop = FALSE) %>%
  mutate(sample = paste(day, sample_type, condition, sampletag_Ms4a3)) %>%
  group_by(sample) %>%
  filter(!(day == "Mock" & (condition %in% c("PR8", "MuHV4", "PVM", "MAV1")))) %>%
  mutate(percentage_cells = ifelse(n < 20, NA, n)) %>%
  mutate(condition = as.character(condition)) %>%
  na.omit() %>%
  group_by(seurat_clusters) %>%
  ggplot(aes(condition, seurat_clusters, fill = percentage_cells)) +
  geom_tile() + 
  ggh4x::facet_nested_wrap(
    vars(sampletag_Ms4a3, sample_type, day), nrow = 1, drop = TRUE, scales = "free_x"
  ) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  scale_fill_gradient2(low = "blue", mid = "green", high = "red", limits = c(50, 1000), breaks = seq(50, 1000, by = 50), labels = scales::comma_format(scale = 1)) +
  ggtitle("Proportions of cells from a given sample distributed across clusters",
          "Ms4a3 classification as sorted")
# ggsave(filename = "proportions of cells from a given sample distributed across clusters.png",path = "../../Documents/machiels_lab_viral/Case_study_annotation/exploration",  height = 20, width = 30, units = "cm")

ggsave(filename = "proportions of cells from a given sample distributed across clusters.png",path = "../../Documents/machiels_lab_viral/Case_study_annotation/exploration",  height = 20, width = 30, units = "cm")



AM %>%
  as_tibble() %>%
  # Convert obj.v5 to a tibble
  separate(orig.ident, sep = "__", into = c("day", "sample_type"), remove = FALSE) %>%
  mutate(day = str_replace_all(day, c("viral.experiment.1" = "d60", "viral.experiment.2" = "d8"))) %>%
  mutate(day = ifelse(condition == "Mock", condition, day)) %>%
  mutate(all = as_factor(all)
  ) %>%
  dplyr::group_by(day, sample_type, condition, sampletag_Ms4a3, seurat_clusters, .cell) %>%
  dplyr::count(day, sample_type, condition, sampletag_Ms4a3, seurat_clusters, .drop = FALSE) %>%
  mutate(sample = paste(day, sample_type, condition, sampletag_Ms4a3)) %>%
  group_by(sample) %>%
  filter(!(day == "Mock" & (condition %in% c("PR8", "MuHV4", "PVM", "MAV1")))) %>%
  mutate(percentage_cells = ifelse(n < 20, NA, n)) %>%
  mutate(condition = as.character(condition)) %>%
  na.omit() %>%
  group_by(seurat_clusters) %>%
  ggplot(aes(condition, seurat_clusters, fill = percentage_cells)) +
  geom_tile() + 
  ggh4x::facet_nested_wrap(
    vars(sampletag_Ms4a3, sample_type, day), nrow = 1, drop = TRUE, scales = "free_x"
  ) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  scale_fill_gradient2(low = "blue", mid = "green", high = "red", limits = c(50, 1000), breaks = seq(50, 1000, by = 50), labels = scales::comma_format(scale = 1)) +
  +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 0),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16)
  ) +
  labs(
    title = "Number of Cells in Each Cluster by Sample Type and Sample",
    x = "Condition",
    y = "Cluster",
    fill = "Cell Count"
  )
  ggtitle("Proportions of cells from a given sample distributed across clusters",
          "Ms4a3 classification as sorted")







  
  # Assuming AM is your Seurat object
  metadata <- AM@meta.data
  AM_nomock = AM %>% filter(!str_detect(virus, "Mock"))
  AM_nomock$day_mock =  factor(AM_nomock$day_mock,     levels = c("Mock", "d8","d60"))
  # Create a table of cell counts for each combination of cluster and metadata columns
  Idents(AM_nomock)
  cluster_counts <- table(Idents(AM_nomock), AM_nomock$sample_type, AM_nomock$sampletag_Ms4a3, AM_nomock$day_mock, AM_nomock$virus)
  
  # Convert the table to a data frame for plotting
  cluster_counts_df <- as.data.frame(cluster_counts)
  colnames(cluster_counts_df) <- c("Cluster", "Tissue", "Ms4a3", "Day", "Condition", "Cell_Count")
  write_xlsx(cluster_counts_df, "cluster_Counts_am.xlsx")
  
  # Plot using ggplot2
  plot <- ggplot(cluster_counts_df, aes(Condition, Cluster, fill = Cell_Count)) +
    geom_tile() + 
    ggh4x::facet_nested_wrap(
      vars(Ms4a3, Tissue, Day), nrow = 1, drop = TRUE, scales = "free_x"
    ) + 
    theme_bw() + 
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      strip.text.y = element_text(angle = 0, hjust = 0),
      axis.text.y = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      plot.title = element_text(hjust = 0.5, size = 16)
    ) + 
    scale_fill_gradient2(
      low = "green", mid = "blue", high = "red", 
      limits = c(25, 1000), breaks = seq(25, 1000, by = 25), midpoint = 500,
      labels = scales::comma_format(scale = 1),
      guide = guide_colourbar(barwidth = 5, barheight = 25)
    ) +
    labs(
      title = "Number of Cells in Each Cluster by Sample Type and Sample",
      x = "Condition",
      y = "Cluster",
      fill = "Cell Count"
    ) + 
    ggtitle("Proportions of cells from a given sample distributed across clusters",
            "Ms4a3 classification as sorted")
  plot
  # Save the plot
  ggsave(
    filename = "proportions_of_cells_from_a_given_sample_distributed_across_clusters.png",
    plot = plot,
    path = "../../Documents/machiels_lab_viral/Case_study_annotation/AM",  
    height = 20, width = 30, units = "cm"
  )
  
  # Assuming obj.v5 is another Seurat object
  metadata <- obj.v5@meta.data
  obj.v5_nomock = obj.v5 %>% filter(!str_detect(virus, "Mock"))
  obj.v5_nomock$day_mock =  factor(obj.v5_nomock$day_mock,levels = c("Mock", "d8","d60"))
  # Create a table of cell counts for each combination of cluster and metadata columns
  cluster_counts <- table(Idents(obj.v5_nomock), obj.v5_nomock$sample_type, obj.v5_nomock$sampletag_Ms4a3, obj.v5_nomock$day_mock, obj.v5_nomock$virus)
  
  # Convert the table to a data frame for plotting
  cluster_counts_df <- as.data.frame(cluster_counts)
  colnames(cluster_counts_df) <- c("Cluster", "Tissue", "Ms4a3", "Day", "Condition", "Cell_Count")
  write_xlsx(cluster_counts_df, "cluster_Counts_full.xlsx")
  
  # Plot using ggplot2
  plot <- ggplot(cluster_counts_df, aes(Condition, Cluster, fill = Cell_Count)) +
    geom_tile() + 
    ggh4x::facet_nested_wrap(
      vars(Ms4a3, Tissue, Day), nrow = 1, drop = TRUE, scales = "free_x"
    ) + 
    theme_bw() + 
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      strip.text.y = element_text(angle = 0, hjust = 0),
      axis.text.y = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      plot.title = element_text(hjust = 0.5, size = 16)
    ) + 
    scale_fill_gradient2(
      low = "green", mid = "blue", high = "red", 
      limits = c(25, 1000), breaks = seq(25, 1000, by = 25), midpoint = 500,
      labels = scales::comma_format(scale = 1),
      guide = guide_colourbar(barwidth = 5, barheight = 25)
    ) +
    labs(
      title = "Number of Cells in Each Cluster by Sample Type and Sample",
      x = "Condition",
      y = "Cluster",
      fill = "Cell Count"
    ) + 
    ggtitle("Proportions of cells from a given sample distributed across clusters",
            "Ms4a3 classification as sorted")
  plot
  
  # Save the plot
  ggsave(
    filename = "proportions_of_cells_from_a_given_sample_distributed_across_clusters.png",
    plot = plot,
    path = "../../Documents/machiels_lab_viral/Case_study_annotation/exploration/recluster/",  
    height = 20, width = 30, units = "cm"
  )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Read Excel files and get top 5 genes for each dataset
  MuHV4_Ms4a3 <- read_xlsx("../../Documents/machiels_lab_viral/CaseStudy_DEA/DEA_tables/markers_MuHV4_pos_vs_all_filtered.xlsx")
  top_5_genes_muhv4_ms4a3 <- MuHV4_Ms4a3 %>% arrange(desc(avg_log2FC)) %>% head(5) %>% mutate(source = "top_5_genes_muhv4_ms4a3")
  
  mav1_Ms4a3 <- read_xlsx("../../Documents/machiels_lab_viral/CaseStudy_DEA/DEA_tables/markers_MAV1_pos_vs_all_filtered.xlsx")
  top_5_genes_mav1_ms4a3 <- mav1_Ms4a3 %>% arrange(desc(avg_log2FC)) %>% head(5) %>% mutate(source = "top_5_genes_mav1_ms4a3")
  PR8_Ms4a3 <- read_xlsx("../../Documents/machiels_lab_viral/CaseStudy_DEA/DEA_tables/markers_PR8_pos_vs_all_filtered.xlsx")
  top_5_genes_PR8_ms4a3 <- PR8_Ms4a3 %>% arrange(desc(avg_log2FC)) %>% head(5) %>% mutate(source = "top_5_genes_PR8_ms4a3")
  PVM_Ms4a3 <- read_xlsx("../../Documents/machiels_lab_viral/CaseStudy_DEA/DEA_tables/markers_PVM_pos_vs_all_filtered.xlsx")
  top_5_genes_PVM_ms4a3 <- PVM_Ms4a3 %>% arrange(desc(avg_log2FC)) %>% head(5) %>% mutate(source = "top_5_genes_PVM_ms4a3")
  bottom_5_genes_muhv4_ms4a3 <- MuHV4_Ms4a3 %>% arrange(avg_log2FC) %>% head(5) %>% mutate(source = "bottom_5_genes_muhv4_ms4a3")
  bottom_5_genes_mav1_ms4a3 <- mav1_Ms4a3 %>% arrange(avg_log2FC) %>% head(5) %>% mutate(source = "bottom_5_genes_mav1_ms4a3")
  bottom_5_genes_PR8_ms4a3 <- PR8_Ms4a3 %>% arrange(avg_log2FC) %>% head(5) %>% mutate(source = "bottom_5_genes_PR8_ms4a3")
  bottom_5_genes_PVM_ms4a3 <- PVM_Ms4a3 %>% arrange(avg_log2FC) %>% head(5) %>% mutate(source = "bottom_5_genes_PVM_ms4a3")
  # Assuming you have similar code for bottom 5 genes and other datasets
  # Combine the dataframes
  markers <- bind_rows(
    top_5_genes_PVM_ms4a3,
    bottom_5_genes_PVM_ms4a3,
    top_5_genes_PR8_ms4a3,
    bottom_5_genes_PR8_ms4a3,
    top_5_genes_mav1_ms4a3,
    bottom_5_genes_mav1_ms4a3,
    top_5_genes_muhv4_ms4a3,
    bottom_5_genes_muhv4_ms4a3
  )
  genes = markers$gene
  genes = unique(genes)
  seurat_d60 = AM %>% filter(str_detect(day, "d60"))
  pallet = viridis(n = 10, option = "D")
  Idents(seurat_d60) =  "sampletag_name"
  DefaultAssay(seurat_d60) = "RNA"
  dotplot_celltypes <- DotPlot_scCustom(seurat_object = seurat_d60, features = genes, colors_use = c("green", "red"), group.by = "sampletag_name", flip_axes = FALSE)
  
  # Customize the plot
  dotplot_celltypes <- dotplot_celltypes +
    theme(
      axis.text.x = element_text(size = 20, angle = 90),
      axis.text.y = element_text(size = 25),
      legend.text = element_text(size = 15),  # Increase legend text size
      legend.title = element_text(size = 20),  # Increase legend title size
      axis.title.x = element_text(size = 20),  # Increase x-axis title size
      axis.title.y = element_text(size = 25)   # Increase y-axis title size
    ) +
    guides(
      size = guide_legend(override.aes = list(size = 15))  # Increase legend key size
    ) +
    scale_size(range = c(5, 15)) 

  ggsave("dotplot_sampletagname.png", width = 40, height = 20, bg = 'white',  path =  "../../Documents/machiels_lab_viral/Case_study_annotation/exploration/recluster/")
  
  
  
  
  
  
  
  
  
  
  
  
  umap_plot <- DimPlot(obj.v5, reduction = "umap_harmony",  label = TRUE, pt.size = 0.5) + NoLegend()
  
  # Get the number of cells per cluster
  cell_counts <- table(Idents(obj.v5))
  cell_counts <- as.data.frame(cell_counts)
  colnames(cell_counts) <- c("cluster", "count")
  
  # Extract UMAP coordinates for cluster labels
  umap_coords <- Embeddings(obj.v5, "umap_harmony")
  umap_coords <- as.data.frame(umap_coords)
  umap_coords$cluster <- Idents(obj.v5)
  
  # Calculate the center of each cluster
  cluster_centers <- aggregate(umap_coords[,1:2], list(umap_coords$cluster), mean)
  colnames(cluster_centers) <- c("cluster", "UMAP_1", "UMAP_2")
  cluster_centers <- merge(cluster_centers, cell_counts, by="cluster")
  
  # Add cluster size to the plot
  umap_plot <- umap_plot + 
    geom_text(data = cluster_centers, aes(x = UMAP_1, y = UMAP_2, label = count), color = "black", size = 5)
  
  # Display the plot
  print(umap_plot)
  ggsave("umap_counts_fullobj.png", path = "../../Documents/machiels_lab_viral/Case_study_annotation/exploration/recluster/", bg = "white")
  
  
 
  
  
  
  
  
  
  
  
  
  
  
  
   
  seurat_full_d60 = obj.v5 %>% filter(str_detect(day, "d60"))
  seurat_full_d8 = obj.v5 %>% filter(str_detect(day, "d8"))  
  seurat_full_d60_bal = seurat_full_d60 %>% filter(str_detect(sample_type, "bal"))
  seurat_full_d60_lung = seurat_full_d60 %>% filter(str_detect(sample_type, "lung"))
  seurat_full_d8_bal = seurat_full_d8 %>% filter(str_detect(sample_type, "bal"))
  seurat_full_d8_lung = seurat_full_d8 %>%  filter(str_detect(sample_type, "lung"))
  
listSeuratobj = list (obj.v5, seurat_full_d60_bal, seurat_full_d60_lung, seurat_full_d8_bal, seurat_full_d8_lung)



for (i in seq_along(listSeuratobj)) {
  obj <- listSeuratobj[[i]]
  obj = obj %>% filter(!str_detect(condition, "Mock"))
  
  # Create UMAP plot
  umap_plot <- DimPlot(obj, reduction = "umap_harmony", split.by = "condition" ,label = TRUE, pt.size = 0.5) + NoLegend()
  
  # Get the number of cells per cluster
  cell_counts <- table(Idents(obj))
  cell_counts <- as.data.frame(cell_counts)
  colnames(cell_counts) <- c("cluster", "count")
  
  # Extract UMAP coordinates for cluster labels
  umap_coords <- Embeddings(obj, "umap_harmony")
  umap_coords <- as.data.frame(umap_coords)
  umap_coords$cluster <- Idents(obj)
  
  # Calculate the center of each cluster
  cluster_centers <- aggregate(umap_coords[,1:2], list(umap_coords$cluster), mean)
  colnames(cluster_centers) <- c("cluster", "UMAP_1", "UMAP_2")
  cluster_centers <- merge(cluster_centers, cell_counts, by="cluster")
  
  # Add cluster size to the plot
  umap_plot <- umap_plot + 
    geom_text(data = cluster_centers, aes(x = UMAP_1, y = UMAP_2, label = count), color = "black", size = 5)
  
  # Display the plot
  print(umap_plot)
  
  # Save the plot with a dynamic file name
  ggsave(filename = paste0("umap_counts_obj_", i, ".png"), plot = umap_plot, 
         path = "../../Documents/machiels_lab_viral/Case_study_annotation/exploration/recluster/", 
         bg = "white",
         width = 20, height = 10)
}





seurat_AM_d60 = AM %>% filter(str_detect(day, "d60"))
seurat_AM_d8 = AM %>% filter(str_detect(day, "d8"))  
seurat_AM_d60_bal = seurat_AM_d60 %>% filter(str_detect(sample_type, "bal"))
seurat_AM_d60_lung = seurat_AM_d60 %>% filter(str_detect(sample_type, "lung"))
seurat_AM_d8_bal = seurat_AM_d8 %>% filter(str_detect(sample_type, "bal"))
seurat_AM_d8_lung = seurat_AM_d8 %>%  filter(str_detect(sample_type, "lung"))

listSeuratobj = list (AM, seurat_AM_d60_bal, seurat_AM_d60_lung, seurat_AM_d8_bal, seurat_AM_d8_lung)
# Load necessary libraries
library(Seurat)
library(ggplot2)

# Assuming listSeuratobj is already defined
for (i in seq_along(listSeuratobj)) {
  obj <- listSeuratobj[[i]]
  obj = obj %>% filter(!str_detect(condition, "Mock"))
  
  # Create UMAP plot
  umap_plot <- DimPlot(obj, reduction = "umap_harmony", split.by = "condition" ,label = TRUE, label.size = 8, label.color = "blue",  pt.size = 0.5) + NoLegend()
  
  # Get the number of cells per cluster
  cell_counts <- table(Idents(obj))
  cell_counts <- as.data.frame(cell_counts)
  colnames(cell_counts) <- c("cluster", "count")
  
  # Extract UMAP coordinates for cluster labels
  umap_coords <- Embeddings(obj, "umap_harmony")
  umap_coords <- as.data.frame(umap_coords)
  umap_coords$cluster <- Idents(obj)
  
  # Calculate the center of each cluster
  cluster_centers <- aggregate(umap_coords[,1:2], list(umap_coords$cluster), mean)
  colnames(cluster_centers) <- c("cluster", "UMAP_1", "UMAP_2")
  cluster_centers <- merge(cluster_centers, cell_counts, by="cluster")
  
  # Add cluster size to the plot
  umap_plot <- umap_plot + 
    geom_text_repel(data = cluster_centers, aes(x = UMAP_1, y = UMAP_2, label = count), color = "black", size = 5, point.padding = 8)
  
  
  # Display the plot
  print(umap_plot)
  
  # Save the plot with a dynamic file name
  ggsave(filename = paste0("umap_counts_obj_", i, ".png"), plot = umap_plot, 
         path = "../../Documents/machiels_lab_viral/Case_study_annotation/AM/", 
         bg = "white",
         width = 20, height = 10)
}






cluster5 = subset(obj.v5, idents = c("5"))
cluster5 = cluster5 %>% filter(str_detect(day, "d60"))
cluster5 |>
  join_features(c("Pparg", "Fabp4", "Itgax", "Mafb", "Ccr2"),assay="RNA", slot="scale.data")  |> 
  group_by(.feature,day, condition, .cell) |>
  summarise(mean_scaled_dsb=mean(.abundance_RNA)) |>
  group_by(.feature) |> 
  mutate(zscore_mean_scaled_dsb=(mean_scaled_dsb-mean(mean_scaled_dsb))/sd(mean_scaled_dsb))|> na.omit() |> 
  mutate(zscore_mean_scaled_dsb=ifelse(zscore_mean_scaled_dsb>2,2,
                                       ifelse(zscore_mean_scaled_dsb<(-2),-2,zscore_mean_scaled_dsb))) |> 
  ggplot(aes(.feature,day, fill=zscore_mean_scaled_dsb)) +
  geom_tile()+ scale_fill_viridis() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +ggtitle("heatmap_cluster5") +   ggh4x::facet_nested_wrap(
    vars(condition), nrow = 1 )
ggsave("hetmap_markers_cluster5.png", path  = "../../Documents/machiels_lab_viral/Case_study_annotation/exploration/recluster")


# Load necessary libraries
library(ggplot2)
library(reshape2)

# Create the data frame
lung_d8 <- data.frame(
  Sample = c("Mock_Ms4a3_neg", "Mock_Ms4a3_pos", "MuHV4_Ms4a3_neg", "MuHV4_Ms4a3_pos",
             "MAV1_Ms4a3_neg", "MAV1_Ms4a3_pos", "PR8_Ms4a3_neg", "PR8_Ms4a3_pos",
             "PVM_Ms4a3_neg", "PVM_Ms4a3_pos"),
  lungd8_0 = c(42177759.26, 39252722.65, 22586086.96, 5664919.21, 24305505.38, 
         13562512.87, 86594355.83, 25799869.35, 68059675.68, 7176017.57),
  lungd8_1 = c(1421722.222, 2156743.003, 32206086.96, 442241359.7, 31249935.48, 
         486555149.3, 52317423.31, 771959248.8, 25522378.38, 531025300.1),
  lungd8_2 = c(0, 1294045.802, 9620000, 27569273.49, 2777772.043, 41252643.32, 
         6314171.779, 17652542.19, 15799567.57, 38272093.7),
  lungd8_3 = c(3791259.259, 39684071.25, 10874782.61, 145399593.1, 9027759.14, 
         181963714.4, 8118220.859, 133073010.3, 19445621.62, 219267203.5),
  lungd8_4_0 = c(0, 431348.6005, 7528695.652, 10574515.86, 8333316.129, 
           7346361.14, 13530368.1, 5431551.443, 15799567.57, 8770688.141),
  lungd8_4_1 = c(4265166.667, 1294045.802, 10874782.61, 4909596.649, 2777772.043, 
           2825523.515, 4510122.699, 2715775.721, 8507459.459, 3189341.142),
  lungd8_4_2 = c(0, 0, 18821739.13, 14351128.67, 10416645.16, 
           9041675.249, 11726319.02, 8826271.094, 20660972.97, 26312064.42),
  lungd8_4_3 = c(947814.8148, 4744834.606, 6692173.913, 5664919.21, 2777772.043, 
           3955732.921, 1804049.08, 1357887.861, 6076756.757, 2392005.857),
  lungd8_4_4 = c(947814.8148, 862697.201, 1254782.609, 3776612.807, 0, 
           1130209.406, 902024.5399, 3394719.652, 3646054.054, 1594670.571),
  lungd8_5 = c(947814.8148, 23724173.03, 2927826.087, 71377982.05, 0, 
         79114658.43, 1804049.08, 14257822.54, 0, 40664099.56),
  lungd8_6 = c(1895629.63, 2156743.003, 418260.8696, 2643628.965, 2083329.032, 
         6781256.437, 3608098.16, 6110495.373, 3646054.054, 5581346.999),
  lungd8_7 = c(0, 0, 1254782.609, 3021290.245, 5555544.086, 
         9606779.952, 902024.5399, 3394719.652, 4861405.405, 9568023.426),
  lungd8_8_0 = c(1895629.63, 862697.201, 418260.8696, 0, 0, 
           565104.7031, 902024.5399, 678943.9303, 0, 797335.2855),
  lungd8_8_1 = c(473907.4074, 0, 1673043.478, 1888306.403, 2777772.043, 
           3390628.218, 902024.5399, 1357887.861, 0, 2392005.857),
  lungd8_8_2 = c(0, 0, 0, 0, 0, 
           0, 0, 0, 0, 0)
)

# Melt the data for ggplot2
melted_data_lung8 <- melt(lung_d8, id.vars = "Sample")

# Create the plot
ggplot(melted_data_lung8, aes(x = variable, y = value, fill = Sample)) + 
  geom_bar(stat = "identity", position = "dodge") +

  labs(title = "Data Visualization with Logarithmic Scale", x = "Clusters", y = "amount cells") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_brewer(palette = "Paired", name = "Samples")
ggsave("barplot_d8_lung.png", bg = "white", height = 10, width = 20)



data <- data.frame(
  Sample = c("Mock_Ms4a3_neg", "Mock_Ms4a3_pos", "MuHV4_Ms4a3_neg", "MuHV4_Ms4a3_pos",
             "MAV1_Ms4a3_neg", "MAV1_Ms4a3_pos", "PR8_Ms4a3_neg", "PR8_Ms4a3_pos",
             "PVM_Ms4a3_neg", "PVM_Ms4a3_pos"),
  lungd60_0 = c(23639491.95, 12457425.47, 3433238.095, 3935454.545, 19003820.47, 
         7422295.877, 26351629.93, 8684525.229, 12047211.71, 6818339.286),
  lungd60_1 = c(2183265.18, 2676802.168, 4364285.714, 32270727.27, 3283382.594, 
         13702700.08, 7411395.918, 22837084.86, 3073268.293, 53446982.14),
  lungd60_2 = c(75285.0062, 514769.6477, 2153047.619, 48602863.64, 596978.6535, 
         18384455.94, 720552.381, 7397928.899, 1690297.561, 37830785.71),
  lungd60_3 = c(978705.0805, 29650731.71, 1396571.429, 48406090.91, 1691439.518, 
         59264177.85, 2779273.469, 60148378.44, 983445.8537, 52127303.57),
  lungd60_4_0 = c(1430415.118, 102953.9295, 1803904.762, 590318.1818, 1790935.961, 
           342567.502, 2676337.415, 321649.0826, 1628832.195, 879785.7143),
  lungd60_4_1 = c(3011400.248, 720677.5068, 3898761.905, 1180636.364, 4477339.901, 
           570945.8367, 6484971.429, 643298.1651, 4241110.244, 439892.8571),
  lungd60_4_2 = c(752850.062, 102953.9295, 1222000, 1180636.364, 497482.2113, 
           0, 1132296.599, 0, 983445.8537, 1099732.143),
  lungd60_4_3 = c(2559690.211, 1235447.154, 3549619.048, 1180636.364, 3979857.69, 
           1370270.008, 2470465.306, 643298.1651, 3595723.902, 1979517.857),
  lungd60_4_4 = c(1580985.13, 102953.9295, 698285.7143, 590318.1818, 2785900.383, 
           570945.8367, 3191017.687, 160824.5413, 1598099.512, 219946.4286),
  lungd60_5 = c(301140.0248, 21002601.63, 232761.9048, 59622136.36, 497482.2113, 
         33799993.53, 720552.381, 44548397.94, 430257.561, 40690089.29),
  lungd60_6 = c(2559690.211, 1235447.154, 523714.2857, 3935454.545, 1392950.192, 
         2055405.012, 1441104.762, 964947.2477, 1014178.537, 2859303.571),
  lungd60_7 = c(150570.0124, 102953.9295, 116380.9524, 1967727.273, 497482.2113, 
         685135.004, 308808.1633, 482473.6239, 737584.3902, 1759571.429),
  lungd60_8_0 = c(3914820.322, 1029539.295, 3084095.238, 393545.4545, 5273311.44, 
           1027702.506, 4529186.395, 964947.2477, 5808477.073, 1979517.857),
  lungd60_8_1 = c(1053990.087, 411815.7182, 1338380.952, 590318.1818, 596978.6535, 
           0, 823488.4354, 0, 860515.122, 0),
  lungd60_8_2 = c(301140.0248, 102953.9295, 58190.47619, 590318.1818, 298489.3268, 
           228378.3347, 1132296.599, 643298.1651, 276594.1463, 439892.8571)
)

# Melt the data for ggplot2
lungd60 <- melt(data, id.vars = "Sample")

# Create the plot
ggplot(lungd60, aes(x = variable, y = value, fill = Sample)) + 
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Data Visualization with Logarithmic Scale", x = "Clusters", y = "amount cells") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_brewer(palette = "Paired", name = "Samples")
ggsave("barplot_d60_lung.png", bg = "white", height = 10, width = 15)


data <- data.frame(
  Sample = c("Mock_Ms4a3_neg", "Mock_Ms4a3_pos", "MuHV4_Ms4a3_neg", "MuHV4_Ms4a3_pos",
             "MAV1_Ms4a3_neg", "MAV1_Ms4a3_pos", "PR8_Ms4a3_neg", "PR8_Ms4a3_pos",
             "PVM_Ms4a3_neg", "PVM_Ms4a3_pos"),
  bald8_0 = c(2034166.641, 1349584.615, 54766.98001, 70841.56206, 88868.75246, 
         132782.623, 1189795.631, 736515.4059, 77381.2504, 63943.56707),
  bald8_1 = c(8008.530084, 47353.84615, 293750.1655, 727058.1369, 370286.4686, 
         2250981.609, 945317.0769, 2622842.633, 50121.94628, 1157786.714),
  bald8_2 = c(88093.83092, 78923.07692, 4978.816364, 0, 59245.83497, 
         31614.91024, 81492.85146, 79480.07977, 13189.98586, 13605.01427),
  bald8_3 = c(0, 47353.84615, 0, 41013.53593, 22217.18811, 
         44260.87433, 81492.85146, 243738.9113, 3517.329563, 29931.0314),
  bald8_4_0 = c(0, 15784.61538, 144385.6746, 78298.56859, 185143.2343, 
           145428.5871, 334120.691, 116570.7837, 84415.90952, 43536.04567),
  bald8_4_1 = c(8008.530084, 7892.307692, 0, 1864.251633, 0, 
           0, 40746.42573, 5298.671985, 5275.994345, 1360.501427),
  bald8_4_2 = c(0, 0, 129449.2255, 42877.78756, 59245.83497, 
           69552.80252, 195582.8435, 100674.7677, 69467.25888, 36733.53853),
  bald8_4_3 = c(16017.06017, 15784.61538, 4978.816364, 0, 0, 
           0, 0, 0, 879.3323909, 0),
  bald8_4_4 = c(0, 0, 4978.816364, 11185.5098, 0, 
           12645.96409, 48895.71087, 5298.671985, 2637.997173, 1360.501427),
  bald8_5 = c(0, 173630.7692, 0, 54063.29736, 0, 
         183366.4794, 0, 26493.35992, 0, 28570.52997),
  bald8_6 = c(208221.7822, 126276.9231, 29872.89819, 16778.2647, 74057.29371, 
         50583.85638, 179284.2732, 100674.7677, 53639.27584, 50338.5528),
  bald8_7 = c(8008.530084, 0, 54766.98001, 31692.27776, 548023.9735, 
         505838.5638, 839376.37, 471581.8067, 365802.2746, 340125.3568),
  bald8_8_0 = c(0, 7892.307692, 0, 1864.251633, 0, 
           0, 0, 5298.671985, 879.3323909, 0),
  bald8_8_1 = c(0, 0, 0, 3728.503266, 0, 
           6322.982047, 8149.285146, 10597.34397, 0, 0),
  bald8_8_2 = c(0, 0, 0, 0, 0, 
           0, 0, 0, 0, 1360.501427)
)

# Melt the data for ggplot2
bald8 <- melt(data, id.vars = "Sample")

# Create the plot
ggplot(bald8, aes(x = variable, y = value, fill = Sample)) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_log10() +
  labs(title = "Data Visualization with Logarithmic Scale", x = "Clusters", y = "amount cells") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_brewer(palette = "Paired", name = "Samples")
ggsave("barplot_d8_bal.png", bg = "white", height = 10, width = 15)






data <- data.frame(
  Sample = c("Mock_Ms4a3_neg", "Mock_Ms4a3_pos", "MuHV4_Ms4a3_neg", "MuHV4_Ms4a3_pos",
             "MAV1_Ms4a3_neg", "MAV1_Ms4a3_pos", "PR8_Ms4a3_neg", "PR8_Ms4a3_pos",
             "PVM_Ms4a3_neg", "PVM_Ms4a3_pos"),
  bal_60_0 = c(2968315.433, 1700041.13, 712800, 554400, 2707281.553, 
         1122634.111, 2945939.976, 1391978.369, 2170981.561, 948484.547),
  bal_60_1 = c(14305.13462, 4266.100703, 3771.428571, 93600, 13349.51456, 
         78388.41675, 13495.79832, 21415.05183, 9623.145216, 61372.52951),
  bal_60_2 = c(21457.70193, 51193.20843, 346028.5714, 8582400, 165533.9806, 
         3762644.004, 46271.30852, 1234934.655, 423418.3895, 5936397.4),
  bal_60_3 = c(0, 2133.050351, 0, 21600, 0, 
         0, 0, 4758.900406, 0, 0),
  bal_60_4_0 = c(0, 0, 5657.142857, 0, 0, 
           0, 5783.913565, 2379.450203, 0, 0),
  bal_60_4_1 = c(3576.283655, 0, 16028.57143, 0, 13349.51456, 
           0, 5783.913565, 0, 5773.88713, 0),
  bal_60_4_2 = c(0, 0, 0, 0, 0, 
           0, 1927.971188, 0, 0, 0),
  bal_60_4_3 = c(0, 0, 1885.714286, 0, 0, 
           0, 0, 0, 3849.258087, 0),
  bal_60_4_4 = c(0, 0, 942.8571429, 0, 0, 
           0, 0, 0, 0, 0),
  bal_60_5 = c(0, 19197.45316, 0, 57600, 0, 
         33595.03575, 0, 21415.05183, 0, 72531.17124),
  bal_60_6 = c(239611.0049, 187708.4309, 81085.71429, 849600, 240291.2621, 
         562716.8488, 237140.4562, 240324.4705, 196312.1624, 719732.3916),
  bal_60_7 = c(17881.41827, 6399.151054, 35828.57143, 244800, 37378.64078, 
         75588.83043, 5783.913565, 2379.450203, 1924.629043, 11158.64173),
  bal_60_8_0 = c(0, 0, 942.8571429, 0, 8009.708738, 
           0, 0, 0, 0, 0),
  bal_60_8_1 = c(0, 0, 942.8571429, 0, 0, 
           0, 1927.971188, 0, 0, 0),
  bal_60_8_2 = c(0, 0, 0, 0, 0, 
           0, 0, 0, 0, 0)
)

# Melt the data for ggplot2
bald60 <- melt(data, id.vars = "Sample")

# Create the plot
ggplot(bald60, aes(x = variable, y = value, fill = Sample)) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_log10() +
  labs(title = "Data Visualization with Logarithmic Scale", x = "Clusters", y = "amount cells") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_brewer(palette = "Paired", name = "Samples")
ggsave("barplot_d60_bal.png", bg = "white", height = 10, width = 15)