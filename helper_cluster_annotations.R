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
