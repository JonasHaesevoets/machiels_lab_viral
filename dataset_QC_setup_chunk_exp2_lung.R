#gets called by QC_dashboard_layout.qmd
set.seed(2023)
#set up file paths
file_path <- vector("list") 
file_path$output <- ".\\output\\" 
file_path$intermediate_data<- ".\\intermediate_data\\" 

obj_identifier <- "experiment_2_lung"
file_name_obj <- paste0("seurat_obj_", obj_identifier, "_workflowed.rds")
#file_name_obj <- "seurat_obj_experiment_1_bal_workflowed.rds"


seurat_obj <- read_rds(file = paste0(file_path$intermediate_data,file_name_obj))






de_genes_tbl <- read_csv(
  "intermediate_data/experiment_2_workflowed.Lung_QCmarkers_min.pct_0.4_logfc.threshold_0.25_max.cells.per.ident_300.csv")|>
  relocate(gene) 
file_path_output <- paste0(file_path$intermediate_data,"seurat_obj_experiment_2_lung_afterQCdashboard.rds")


#visualisation parameters
protein_all_cells_sum_tresh <- 100
tags_all_cells_sum_thresh <- 1000
xlim_max <- 70000
ylim_max <- 6
ylim_max_mito <- 6


#cutoff quantiles
upper_mito_thresh_quantile <- 0.98
lower_nCount_RNA_thresh_quantile <- 0.12
upper_nCount_RNA_thresh_quantile <- 0.97

# cut off absolut values
upper_mito_thresh <- seurat_obj$percent_mito |>
        quantile(probs = upper_mito_thresh_quantile)

lower_nCount_RNA_thresh <- seurat_obj$nCount_RNA |>
        quantile(probs = lower_nCount_RNA_thresh_quantile)
upper_nCount_RNA_thresh <- seurat_obj$nCount_RNA |>
        quantile(probs = upper_nCount_RNA_thresh_quantile)
upper_nCount_RNA_thresh_non_scientific <- format(round(upper_nCount_RNA_thresh,digits = 4), scientific = FALSE)

# clusters
deleted_clusters=c(3,6, 8, 9, 10, 14)

# cluster_determination <- c("low_RNA_5"=5,
#                            "low_RNA_mt_high_10"=10,
#                            "contamination_B_NK"=13)


seurat_obj <- seurat_obj |> mutate(orig_cluster_specification=case_when(
         seurat_clusters==3 ~ "3_T_cell",
         seurat_clusters==6 ~ "6_B_cells",
         seurat_clusters==8 ~ "8_T_cell_prolif.",

         seurat_clusters==9 ~ "9_epithelial_cell",
        seurat_clusters==10 ~ "10_NK_cell",
        seurat_clusters==14 ~ "14_basophil_mast_cell",
        TRUE~seurat_clusters

 ))
 del_clusters_string <- paste0(as.character(deleted_clusters))




#genes to gighlight in feature umaps
clus_0_genes <-  c(de_genes_tbl |> filter(cluster==0) |> slice_max(avg_log2FC,n = 5) |> pull(gene), c("H2-Ab1"))
title_feature_umap_clus_0 <- "Macrophage, MHC complex"
caption_feature_umap_clus_0 <-""  # 

clus_1_genes <- c(c("B2m", "Ccl2", "Ctsb", "Cstb"), de_genes_tbl |> filter(cluster==1) |> slice_max(avg_log2FC,n = 5) |> pull()) |> unique()
title_feature_umap_clus_1 <- "Macrophage - cathepsin, cystatin"
caption_feature_umap_clus_1 <- "" # 

clus_2_genes <- c(c("Isg15", "Gbp2", "Irf", "Naaa") , de_genes_tbl |> filter(cluster==2) |> slice_max(avg_log2FC,n = 5))|> unique()
title_feature_umap_clus_2 <- "Macrophage IFN signature"
caption_feature_umap_clus_2 <- ""# 

clus_3_genes <- c(c("Cd3e", "Cd3g", "Cd28"), de_genes_tbl |> filter(cluster==3) |> slice_max(avg_log2FC,n = 5))|> unique()
title_feature_umap_clus_3 <- "T cells"
caption_feature_umap_clus_3 <- ""


clus_4_genes <- c(c("Ccl9", "Vcan", "Plac8", "Ccl6"), de_genes_tbl |> filter(cluster==4) |> slice_max(avg_log2FC,n = 5))|> unique()
title_feature_umap_clus_4 <- "Macrophage, Ccl9"
caption_feature_umap_clus_4 <- ""


clus_5_genes <- c(c("Csf1r", "Itgal", "Nr4a1"), de_genes_tbl |> filter(cluster==5) |> slice_max(avg_log2FC,n = 5))|> unique()
title_feature_umap_clus_5 <- "Macrophage, Csf1r"
caption_feature_umap_clus_5 <- ""


clus_6_genes <- c(c("Cd19", "Ms4a1", "Ighd"), de_genes_tbl |> filter(cluster==6) |> slice_max(avg_log2FC,n = 5))|> unique()
title_feature_umap_clus_6 <- "B cells"
caption_feature_umap_clus_6 <- ""

clus_7_genes <- c( c("H2-Aa", "Clec9a", "Ciita"), de_genes_tbl |> filter(cluster==7) |> slice_max(avg_log2FC,n = 5))|> unique()
title_feature_umap_clus_7 <- "Macrophage, antigen presentation"
caption_feature_umap_clus_7 <- ""

clus_8_genes <- c(c("Cd3d", "Mki67", "Top2a"), de_genes_tbl |> filter(cluster==8) |> slice_max(avg_log2FC,n = 5))|> unique()
title_feature_umap_clus_8 <- "T cell, proliferationg, cytotoxic"
caption_feature_umap_clus_8 <- ""

clus_9_genes <- c(c("Col4a1", "Cldn5", "Ly6c1", "Itga1"), de_genes_tbl |> filter(cluster==9) |> slice_max(avg_log2FC,n = 5))|> unique()
title_feature_umap_clus_9 <- "Epithelial cell, basement layer, OR Macrophage phagocytosed "
caption_feature_umap_clus_9 <-  ""

clus_10_genes <- c(c("Klre1", "Gzma", "Eomes"), de_genes_tbl |> filter(cluster==10) |> slice_max(avg_log2FC,n = 5))|> unique()
title_feature_umap_clus_10 <- "NK cell"
caption_feature_umap_clus_10 <-   ""


clus_11_genes <- c(c("Mrc1", "Cd9", "Ctsk", "Chil3", "Mertk"), de_genes_tbl |> filter(cluster==11) |> slice_max(avg_log2FC,n = 5))|> unique()
title_feature_umap_clus_11 <- "Macrophage"
caption_feature_umap_clus_11 <- ""


clus_12_genes <- c(c("S100a9", "Mmp9", "S100a9", "Il1b", "Ccr1"), de_genes_tbl |> filter(cluster==12) |> slice_max(avg_log2FC,n = 5))|> unique()
title_feature_umap_clus_12 <-  "Macrophage/monocyte inflammatory"
caption_feature_umap_clus_12 <- ""

clus_13_genes <- c("Ccr7", de_genes_tbl |> filter(cluster==13) |> slice_max(avg_log2FC,n = 5))|> unique()
title_feature_umap_clus_13 <- "Macrophage/dendritic"
caption_feature_umap_clus_13 <- "Ccr7"


clus_14_genes <- c(c("Ms4a2", "Il6", "Csf1", "Il13", "Cd63"), de_genes_tbl |> filter(cluster==14) |> slice_max(avg_log2FC,n = 5))|> unique()
title_feature_umap_clus_14 <- "Basophil,  IgE, receptor, (mast cell)"
caption_feature_umap_clus_14 <- ""


# clus_15_genes <- de_genes_tbl |> filter(cluster==15) |> slice_max(avg_log2FC,n = 5)
# title_feature_umap_clus_15 <- ""
# caption_feature_umap_clus_15 <- ""
# 
# clus_16_genes <- de_genes_tbl |> filter(cluster==16) |> slice_max(avg_log2FC,n = 5)
# title_feature_umap_clus_16 <- ""
# caption_feature_umap_clus_16 <- ""
# 
# clus_17_genes <- ""
# title_feature_umap_clus_17 <- ""
# caption_feature_umap_clus_17 <- ""
# 
# clus_18_genes <- ""
# title_feature_umap_clus_18 <- ""
# caption_feature_umap_clus_18 <- ""
# 
# clus_19_genes <- ""
# title_feature_umap_clus_19 <- ""
# caption_feature_umap_clus_19 <- ""
# 
# clus_20_genes <- ""
# title_feature_umap_clus_20 <- "??"
# caption_feature_umap_clus_20 <- ""
# 
# clus_21_genes <- ""
# title_feature_umap_clus_21 <- ""
# caption_feature_umap_clus_21 <- ""


###---evaluation helpers for chunks
ids_of_clusters <- unique(seurat_obj$seurat_clusters)
eval_0 <- 0 %in% ids_of_clusters
eval_1 <- 1 %in% ids_of_clusters
eval_2 <- 2 %in% ids_of_clusters
eval_3 <- 3 %in% ids_of_clusters
eval_4 <- 4 %in% ids_of_clusters
eval_5 <- 5 %in% ids_of_clusters
eval_6 <- 6 %in% ids_of_clusters
eval_7 <- 7 %in% ids_of_clusters
eval_8 <- 8 %in% ids_of_clusters
eval_9 <- 9 %in% ids_of_clusters
eval_10 <- 10 %in% ids_of_clusters
eval_11 <- 11 %in% ids_of_clusters
eval_12 <- 12 %in% ids_of_clusters
eval_13 <- 13 %in% ids_of_clusters
eval_14 <- 14 %in% ids_of_clusters
eval_15 <- 15 %in% ids_of_clusters
eval_16 <- 16 %in% ids_of_clusters
eval_17 <- 17 %in% ids_of_clusters
eval_18 <- 18 %in% ids_of_clusters
eval_19 <- 19 %in% ids_of_clusters
eval_20 <- 20 %in% ids_of_clusters
eval_21 <- 21 %in% ids_of_clusters
eval_22 <- 22 %in% ids_of_clusters
eval_23 <- 23 %in% ids_of_clusters#


eval_10_extra <- FALSE
eval_16_extra <-  FALSE
###---

theme_1 <- theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none", axis.title = element_blank())
