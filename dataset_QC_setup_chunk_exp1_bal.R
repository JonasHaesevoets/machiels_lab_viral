
set.seed(2023)
#set up file paths
file_path <- vector("list") 
file_path$output <- ".\\output\\" 
file_path$intermediate_data<- ".\\intermediate_data\\" 
#file_path$raw_data <- "C:\\Users\\danne\\raw_data\\machiels_lab\\viral\\2023-10-02_output_bal\\Output_bal\\BD-Analysis-BMachiels_Expression_Data_Unfiltered.st.gz"

obj_identifier <- "experiment_1_bal"
file_name_obj <- paste0("seurat_obj_", obj_identifier, "_workflowed.rds")
#file_name_obj <- "seurat_obj_experiment_1_bal_workflowed.rds"

seurat_obj <- read_rds(file = paste0(file_path$intermediate_data,file_name_obj)) |> mutate(sampletag_multiplets=factor(sampletag_multiplets, 
                                                                                                                       levels=c("undeterminded" , "multiplet"  ,   "single_hashtag"))) #level change can be ommited later when intermediate data has been asjusted
#seurat_obj <- subset(x = seurat_obj, downsample = 100)



#seurat_obj <-   read_rds(file = paste0(file_path$intermediate_data,"SAMPLED_2000cells_seurat_obj_experiment_1_combined_bal_raw_dbec_workflowed.rds") )

#de_genes_tbl <- "C:\\Users\\danne\\R_projects\\machiels_lab_viral\\intermediate_data\\experiment_1._bal__QCmarkers_min.pct_0.4_logfc.threshold_0.25_max.cells.per.ident_300.rds" |> read_rds() #|> select(cluster:gene)


de_genes_tbl <- read_csv("intermediate_data/experiment_1_workflowed.Bal_QCmarkers_min.pct_0.4_logfc.threshold_0.25_max.cells.per.ident_300.csv")|>
        relocate(gene) 


file_path_output <- paste0(file_path$intermediate_data,"seurat_obj_experiment_1_bal_afterQCdashboard.rds")


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
deleted_clusters=c(5,10,13)

cluster_determination <- c("low_RNA_5"=5,
                           "low_RNA_mt_high_10"=10,
                           "contamination_B_NK"=13)

seurat_obj <- seurat_obj |> mutate(orig_cluster_specification=case_when(
        seurat_clusters==5 ~ "5_low_RNA_5",
        seurat_clusters==10 ~ "10_low_RNA_mt_high",
        seurat_clusters==13 ~ "13_contamination_B_NK",
        TRUE~seurat_clusters
        
))

del_clusters_string <- paste0(as.character(deleted_clusters))




#genes to gighlight in feature umaps
clus_0_genes <-  c(de_genes_tbl |> filter(cluster==0) |> slice_max(avg_log2FC,n = 5) |> pull(gene), "Cd2")
title_feature_umap_clus_0 <- ""
caption_feature_umap_clus_0 <- "" # 

clus_1_genes <- de_genes_tbl |> filter(cluster==1) |> slice_max(avg_log2FC,n = 5)|> pull(gene)
title_feature_umap_clus_1 <- ""
caption_feature_umap_clus_1 <- "" # 

clus_2_genes <- de_genes_tbl |> filter(cluster==2) |> slice_max(avg_log2FC,n = 5)|> pull(gene)
title_feature_umap_clus_2 <- ""
caption_feature_umap_clus_2 <- "" # 

clus_3_genes <- de_genes_tbl |> filter(cluster==3) |> slice_max(avg_log2FC,n = 5)|> pull(gene)
title_feature_umap_clus_3 <- ""
caption_feature_umap_clus_3 <- "" 


clus_4_genes <- de_genes_tbl |> filter(cluster==4) |> slice_max(avg_log2FC,n = 5)|> pull(gene)
title_feature_umap_clus_4 <- ""
caption_feature_umap_clus_4 <- ""


clus_5_genes <- de_genes_tbl |> filter(cluster==5) |> slice_max(avg_log2FC,n = 5)|> pull(gene)
title_feature_umap_clus_5 <- ""
caption_feature_umap_clus_5 <- ""


clus_6_genes <- de_genes_tbl |> filter(cluster==6) |> slice_max(avg_log2FC,n = 5)|> pull(gene)
title_feature_umap_clus_6 <- ""
caption_feature_umap_clus_6 <- ""

clus_7_genes <- de_genes_tbl |> filter(cluster==7) |> slice_max(avg_log2FC,n = 5)|> pull(gene)
title_feature_umap_clus_7 <- ""
caption_feature_umap_clus_7 <- ""

clus_8_genes <- de_genes_tbl |> filter(cluster==8) |> slice_max(avg_log2FC,n = 5)|> pull(gene)
title_feature_umap_clus_8 <- ""
caption_feature_umap_clus_8 <- " "

clus_9_genes <- de_genes_tbl |> filter(cluster==9) |> slice_max(avg_log2FC,n = 5)|> pull(gene)
title_feature_umap_clus_9 <- ""
caption_feature_umap_clus_9 <-  ""

clus_10_genes <- de_genes_tbl |> filter(cluster==10) |> slice_max(avg_log2FC,n = 5)|> pull(gene)
title_feature_umap_clus_10 <- ""
caption_feature_umap_clus_10 <-""


clus_11_genes <- de_genes_tbl |> filter(cluster==11) |> slice_max(avg_log2FC,n = 5)|> pull(gene)
title_feature_umap_clus_11 <- ""
caption_feature_umap_clus_11 <- ""


clus_12_genes <- de_genes_tbl |> filter(cluster==12) |> slice_max(avg_log2FC,n = 5)|> pull(gene)
title_feature_umap_clus_12 <-  
caption_feature_umap_clus_12 <- 

clus_13_genes <- de_genes_tbl |> filter(cluster==13) |> slice_max(avg_log2FC,n = 5)|> pull(gene)
title_feature_umap_clus_13 <- ""
caption_feature_umap_clus_13 <- ""


clus_14_genes <- de_genes_tbl |> filter(cluster==14) |> slice_max(avg_log2FC,n = 5)|> pull(gene)
title_feature_umap_clus_14 <- ""
caption_feature_umap_clus_14 <- ""


clus_15_genes <- de_genes_tbl |> filter(cluster==15) |> slice_max(avg_log2FC,n = 5)|> pull(gene)
title_feature_umap_clus_15 <- ""
caption_feature_umap_clus_15 <- ""

clus_16_genes <- de_genes_tbl |> filter(cluster==16) |> slice_max(avg_log2FC,n = 5)|> pull(gene)
title_feature_umap_clus_16 <- ""
caption_feature_umap_clus_16 <- ""

clus_17_genes <- ""
title_feature_umap_clus_17 <- ""
caption_feature_umap_clus_17 <- ""

clus_18_genes <- ""
title_feature_umap_clus_18 <- ""
caption_feature_umap_clus_18 <- ""

clus_19_genes <- ""
title_feature_umap_clus_19 <- ""
caption_feature_umap_clus_19 <- ""

clus_20_genes <- ""
title_feature_umap_clus_20 <- "??"
caption_feature_umap_clus_20 <- ""

clus_21_genes <- ""
title_feature_umap_clus_21 <- ""
caption_feature_umap_clus_21 <- ""


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
