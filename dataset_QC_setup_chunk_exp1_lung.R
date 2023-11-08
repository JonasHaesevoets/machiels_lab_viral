
set.seed(2023)
#set up file paths
file_path <- vector("list") 
file_path$output <- ".\\output\\" 
file_path$intermediate_data<- ".\\intermediate_data\\" 
#file_path$raw_data <- "C:\\Users\\danne\\raw_data\\machiels_lab\\viral\\2023-10-02_output_lung\\Output_Lung\\BD-Analysis-BMachiels_Expression_Data_Unfiltered.st.gz"

file_name_obj <- "seurat_obj_experiment_1_combined_lung_raw_dbec_workflowed.rds"

seurat_obj <- read_rds(file = paste0(file_path$intermediate_data,file_name_obj)) |> mutate(sampletag_multiplets=factor(sampletag_multiplets, levels=c("undeterminded" , "multiplet"  ,   "single_hashtag"))) #level change can be ommited later when intermediate data has been asjusted

#seurat_obj <-   read_rds(file = paste0(file_path$intermediate_data,"SAMPLED_2000cells_seurat_obj_experiment_1_combined_lung_raw_dbec_workflowed.rds") )

de_genes_tbl <- "C:\\Users\\danne\\R_projects\\machiels_lab_viral\\intermediate_data\\experiment_1._lung__QCmarkers_min.pct_0.4_logfc.threshold_0.25_max.cells.per.ident_300.rds" |> read_rds() #|> select(cluster:gene)


de_genes_tbl <- read_csv("intermediate_data/experiment_1.Lung_QCmarkers_min.pct_0.4_logfc.threshold_0.25_max.cells.per.ident_300.csv")|>
        relocate(gene) 

file_path_output <- paste0(file_path$intermediate_data,"seurat_obj_experiment_1_combined_lung_raw_dbec_cleaned.rds")


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
deleted_clusters=c(0,5,7,10,15,16)

del_clusters_string <- paste0(as.character(deleted_clusters))




#genes to gighlight in feature umaps
clus_0_genes <- c("Prf1","Il2rb","Gzmb","Gzma","Eomes","Ms4a4b","Il12rb2","H2-Q7", "Cd3g", "Cd8a", "Cd4", "Ifng","Il2ra", "Foxp3","Ncr1", "Rorc", "Ncam1", "Tbx21", "Il7r", "Klrc1", "Klrk1", "Ncr3")
title_feature_umap_clus_0 <- ""
caption_feature_umap_clus_0 <- "likely MHC2 expressing NK cells" # 

clus_1_genes <- c("Ccr2","Ifitm3","Thbs1", "Ms4a4c")
title_feature_umap_clus_1 <- ""
caption_feature_umap_clus_1 <- "likely MHC2 expressing NK cells" # 

clus_2_genes <- c("Ace", "Csf1r")
title_feature_umap_clus_2 <- ""
caption_feature_umap_clus_2 <- "" # 

clus_3_genes <- c("Fpr", "Chil3","Lpl", "Cd9", "Ctsd", "Mertk")
title_feature_umap_clus_3 <- ""
caption_feature_umap_clus_3 <- "scavanger Macrophages" 


clus_4_genes <- c("Pld3", "Ftl1")
title_feature_umap_clus_4 <- ""
caption_feature_umap_clus_4 <- ""


clus_5_genes <- c("Cd79a", "Cd19", "Ms4a1")
title_feature_umap_clus_5 <- ""
caption_feature_umap_clus_5 <- "B lineage"


clus_6_genes <- c("Cd209a", "Flt3", "H2-Ab1", "H2-Aa")
title_feature_umap_clus_6 <- ""
caption_feature_umap_clus_6 <- "MHC-2 high"

clus_7_genes <- c("Cd3e", "Cd8a")
title_feature_umap_clus_7 <- ""
caption_feature_umap_clus_7 <- "CD8 T cells"

clus_8_genes <- c("C1qb", "CD81", "Apoe", "H2-Aa")
title_feature_umap_clus_8 <- ""
caption_feature_umap_clus_8 <- "Apoe alveolar macrophages"

clus_9_genes <- c("H2-Eb1", "Mgl2", "Tnfsf9", "H2-Aa")
title_feature_umap_clus_9 <- ""
caption_feature_umap_clus_9 <-  "CD137L?, Apoe alveolar macrophages"

clus_10_genes <- c("mt-Nd3", "mt-Co3" ,"Lpp")
title_feature_umap_clus_10 <- ""
caption_feature_umap_clus_10 <-"mt high"


clus_11_genes <- c("Clec9a", "Cst3", "Cd24a", "H2-Ab1")
title_feature_umap_clus_11 <- ""
caption_feature_umap_clus_11 <- "Clec9a, MHC2"


clus_12_genes <- c("Mki67", "Cks1b", "Ctsk" , "Tpx2", "Spp1", "H2-Aa")
title_feature_umap_clus_12 <- "Proliferating/Spp1+" 
caption_feature_umap_clus_12 <- "Ki-67, Tpx: Spindle assembly factor required for normal assembly of mitotic spindles"

clus_13_genes <- c("Ccr7", "Il4i1", "Ccl22" , "Socs2", "Relb", "Il12b")
title_feature_umap_clus_13 <- "Ccr7/Il4i1" 
caption_feature_umap_clus_13 <- ""


clus_14_genes <- c("Mki67", "Top2a", "H2-Ab1")
title_feature_umap_clus_14 <- "Proliferating"
caption_feature_umap_clus_14 <- ""


clus_15_genes <- c("Stmn1", "Mcm5", "Top2a")
title_feature_umap_clus_15 <- "Proliferating"
caption_feature_umap_clus_15 <- ""

clus_16_genes <- c( "Cldn5", "Col4a1","Ly6c1", "H2-Aa")
title_feature_umap_clus_16 <- "Macrophage - Epithelial cell doublet/phagocytosed" 
caption_feature_umap_clus_16 <- "Caudin5, collagen 4, tight junction, Col4 producing, MHC2 positive, Epidermal Growth Factor Receptor 5 expressing"

clus_17_genes <- c( "Csf1", "Ms4a2", "Il6","Cd63")
title_feature_umap_clus_17 <- "Csf1 producing, IgE-R+"
caption_feature_umap_clus_17 <- "Ms4a2= IgE-R, Fcer1a, "

clus_18_genes <- c( "S100a9","Mmp9", "S100a8", "H2-Q10")
title_feature_umap_clus_18 <- "Procalcitonin/Mmp9"
caption_feature_umap_clus_18 <- ""

clus_19_genes <- c( "Cidec", "Kcnn3","F7", "Mertk", "Cd24a")
title_feature_umap_clus_19 <- "??"
caption_feature_umap_clus_19 <- "Cidec: Lipid transferase, "

clus_20_genes <- c( "F5", "Acod1", "Cd24a", "Ccr3")
title_feature_umap_clus_20 <- "??"
caption_feature_umap_clus_20 <- ""

clus_21_genes <- c( "F13a1", "Mafb", "Csf3r", "Ccr2", "Gzmb", "S100a4")
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


eval_10_extra <- TRUE
eval_16_extra <-  TRUE
###---

theme_1 <- theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none", axis.title = element_blank())
