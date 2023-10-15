####set up file paths
raw_data_root <- "C:\\Users\\danne\\raw_data\\machiels_lab\\viral\\"
file_path <- vector("list") 
file_path$output <- ".\\output\\" 
file_path$intermediate_data<- ".\\intermediate_data\\" 
file_path$raw_data_root <- raw_data_root

seurat_obj <- read_rds(file = paste0(file_path$intermediate_data,"seurat_obj_experiment_1_combined_bal_raw_dbec.rds"))


seurat_basic_workflow <- function(seurat_obj){
        #future::plan(multisession, workers = 4)
        #future::plan()
        #seurat_obj <- seurat_obj |> separate(.cell, into=c("discard","origin"), remove = F, sep = "__")
        seurat_obj <- seurat_obj |> PercentageFeatureSet("^Rp[sl]", col.name = "percent_ribo") # Check if all genes are included
        seurat_obj[["percent_mito"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Mt")
        seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
        seurat_obj <- ScaleData(object = seurat_obj)
        seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
        all.genes <- rownames(seurat_obj)
        seurat_obj <- ScaleData(seurat_obj, features = all.genes)
        seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
        seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
        seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
        seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
        return(seurat_obj)
}






seurat_obj <- seurat_basic_workflow(seurat_obj)
write_rds(seurat_obj,file = paste0(file_path$intermediate_data,"seurat_obj_experiment_1_combined_bal_raw_dbec_workflowed.rds"))
gc()

#find markers for every cluster compared to all remaining cells, report only the positive
# ones
min_pct = 0.4
logfc_threshold = 0.25
max_cells_ = 300


seurat_obj.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.4, logfc.threshold = 0.25,max.cells.per.ident = 300 )|> as_tibble(rownames = "gene")

write_rds(experiment_1.markers, file = paste0(file_path$intermediate_data, "experiment_1.",
                                              "_QCmarkers_min.pct_",min_pct,
                                              "_logfc.threshold_",logfc_threshold,
                                              "_max.cells.per.ident_",max_cells_,
                                              ".rds"))


# experiment_1.markers_top <- experiment_1.markers %>%
#     group_by(cluster) %>%
#     slice_max(n = 2, order_by = avg_log2FC)
# 
# write_rds(experiment_1.markers, file = paste0(file_path$intermediate_data, "experiment_1._QCmarkers_min.pct_0.4_logfc.threshold_0.25_max.cells.per.ident_300.rds"))
# 
# 
# 
# experiment_1.markers_top <- read_rds(file = paste0(file_path$intermediate_data, "experiment_1._QCmarkers_min.pct_0.4_logfc.threshold_0.25_max.cells.per.ident_300.rds")) %>%
#         group_by(cluster) %>%
#         slice_max(n = 2, order_by = avg_log2FC)
# 
# experiment_1.markers_top_4 <- read_rds(file = paste0(file_path$intermediate_data, "experiment_1._QCmarkers_min.pct_0.4_logfc.threshold_0.25_max.cells.per.ident_300.rds")) %>%
#         group_by(cluster) %>%
#         slice_max(n = 4, order_by = avg_log2FC)
# 
# 
# levels_diff_exp_features_4 <-   pull(experiment_1.markers_top_4, gene) |> unique()
# experiment_1.markers_top_4 <- experiment_1.markers_top_4 |> mutate(gene=factor(gene, levels=levels_diff_exp_features_4))

#write_rds(experiment_1.markers, file = paste0(file_path$intermediate_data, "experiment_1._QCmarkers_min.pct_0.4_logfc.threshold_0.25_max.cells.per.ident_300.rds"))
