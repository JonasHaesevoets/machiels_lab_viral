####set up file paths
file_path <- vector("list") 
file_path$output <- ".\\output\\" 
file_path$intermediate_data<- ".\\intermediate_data\\" 
file_path$raw_data_root <- "C:\\Users\\danne\\raw_data\\machiels_lab\\viral\\"

#for read: seurat object path
#for write processed seurt object path
#sample_tag_calls path
#write path markers table


file_names_tibble <- tribble(
    ~raw_seurat_obj_path,
    ~path_sample_tag_calls,
    ~write_processed_seurat_obj_path,
    ~write_markers_path,#most parts of file name added by specificing parameters
    
    
    paste0(file_path$intermediate_data,"seurat_obj_experiment_1_combined_bal_raw_dbec.rds"),
    paste0(file_path$raw_data_root,"2023-10-02_output_bal\\Output_BAL\\","BD-Analysis-BMachiels_Sample_Tag_Calls.csv"),
    paste0(file_path$intermediate_data,"seurat_obj_experiment_1_combined_bal_raw_dbec_workflowed.rds"),
    paste0(file_path$intermediate_data, "experiment_1.")
    
    
)



seurat_basic_workflow <- function(seurat_obj, sample_tag_calls){
    #future::plan(multisession, workers = 4)
    #future::plan()
    #seurat_obj <- seurat_obj |> separate(.cell, into=c("discard","origin"), remove = F, sep = "__")
    seurat_obj <- seurat_obj |> PercentageFeatureSet("^Rp[sl]", col.name = "percent_ribo") # Check if all genes are included
    seurat_obj[["percent_mito"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Mt")
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(seurat_obj)
    seurat_obj <- ScaleData(seurat_obj, features = all.genes)
    seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
    
    seurat_obj$sample_tag_name <- sample_tag_calls  |>
        right_join(tibble(Cell_Index=colnames(seurat_obj))) |>
        pull(Sample_Name)
    
    seurat_obj <- seurat_obj |>  mutate(sampletag_multiplets=case_when(
        sample_tag_name=="Multiplet" ~ "multiplet",
        sample_tag_name=="Undetermined" ~"undeterminded",
        TRUE ~ "single_hashtag")
    )
    
    seurat_obj <- seurat_obj |>  mutate(sampletag_Ms4a3=case_when(
        sample_tag_name=="Multiplet" ~ "multiplet",
        sample_tag_name=="Undetermined" ~"undeterminded",
        str_detect(sample_tag_name, pattern="\\+") ~ "Ms4a3_pos",
        str_detect(sample_tag_name, pattern="\\-") ~ "Ms4a3_neg",
        TRUE ~ "single_hashtag")
    )
    
    seurat_obj$condition <- seurat_obj |> pull(sample_tag_name) |> str_split_i(i=1, pattern = "_")
    
    return(seurat_obj)
}

line=1 # add  for loop here later

path_raw_specific_dataset <- "2023-10-02_output_bal\\Output_BAL\\"

seurat_obj <- read_rds(file = file_names_tibble[[line,"raw_seurat_obj_path"]])

sample_tag_calls <- read_csv(file_names_tibble[[line,"path_sample_tag_calls"]], skip = 7) |>
    mutate(Cell_Index=as.character(Cell_Index))#paste0("bal_",Cell_Index)





#######
seurat_obj <- seurat_basic_workflow(seurat_obj, sample_tag_calls)
write_rds(seurat_obj,
          file = file_names_tibble[[line,"write_processed_seurat_obj_path"]] )
gc()
#####



#find markers for every cluster compared to all remaining cells, report only the positive
# ones
min_pct = 0.4
logfc_threshold = 0.25
max_cells_ = 300


seurat_obj.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE,
                                     min.pct = min_pct,
                                     logfc.threshold = logfc_threshold,
                                     max.cells.per.ident = max_cells_ )|> as_tibble(rownames = "gene")

write_rds(seurat_obj.markers, file = paste0(file_names_tibble[[line,"write_markers_path"]],
                                              "_QCmarkers_min.pct_",min_pct,
                                              "_logfc.threshold_",logfc_threshold,
                                              "_max.cells.per.ident_",max_cells_,
                                              ".rds"))

write_rds(seurat_obj.markers, file = paste0(file_names_tibble[[line,"write_markers_path"]],
                                            "_QCmarkers_min.pct_",min_pct,
                                            "_logfc.threshold_",logfc_threshold,
                                            "_max.cells.per.ident_",max_cells_,
                                            ".csv"))


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
