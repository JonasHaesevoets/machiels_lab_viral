library(easypackages)
libraries("tidyverse", "Seurat", "tidyseurat", "tibble")

####set up file paths

file_path <- vector("list") 
file_path$output <- "../../Documents/machiels_lab_viral/output"
file_path$intermediate_data<- "../../Documents/machiels_lab_viral/intermediate_data/"
file_path$raw_data_root <- "../../Documents/machiels_lab_viral/raw_data/machiels_lab/viral/"

#for read: seurat object path
#for write processed seurt object path
#sample_tag_calls path
#write path markers table


file_names_tibble <- tribble(
    ~library_name,
    ~raw_seurat_obj_path,
    ~path_sample_tag_calls,
    ~write_processed_seurat_obj_path,
    ~write_markers_path,#most parts of file name added by specificing parameters
    
    
    "Lung_d8",#library_name
    paste0(file_path$intermediate_data,"seurat_obj_experiment_2_combined_lung_raw_dbec.rds"),#raw_seurat_obj_path
    paste0(file_path$raw_data_root,"output_lung_d8\\","BD-Analysis-BMachiels-Lung_Sample_Tag_Calls.csv"),#path_sample_tag_calls
    paste0(file_path$intermediate_data,"seurat_obj_experiment_2_lung_workflowed.rds"),#write_processed_seurat_obj_path
    paste0(file_path$intermediate_data, "experiment_2_workflowed.Lung"),#write_markers_path

    
    "BAL_d8",#library_name
    paste0(file_path$intermediate_data,"seurat_obj_experiment_2_combined_bal_raw_dbec.rds"),#raw_seurat_obj_path
    paste0(file_path$raw_data_root,"output_bal_d8\\","BD-Analysis-BMachiels-Bal_Sample_Tag_Calls.csv"),#path_sample_tag_calls
    paste0(file_path$intermediate_data,"seurat_obj_experiment_2_bal_workflowed.rds"),#write_processed_seurat_obj_path
    paste0(file_path$intermediate_data, "experiment_2_workflowed.BAL"),#write_markers_path
    
    
    "Lung",#library_name
    paste0(file_path$intermediate_data,"seurat_obj_experiment_1_combined_lung_raw_dbec.rds"),#raw_seurat_obj_path
    paste0(file_path$raw_data_root,"2023-10-02_output_lung\\Output_Lung\\","BD-Analysis-BMachiels_Sample_Tag_Calls.csv"),#path_sample_tag_calls
    paste0(file_path$intermediate_data,"seurat_obj_experiment_1_lung_workflowed.rds"),#write_processed_seurat_obj_path
    paste0(file_path$intermediate_data, "experiment_1_workflowed.Lung"),#write_markers_path


    "BAL",#library_name
    paste0(file_path$intermediate_data,"seurat_obj_experiment_1_combined_bal_raw_dbec.rds"),#raw_seurat_obj_path
    paste0(file_path$raw_data_root,"2023-10-02_output_bal\\Output_BAL\\","BD-Analysis-BMachiels_Sample_Tag_Calls.csv"),#path_sample_tag_calls
    paste0(file_path$intermediate_data,"seurat_obj_experiment_1_bal_workflowed.rds"),#write_processed_seurat_obj_path
    paste0(file_path$intermediate_data, "experiment_1_workflowed.BAL")#write_markers_path



    
)


path_raw_specific_dataset <- "2023-10-02_output_bal\\Output_BAL\\"


seurat_basic_workflow <- function(seurat_obj, sample_tag_calls){
  
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
    
    #renaming of the sampletag
    seurat_obj$sampletag_name <- sample_tag_calls  |>
        right_join(tibble(Cell_Index=colnames(seurat_obj))) |> 
        mutate(Sample_Name=str_replace(Sample_Name,pattern="\\-", replacement = "_")) |> 
       mutate(Sample_Name=str_replace_all(Sample_Name,pattern="\\+", replacement = "_pos")) |> 
       mutate(Sample_Name=str_replace_all(Sample_Name,pattern="\\-", replacement = "_neg")) |> 
        pull(Sample_Name) |>  as_factor()
    
    seurat_obj <- seurat_obj |>  mutate(sampletag_multiplets=case_when(
        sampletag_name=="Multiplet" ~ "multiplet",
        sampletag_name=="Undetermined" ~"undeterminded",
        TRUE ~ "single_hashtag")
    ) |> 
        mutate(sampletag_multiplets=factor(sampletag_multiplets, levels=c("undeterminded" , "multiplet"  ,   "single_hashtag"))) 
    #create separate groups based on Ms4a3
    seurat_obj <- seurat_obj |>  mutate(sampletag_Ms4a3=case_when(
        sampletag_name=="Multiplet" ~ "multiplet",
        sampletag_name=="Undetermined" ~"undeterminded",
        str_detect(sampletag_name, pattern="pos") ~ "Ms4a3_pos",
        str_detect(sampletag_name, pattern="neg") ~ "Ms4a3_neg",
        TRUE ~ "single_hashtag")
    )|> 
        mutate(sampletag_Ms4a3=as_factor(sampletag_Ms4a3))
    
    seurat_obj$condition <- seurat_obj |> pull(sampletag_name) |> str_split_i(i=1, pattern = "_") |> as_factor()
    seurat_obj$virus <- tibble(sampletag_name=seurat_obj |> pull(sampletag_name)) |> separate(sampletag_name,into = c("virus")) |> pull("virus")
    
    return(seurat_obj)
}


for (line in 1:nrow(file_names_tibble)) {

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

seurat_obj.markers <- FindAllMarkers(seurat_obj,
                                     only.pos = TRUE,
                                     min.pct = min_pct,
                                     logfc.threshold = logfc_threshold,
                                     max.cells.per.ident = max_cells_ ) |>
    as_tibble()

write_rds(seurat_obj.markers, file = paste0(file_names_tibble[[line,"write_markers_path"]],
                                              "_QCmarkers_min.pct_",min_pct,
                                              "_logfc.threshold_",logfc_threshold,
                                              "_max.cells.per.ident_",max_cells_,
                                              ".rds"))

write_csv(seurat_obj.markers, file = paste0(file_names_tibble[[line,"write_markers_path"]],
                                            "_QCmarkers_min.pct_",min_pct,
                                            "_logfc.threshold_",logfc_threshold,
                                            "_max.cells.per.ident_",max_cells_,
                                            ".csv"))


}


