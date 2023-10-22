#this script loads rhapsody pipeline mutlimodal(tracripts, Abseq, cell hashtags) count data
#into seprate Seurat files per library and also combines seurat files per experiment and stores them in
#a "intermediate_data" file in the data analysis project folder as .rds
set.seed(2023)
#####
#load packeges
library("easypackages")
libraries("Seurat", "tidyverse", "here", "janitor", "forcats","tidyseurat", "Matrix", "vroom", "tidyfst")
####
#below file paths are to be specified

raw_data_root <- "C:\\Users\\danne\\raw_data\\machiels_lab\\viral\\"

#name input and output file names
#dont use "_" for exeriment names because Seurat can't deal with this
file_names_tbl <- tribble(
~name_for_merged_seurat_file,~cell_id_A, ~cell_id_B,  ~dbec_file_path_A,     ~name_for_seurat_file_A,     ~dbec_file_path_B,     ~name_for_seurat_file_B,
"viral.experiment.1", "lung", "bal", "2023-10-02_output_lung\\Output_Lung\\Combined_BD-Analysis-BMachiels_DBEC_MolsPerCell.csv", "seurat_obj_experiment_1_combined_lung_raw_dbec", 
"2023-10-02_output_bal\\Output_BAL\\Combined_BD-Analysis-BMachiels_DBEC_MolsPerCell.csv", "seurat_obj_experiment_1_combined_bal_raw_dbec"

)

####set up file paths
file_path <- vector("list") 
file_path$output <- ".\\output\\" 
file_path$intermediate_data<- ".\\intermediate_data\\" 
file_path$raw_data_root <- raw_data_root


#### functions
read_rhapsody_multi_assay_tibble_based <- function(dbec_counts_path, project_name){##reads in rhapsody combined file into seurat mutltiassay
        counts <- vroom(dbec_counts_path, skip = 7) 
        barcodes <- counts |> pull(Cell_Index)
        
        #colnames(counts) <- features
        #counts <- as(counts, "sparseMatrix")
        print("count table loaded")

        protein <- counts[grep(names(counts),pattern = "pAbO")] |> clean_names() |>  #get rid of special chraacters that would interfere with downstream functions in feature names
                tidyfst::t_dt() #efficiently transposes dataframes 
        colnames(protein) <- barcodes # reapply cell index
        
        
        transcriptome <- counts[grep(names(counts),pattern = c("pAbO|Cell"), invert = T)] 
        print("transposing transcriptome tibble")
        transcriptome <- transcriptome |> tidyfst::t_dt()
        colnames(transcriptome) <- barcodes
        rm(counts)# freeing up memory
        gc() #freeing up memory
        
        print("creating seurat object")
        seurat <- CreateSeuratObject(counts = transcriptome, project=project_name,)
        
        print("adding assay objects")
        seurat[['protein']] <- CreateAssayObject(counts = protein)
        return(seurat)
}



#iterate through rows of "file_names_tbl" and use the specified names to get and store the data files
for (line in nrow(file_names_tbl)) {
        seurat_path_A <- file_names_tbl[[line,"dbec_file_path_A"]]
        seurat_path_B <- file_names_tbl[[line,"dbec_file_path_B"]]
        proj_name_A <- paste0(file_names_tbl[[line,"name_for_merged_seurat_file"]],
                              "__",
                              file_names_tbl[[line,"cell_id_A"]])
        proj_name_B <- paste0(file_names_tbl[[line,"name_for_merged_seurat_file"]],
                              "__",
                              file_names_tbl[[line,"cell_id_B"]])

       seurat_A <-  read_rhapsody_multi_assay_tibble_based(paste0(file_path$raw_data_root,seurat_path_A), project_name = proj_name_A)
       seurat_B <-  read_rhapsody_multi_assay_tibble_based(paste0(file_path$raw_data_root,seurat_path_B),project_name = proj_name_A )
       
       seurat_A |> write_rds(file = paste0(file_path$intermediate_data, file_names_tbl[[line,"name_for_seurat_file_A"]], ".rds"))
       seurat_B |> write_rds(file = paste0(file_path$intermediate_data, file_names_tbl[[line,"name_for_seurat_file_B"]], ".rds"))
       gc()
       cell_ids <- c(file_names_tbl[[line,"cell_id_A"]],file_names_tbl[[line,"cell_id_B"]])
       
       experiment_1 <- merge(seurat_A,seurat_B,  add.cell.ids = cell_ids, project = file_names_tbl[[line,"name_for_merged_seurat_file"]])
       rm(seurat_A,seurat_B)
       gc()
       write_rds(experiment_1, file = paste0(file_path$intermediate_data, "seurat_obj_", file_names_tbl[[line,"name_for_merged_seurat_file"]], ".rds"))
}

