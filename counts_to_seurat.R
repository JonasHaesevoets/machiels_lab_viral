#this script loads rhapsody pipeline multimodal count data (tracripts, Abseq, cell hashtags)
#into seprate Seurat files per library and also combines seurat files per experiment and stores them in
#a "intermediate_data" file in the data analysis project folder as .rds
set.seed(2023)
#####
#load packeges
library("easypackages")
libraries("Seurat", "tidyverse", "janitor", "forcats","tidyseurat", "Matrix", "vroom", "tidyfst")
####
#below file paths are to be specified

raw_data_root <- "C:\\Users\\danne\\raw_data\\machiels_lab\\viral\\"

#name input and output file names
#dont use "_" for exeriment names because Seurat can't deal with this
file_names_tbl <- tribble(
        ~name_for_merged_seurat_file,
        ~raw_data_root,
        ~cell_id_A,
        ~cell_id_B, 
        ~dbec_file_path_A,
        ~name_for_seurat_file_A,  
        ~dbec_file_path_B,
        ~name_for_seurat_file_B,
        ~sample_tag_reads_per_cell_A,
        ~sample_tag_reads_per_cell_B,
        
        #experiment 2
        "viral.experiment.2", #name_for_merged_seurat_file
        "C:\\Users\\danne\\raw_data\\machiels_lab\\viral\\",#absolute path raw_data_root
        "lung", #cell_id_A
        "bal",#cell_id_B
        "output_lung_d8\\Combined_BD-Analysis-BMachiels-Lung_DBEC_MolsPerCell.csv",#dbec_file_path_A this is where the relevant count data is stored, "dbec" corrected,contains both ab-seq and transcriptome
        "seurat_obj_experiment_2_combined_lung_raw_dbec", #~name_for_seurat_file_A; file path for the intermediate data output of the script
        "output_bal_d8\\Combined_BD-Analysis-BMachiels-Bal_DBEC_MolsPerCell.csv",#~dbec_file_path_B,
        "seurat_obj_experiment_2_combined_bal_raw_dbec",#name_for_seurat_file_B
        "output_lung_d8\\BD-Analysis-BMachiels-Lung_Sample_Tag_ReadsPerCell.csv",#~sample_tag_reads_per_cell_A,
        "output_bal_d8\\BD-Analysis-BMachiels-Bal_Sample_Tag_ReadsPerCell.csv",#sample_tag_reads_per_cell_B,
        
        #experiment 1
        "viral.experiment.1", #name_for_merged_seurat_file
        "C:\\Users\\danne\\raw_data\\machiels_lab\\viral\\",#absolute path raw_data_root
        "lung", #cell_id_A
        "bal",#cell_id_B
        "2023-10-02_output_lung\\Output_Lung\\Combined_BD-Analysis-BMachiels_DBEC_MolsPerCell.csv",#dbec_file_path_A this is where the relevant count data is stored, "dbec" corrected,contains both ab-seq and transcriptome
        "seurat_obj_experiment_1_combined_lung_raw_dbec", #~name_for_seurat_file_A; file path for the intermediate data output of the script
        "2023-10-02_output_bal\\Output_BAL\\Combined_BD-Analysis-BMachiels_DBEC_MolsPerCell.csv",#~dbec_file_path_B,
        "seurat_obj_experiment_1_combined_bal_raw_dbec",#name_for_seurat_file_B
        "2023-10-02_output_lung\\Output_Lung\\BD-Analysis-BMachiels_Sample_Tag_ReadsPerCell.csv",#~sample_tag_reads_per_cell_A,
        "2023-10-02_output_bal\\Output_BAL\\BD-Analysis-BMachiels_Sample_Tag_ReadsPerCell.csv"#sample_tag_reads_per_cell_B,

)

####set up file paths
file_path <- vector("list") 
file_path$output <- ".\\output\\" 
file_path$intermediate_data<- ".\\intermediate_data\\" 
#file_path$raw_data_root <- raw_data_root


#### functions

##reads in rhapsody mutltiassay  csv file into seurat object
read_rhapsody_multi_assay_tibble_based <- function(dbec_counts_path, project_name, sample_tag_reads){
        counts <- vroom(dbec_counts_path, skip = 7) 
        barcodes <- counts |> pull(Cell_Index)
        
        #colnames(counts) <- features
        #counts <- as(counts, "sparseMatrix")
        print("count table loaded")

        protein <- counts[grep(names(counts),pattern = "pAbO")] |> clean_names() |>  #get rid of special characters that would interfere with downstream functions in feature names
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
        
        print("add sample tags")
        sample_tag_reads <- read_csv(sample_tag_reads, skip = 7) |>
                mutate(Cell_Index=as.character(Cell_Index)) |> clean_names()
        
        colnames(sample_tag_reads) <- str_remove_all(colnames(sample_tag_reads), "_mm_st_ab_o")
        barcodes <- pull(sample_tag_reads,cell_index)
        sample_tag_reads <- sample_tag_reads |> select(-cell_index) |> tidyfst::t_dt()
        colnames(sample_tag_reads) <- barcodes
        seurat[['sampletags']] <- CreateAssayObject(counts = sample_tag_reads)
        
        return(seurat)
}


##main

#iterate through rows of "file_names_tbl" and use the specified file paths of each row to get input files and store the output data files 
for (line in 1:nrow(file_names_tbl)) {
        #define absolute paths
        seurat_path_A_abolute <- paste0( file_names_tbl[[line,"raw_data_root"]],file_names_tbl[[line,"dbec_file_path_A"]])#Lung 
        seurat_path_B <- file_names_tbl[[line,"dbec_file_path_B"]]#BAL
        proj_name_A <- paste0(file_names_tbl[[line,"name_for_merged_seurat_file"]],
                              "__",
                              file_names_tbl[[line,"cell_id_A"]])
        proj_name_B <- paste0(file_names_tbl[[line,"name_for_merged_seurat_file"]],
                              "__",
                              file_names_tbl[[line,"cell_id_B"]])
        
        sample_tag_reads_A <- paste0( file_names_tbl[[line,"raw_data_root"]],
                                      file_names_tbl[[line,"sample_tag_reads_per_cell_A"]])
        print("Seurat A")

        seurat_A <-  read_rhapsody_multi_assay_tibble_based(dbec_counts_path= seurat_path_A_abolute,#absolute path
                                                            project_name = proj_name_A,
                                                            sample_tag_reads= sample_tag_reads_A)
        print("Seurat B")
        seurat_B <-  read_rhapsody_multi_assay_tibble_based(dbec_counts_path= paste0( file_names_tbl[[line,"raw_data_root"]],seurat_path_B
                                                                                      ),
                                                            project_name = proj_name_B,
                                                            sample_tag_reads= paste0( file_names_tbl[[line,"raw_data_root"]],
                                                                                      file_names_tbl[[line,"sample_tag_reads_per_cell_B"]]
                                                                                      )
                                                            )
       

       seurat_A |> write_rds(file = paste0(file_path$intermediate_data, file_names_tbl[[line,"name_for_seurat_file_A"]], ".rds"))
       seurat_B |> write_rds(file = paste0(file_path$intermediate_data, file_names_tbl[[line,"name_for_seurat_file_B"]], ".rds"))
       gc()
       cell_ids <- c(file_names_tbl[[line,"cell_id_A"]],file_names_tbl[[line,"cell_id_B"]])

       experiment<- merge(x = seurat_A,y=seurat_B,  add.cell.ids = cell_ids, project = file_names_tbl[[line,"name_for_merged_seurat_file"]])
       rm(seurat_A,seurat_B)
       gc()
       write_rds(x = experiment, file = paste0(file_path$intermediate_data, "seurat_obj_", file_names_tbl[[line,"name_for_merged_seurat_file"]], ".rds"))
}

