library(Seurat)
library(tidyverse)
library(tidyseurat)
library(Matrix)



# read in all sequencing counts of all Rhapsody wells
# this means also those wells that the rhapsody pipeline does not call as containing a cell

## read in the empty
unfiltered_csv_paths <- c("..\\..\\Desktop/Analysis_Lucia/raw_data\\machiels_lab\\viral\\output_lung_d8\\BD-Analysis-BMachiels-Lung_DBEC_MolsPerCell_Unfiltered.csv.gz",
                          "..\\..\\Desktop/Analysis_Lucia/raw_data\\machiels_lab\\viral\\output_bal_d8\\BD-Analysis-BMachiels-BAL_DBEC_MolsPerCell_Unfiltered.csv.gz",
                          "..\\..\\Desktop/Analysis_Lucia/raw_data\\machiels_lab\\viral\\2023-10-02_output_lung\\Output_Lung\\BD-Analysis-BMachiels_DBEC_MolsPerCell_Unfiltered.csv.gz",
                          "../../Desktop/Analysis_Lucia/raw_data\\machiels_lab\\viral\\2023-10-02_output_bal\\Output_BAL\\BD-Analysis-BMachiels_DBEC_MolsPerCell_Unfiltered.csv.gz")

exp_name <- c("lung_d8","bal_d8",
              "lung_d60", "bal_d60")

for (i in seq_along(unfiltered_csv_paths)) {
  # we read only those proteins who were actually stained
  counts_fread = fread(unfiltered_csv_paths[i], select = c("Cell_Index",
                                                           "Siglec-F|Siglecf|AMM2013|pAbO",
                                                           "I-A_I-E|H2-Ab_Ad_Aq_Ed_Ek|AMM2019|pAbO",
                                                           "CD274|Cd274|AMM2038|pAbO",
                                                           "CD11c:HL3|Itgax|AMM2008|pAbO",
                                                           "Ly-6G|Ly6g|AMM2009|pAbO",
                                                           "Ly-6A_Ly-6E|Ly6a_Ly6e|AMM2026|pAbO"), 
                       showProgress = T)
  #rename the proteins
  colnames(counts_fread) <- c("cell_index","Siglecf_AbSeq","H2_ia_ie_AbSeq","Cd274_AbSeq","Cd11c","Ly_6g_AbSeq","Ly_6a_AbSeq")
  write_csv(counts_fread,paste0("../../Desktop/Analysis_Lucia/R_Projects/machiels_lab_viral-main/intermediate_data/",exp_name[i],"unfiltered_prot_counts_fread.csv") )
}
  
###################


# coming from the specific scripts for each of the data specific scripts which were ran in the QC dashboard
setup_chunks <- c(
  "dataset_QC_setup_chunk_exp1_bal.R",
  "dataset_QC_setup_chunk_exp2_bal.R",
  "dataset_QC_setup_chunk_exp1_lung.R",
  "dataset_QC_setup_chunk_exp2_lung.R")

# the counts also containing empty wells for the proteins
unfiltered_prot_counts <- c(
  "bal_d60unfiltered_prot_counts_fread.csv",
  "bal_d8unfiltered_prot_counts_fread.csv",
  "lung_d60unfiltered_prot_counts_fread.csv",
  "lung_d8unfiltered_prot_counts_fread.csv")

dataset_name <- c(
  "bal_d60",
  "bal_d8",
  "lung_d60",
  "lung_d8")

unfiltered_prot_counts <- unfiltered_prot_counts |> map_vec(\(x) paste0(".\\intermediate_data\\", x))

plot_list <- list()
for (i in seq_along(unfiltered_prot_counts)) {
  print(unfiltered_prot_counts[i])
  all_BD_protein_counts <- read_csv(unfiltered_prot_counts[i])
  
  protein_counts_cell_is_col <- all_BD_protein_counts |>  tidyfst::t_dt() 
  colnames(protein_counts_cell_is_col) <- protein_counts_cell_is_col[1,] 
  protein_counts_cell_is_col <- protein_counts_cell_is_col[-1,]
  protein_sum_tbl <- tibble(cell_sum=protein_counts_cell_is_col |>
                              colSums(), cell=colnames(protein_counts_cell_is_col))
  
  cutoff <- 6000
  protein_sum_tbl <- protein_sum_tbl 
  
  source(file=setup_chunks[i]) 
  print("in memory: seurat")
  seurat_obj <- seurat_obj |>
    mutate( kept_cell= case_when(
      (seurat_clusters %in% deleted_clusters) ~ "deleted_clusters",
      nCount_RNA>upper_nCount_RNA_thresh~ "too_much_RNA",
      nCount_RNA<lower_nCount_RNA_thresh~"too_little_RNA",
      percent_mito>upper_mito_thresh~"too_much_mito",
      sampletag_multiplets!="single_hashtag"~ "no_singlet",
      TRUE ~"keep"))
  
  meta_data_protein_sum_tbl <- seurat_obj |>
    as_tibble() |>
    rename(cell=".cell") |>
    full_join(protein_sum_tbl ) |>
    mutate(kept_cell=if_else(
      is.na(kept_cell),"empty_by_BD",kept_cell)
    )
  
  
  seurat_obj <- NULL
  
  
  ######
  #perform dsb normalization
  

  
  cells <- meta_data_protein_sum_tbl |> filter(kept_cell=="keep")|> pull(cell)
  cells_tbl <- all_BD_protein_counts |> filter(cell_index %in% cells) 
  cells_matrix <- t(as.matrix(cells_tbl|> select(-cell_index)))
  colnames(cells_matrix) <- cells_tbl |> pull(cell_index)
  
  empty_droplets <- meta_data_protein_sum_tbl |> filter(kept_cell=="empty_by_BD")|> pull(cell)
  empty_droplets_tbl <- all_BD_protein_counts |> filter(cell_index %in% empty_droplets) 
  empty_droplets_matrix <- t(as.matrix(empty_droplets_tbl|> select(-cell_index)))
  colnames(empty_droplets_matrix) <- empty_droplets_tbl |> pull(cell_index)
  
  cells.dsb.norm = DSBNormalizeProtein(
    cell_protein_matrix = cells_matrix, #
    empty_drop_matrix = empty_droplets_matrix, 
    denoise.counts = F, 
    use.isotype.control = F
  )
  
  write_rds(cells.dsb.norm,paste0(".\\intermediate_data\\dsb_matrix_",dataset_name[i],".rds"))
  
}


# write_rds(plot_list,paste0(".\\intermediate_data\\protein_reads_QC_plot_list.rds"))







##################

#path_1 <- "C:/Users/danne/R_projects/machiels_lab_viral/intermediate_data/seurat_obj_experiment_1_2_merged_v5_split_integrated.rds"
path_1 <-"..\\..\\R_projects\\machiels_lab_viral-main\\intermediate_data\\seurat_obj_experiment_1_2_integrated.rds"
obj.v5 <- read_rds(path_1)




exp_1_lung <- "intermediate_data/dsb_matrix_lung_d60.rds" |> read_rds() |> t() |> as_tibble(rownames="cell") |> mutate(cell=paste0("exp_1_lung_",cell))
exp_2_lung <- "intermediate_data/dsb_matrix_lung_d8.rds" |> read_rds() |> t() |> as_tibble(rownames="cell") |> mutate(cell=paste0("exp_2_lung_",cell))
exp_1_bal <-  "intermediate_data/dsb_matrix_bal_d60.rds" |> read_rds() |> t() |> as_tibble(rownames="cell") |> mutate(cell=paste0("exp_1_bal_",cell))
exp_2_bal <-"intermediate_data/dsb_matrix_bal_d8.rds" |> read_rds() |> t() |> as_tibble(rownames="cell") |> mutate(cell=paste0("exp_2_bal_",cell))

dsb_all <- bind_rows(exp_1_lung,
                     exp_2_lung,
                     exp_1_bal,
                     exp_2_bal)

dsb_features <- colnames(dsb_all)[-1]

# reduce protein count matrix to cells that are present in seurat object
dsb_all <- dsb_all |> filter(cell %in% (obj.v5 |> colnames()))
#append protein object bei "empty" cells with missing cell name
missing_cellnames <- setdiff(colnames(obj.v5), pull(dsb_all, cell))
dsb_all <- dsb_all |> filter(!is.na(cell)) |> bind_rows(tibble(cell=missing_cellnames))

dsb_all <- dsb_all[match(obj.v5 |> colnames(), dsb_all$cell),] 


dsb_all_t <- dsb_all |> select(-cell)|> t() |> as.sparse()
colnames(dsb_all_t) <- pull(dsb_all, cell)

dsb_all_seur <- CreateAssay5Object(counts=dsb_all_t )

proteins <- c("Siglecf-AbSeq"="siglec-f-siglecf-amm2013-p-ab-o",
              "H2-ia-ie-AbSeq"="i-a-i-e-h2-ab-ad-aq-ed-ek-amm2019-p-ab-o",
              "Cd274-AbSeq"="cd274-cd274-amm2038-p-ab-o",
              "Cd11c"="cd11c-hl3-itgax-amm2008-p-ab-o",
              "Ly-6g-AbSeq"="ly-6g-ly6g-amm2009-p-ab-o",
              "Ly-6a-AbSeq"="ly-6a-ly-6e-ly6a-ly6e-amm2026-p-ab-o")


adt_counts <- obj.v5@assays$protein$counts[names(proteins),]
rownames(adt_counts) <- names(proteins)


dsb_all_seur <- CreateAssay5Object(counts=adt_counts,  data=dsb_all_t )
#obj.v5[["RNA"]] <- JoinLayers(obj.v5[["RNA"]])#????????
options(Seurat.object.assay.version = "v5")

# 
obj.v5[["adt"]] <- dsb_all_seur


#change protein assay to v5 assay
DefaultAssay(obj.v5) <- "adt"
obj.v5 <- JoinLayers(obj.v5)


dsb <- obj.v5@assays$adt$data |> t() |>
        as_tibble(rownames=".cell") |>
        pivot_longer(cols = "Siglecf-AbSeq":"Ly-6a-AbSeq", names_to = "marker") |>
        mutate(method="prot_counts")


dsb <- dsb |> left_join(
        obj.v5@meta.data |>
                as_tibble(rownames = ".cell") |>
                select(.cell, orig.ident))

dsb <- dsb |>
        mutate(dsb_zero=abs(value))

dsb <-dsb |>
        group_by( marker, orig.ident) |>
        mutate(method_marker_0.9_quantile_dsb_zero = quantile(dsb_zero, probs = 0.9, na.rm = T)) |>
        mutate(method_marker_max_dsb_zero = quantile(dsb_zero, probs = 1, na.rm = T)) |>
        mutate(quantile_0.9_scaled_dsb_zero = dsb_zero / method_marker_0.9_quantile_dsb_zero) |>
        mutate(quantile_max_scaled_dsb_zero = dsb_zero / method_marker_max_dsb_zero) |> ungroup()

scaled_data_tbl <- dsb |>  select(.cell,marker,quantile_max_scaled_dsb_zero) |>
        pivot_wider(names_from = marker, values_from = quantile_max_scaled_dsb_zero)


scaled_data_mtx <- scaled_data_tbl |> select(-.cell) |> as.matrix()
rownames(scaled_data_mtx) <- pull(scaled_data_tbl, .cell)
#arrange according to 
scaled_data_mtx <- scaled_data_mtx[match(obj.v5 |> colnames(), scaled_data_tbl$.cell),] 
scaled_data_mtx_t <- scaled_data_mtx |> t()

obj.v5@assays$adt$scale.data <- scaled_data_mtx_t

obj.v5[["protein"]] <- NULL




write_rds(x = obj.v5,
          file ="../../Desktop/Analysis_Lucia/R_projects\\machiels_lab_viral-main\\intermediate_data\\seurat_obj_central.rds" )




