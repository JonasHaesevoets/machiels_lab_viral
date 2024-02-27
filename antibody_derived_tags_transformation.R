library(Seurat)
library(tidyverse)
library(tidyseurat)
library(Matrix)

#path_1 <- "C:/Users/danne/R_projects/machiels_lab_viral/intermediate_data/seurat_obj_experiment_1_2_merged_v5_split_integrated.rds"
path_1 <-"C:\\Users\\danne\\R_projects\\machiels_lab_viral\\intermediate_data\\seurat_obj_experiment_1_2_integrated.rds"
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
          file ="C:\\Users\\danne\\R_projects\\machiels_lab_viral\\intermediate_data\\seurat_obj_central.rds" )




