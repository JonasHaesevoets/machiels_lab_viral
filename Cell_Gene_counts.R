# load in objects
seurat_obj_experiment_1_bal_workflowed <- read_rds("../../Documents/machiels_lab_viral/intermediate_data/seurat_obj_experiment_1_bal_workflowed.rds")
seurat_obj_experiment_2_bal_workflowed <- read_rds("../../Documents/machiels_lab_viral/intermediate_data/seurat_obj_experiment_2_bal_workflowed.rds")
seurat_obj_experiment_1_lung_workflowed <- read_rds("../../Documents/machiels_lab_viral/intermediate_data/seurat_obj_experiment_1_lung_workflowed.rds")
seurat_obj_experiment_2_lung_workflowed <- read_rds("../../Documents/machiels_lab_viral/intermediate_data/seurat_obj_experiment_2_lung_workflowed.rds")
seurat_obj_experiment_1_bal_afterQCdashboard = read_rds("seurat_obj_experiment_1_bal_afterQCdashboard.rds")
seurat_obj_experiment_2_bal_afterQCdashboard = read_rds("seurat_obj_experiment_2_bal_afterQCdashboard.rds")
seurat_obj_experiment_2_lung_afterQCdashboard = read_rds("seurat_obj_experiment_2_lung_afterQCdashboard.rds")
seurat_obj_experiment_1_lung_afterQCdashboard = read_rds("seurat_obj_experiment_1_lung_afterQCdashboard.rds")
# genes
genes_d60_bal = nrow(seurat_obj_experiment_1_bal_workflowed)
genes_d60_lung = nrow(seurat_obj_experiment_1_lung_workflowed)
genes_d8_bal = nrow(seurat_obj_experiment_2_bal_workflowed)
genes_d8_lung = nrow(seurat_obj_experiment_2_lung_workflowed)

# genes after qc
genes_d60_bal_qc = nrow(seurat_obj_experiment_1_bal_afterQCdashboard)
genes_d60_lung_qc = nrow(seurat_obj_experiment_1_lung_afterQCdashboard)
genes_d8_bal_qc = nrow(seurat_obj_experiment_2_bal_afterQCdashboard)
genes_d8_lung_qc = nrow(seurat_obj_experiment_2_lung_afterQCdashboard)



# cells
cells_d60_bal = ncol(seurat_obj_experiment_1_bal_workflowed)
cells_d60_lung = ncol(seurat_obj_experiment_1_lung_workflowed)
cells_d8_bal = ncol(seurat_obj_experiment_2_bal_workflowed)
cells_d8_lung = ncol(seurat_obj_experiment_2_lung_workflowed)
# cells after qc

cells_d60_bal_qc = ncol(seurat_obj_experiment_1_bal_afterQCdashboard)
cells_d60_lung_qc = ncol(seurat_obj_experiment_1_lung_afterQCdashboard)
cells_d8_bal_qc = ncol(seurat_obj_experiment_2_bal_afterQCdashboard)
cells_d8_lung_qc = ncol(seurat_obj_experiment_2_lung_afterQCdashboard)


## percentage filtered cells per case

percentage_d60_bal = cells_d60_bal_qc / cells_d60_bal
percentage_d60_lung = cells_d60_lung_qc / cells_d60_lung
percentage_d8_bal = cells_d8_bal_qc / cells_d8_bal
percentage_d8_lung = cells_d8_lung_qc / cells_d8_lung

amountCellsBeforeAfterQC <- data.frame(
  sample = c("bal_d60", "lung_d60", "bal_d8", "lung_d8"),
  genes = c(genes_d60_bal, genes_d60_lung, genes_d8_bal, genes_d8_lung_qc)
  amount_cells = c(18000, 14000, 11000, 11500),
  amount_cells_after_qc = c(13372, 7485, 5533, 6858),
  percentage_filtered = c(percentage_d60_bal, percentage_d60_lung, percentage_d8_bal, percentage_d8_lung)
)