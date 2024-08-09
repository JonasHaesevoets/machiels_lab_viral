easypackages::libraries("tidyverse", "RColorBrewer", "fgsea", "GSA", "readxl","writexl", "clusterProfiler", "enrichplot", "ggplot2", "Seurat")

# load in objects
seurat_obj_experiment_1_bal_workflowed <- read_rds("../../Documents/machiels_lab_viral/intermediate_data/seurat_obj_experiment_1_combined_bal_raw_dbec.rds")
seurat_obj_experiment_2_bal_workflowed <- read_rds("../../Documents/machiels_lab_viral/intermediate_data/seurat_obj_experiment_2_combined_bal_raw_dbec.rds")
seurat_obj_experiment_1_lung_workflowed <- read_rds("../../Documents/machiels_lab_viral/intermediate_data/seurat_obj_experiment_1_combined_lung_raw_dbec.rds")
seurat_obj_experiment_2_lung_workflowed <- read_rds("../../Documents/machiels_lab_viral/intermediate_data/seurat_obj_experiment_2_combined_lung_raw_dbec.rds")
seurat_obj_experiment_1_bal_afterQCdashboard = read_rds("../../Documents/machiels_lab_viral/intermediate_data/seurat_obj_experiment_1_bal_afterQCdashboard.rds")
seurat_obj_experiment_2_bal_afterQCdashboard = read_rds("../../Documents/machiels_lab_viral/intermediate_data/seurat_obj_experiment_2_bal_afterQCdashboard.rds")
seurat_obj_experiment_2_lung_afterQCdashboard = read_rds("../../Documents/machiels_lab_viral/intermediate_data/seurat_obj_experiment_2_lung_afterQCdashboard.rds")
seurat_obj_experiment_1_lung_afterQCdashboard = read_rds("../../Documents/machiels_lab_viral/intermediate_data/seurat_obj_experiment_1_lung_afterQCdashboard.rds")
seurat_obj_full = read_rds("../../Documents/machiels_lab_viral/intermediate_data/seurat_obj_central_reclustered.rds")
seurat_obj_am = read_rds("../../Documents/machiels_lab_viral/intermediate_data/seurat_obj_central_am.rds")
seurat_d8_am = subset(seurat_obj_am, subset = day == "d8")
seurat_d60_am = subset(seurat_obj_am, subset = day == "d60")
seurat_d8_bal = subset(seurat_d8_am, subset = sample_type == "bal")
seurat_d8_lung = subset(seurat_d8_am, subset = sample_type == "lung")
seurat_d60_bal = subset(seurat_d60_am, subset = sample_type == "bal")
seurat_d60_lung = subset(seurat_d60_am, subset = sample_type == "lung")
seurat_d8_full = subset(seurat_obj_full, subset = day == "d8")
seurat_d60_full = subset(seurat_obj_full, subset = day == "d60")
seurat_d8_bal_full = subset(seurat_d8_full, subset = sample_type == "bal")
seurat_d8_lung_full = subset(seurat_d8_full, subset = sample_type == "lung")
seurat_d60_bal_full = subset(seurat_d60_full, subset = sample_type == "bal")
seurat_d60_lung_full = subset(seurat_d60_full, subset = sample_type == "lung")
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
  genes = c(genes_d60_bal, genes_d60_lung, genes_d8_bal, genes_d8_lung_qc),
  amount_cells = c(cells_d60_bal, cells_d60_lung, cells_d8_bal, cells_d8_lung),
  amount_cells_after_qc = c(cells_d60_bal_qc, cells_d60_lung_qc, cells_d8_bal_qc, cells_d8_lung_qc),
  percentage_filtered = c(percentage_d60_bal, percentage_d60_lung, percentage_d8_bal, percentage_d8_lung)
)




# amount cells in d60 bal
d60_bal = as.data.frame(table(seurat_obj_experiment_1_bal_workflowed$sampletag_name))
d60_bal_qc = as.data.frame(table(seurat_obj_experiment_1_bal_afterQCdashboard$sampletag_name))
d60_bal_am = as.data.frame(table(seurat_d60_bal$sampletag_name))
d60_bal_full = as.data.frame(table(seurat_d60_bal_full$sampletag_name))
write_xlsx(d60_bal, "d60_bal.xlsx")
write_xlsx(d60_bal_qc, "d60_bal_qc.xlsx")
write_xlsx(d60_bal_am, "d60_bal_am.xlsx")
write_xlsx(d60_bal_full, "d60_bal_full.xlsx")
# amount cells in d60 lung
d60_lung = as.data.frame(table(seurat_obj_experiment_1_lung_workflowed$sampletag_name))
d60_lung_qc = as.data.frame(table(seurat_obj_experiment_1_lung_afterQCdashboard$sampletag_name))
d60_lung_am = as.data.frame(table(seurat_d60_lung$sampletag_name))
d60_lung_full = as.data.frame(table(seurat_d60_lung_full$sampletag_name))
write_xlsx(d60_lung, "d60_lung.xlsx")
write_xlsx(d60_lung_qc, "d60_lung_qc.xlsx")
write_xlsx(d60_lung_am, "d60_lung_am.xlsx")
write_xlsx(d60_lung_full, "d60_lung_full.xlsx")
# amount cells in d8 bal
d8_bal = as.data.frame(table(seurat_obj_experiment_2_bal_workflowed$sampletag_name))
d8_bal_qc = as.data.frame(table(seurat_obj_experiment_2_bal_afterQCdashboard$sampletag_name))
d8_bal_am = as.data.frame(table(seurat_d8_bal$sampletag_name))
d8_bal_full = as.data.frame(table(seurat_d8_bal_full$sampletag_name))
write_xlsx(d8_bal, "d8_bal.xlsx")
write_xlsx(d8_bal_qc, "d8_bal_qc.xlsx")
write_xlsx(d8_bal_am, "d8_bal_am.xlsx")
write_xlsx(d8_bal_full, "d8_bal_full.xlsx")
# amount cells in d8 lung
d8_lung = as.data.frame(table(seurat_obj_experiment_2_lung_workflowed$sampletag_name))
d8_lung_qc = as.data.frame(table(seurat_obj_experiment_2_lung_afterQCdashboard$sampletag_name))
d8_lung_am = as.data.frame(table(seurat_d8_lung$sampletag_name))
d8_lung_full = as.data.frame(table(seurat_d8_lung_full$sampletag_name))
write_xlsx(d8_lung, "d8_lung.xlsx")
write_xlsx(d8_lung_qc, "d8_lung_qc.xlsx")
write_xlsx(d8_lung_am, "d8_lung_am.xlsx")
write_xlsx(d8_lung_full, "d8_lung_full.xlsx")










# amount cells in d60 bal
d60_bal_clusters = as.data.frame(table(seurat_obj_experiment_1_bal_workflowed$seurat_clusters))
d60_bal_clusters_qc = as.data.frame(table(seurat_obj_experiment_1_bal_afterQCdashboard$seurat_clusters))
d60_bal_clusters_am = as.data.frame(table(seurat_d60_bal$seurat_clusters))
d60_bal_clusters_full = as.data.frame(table(seurat_d60_bal_full$seurat_clusters))
write_xlsx(d60_bal_clusters, "d60_bal_clusters.xlsx")
# write_xlsx(d60_bal_clusters_qc, "d60_bal_clusters_qc.xlsx")
write_xlsx(d60_bal_clusters_am, "d60_bal_clusters_am.xlsx")
write_xlsx(d60_bal_clusters_full, "d60_bal_clusters_full.xlsx")
# amount cells in d60 lung
d60_lung_clusters = as.data.frame(table(seurat_obj_experiment_1_lung_workflowed$seurat_clusters))
 d60_lung_clusters_qc = as.data.frame(table(seurat_obj_experiment_1_lung_afterQCdashboard$seurat_clusters))
d60_lung_clusters_am = as.data.frame(table(seurat_d60_lung$seurat_clusters))
d60_lung_clusters_full = as.data.frame(table(seurat_d60_lung_full$seurat_clusters))
write_xlsx(d60_lung_clusters, "d60_lung_clusters.xlsx")
 write_xlsx(d60_lung_clusters_qc, "d60_lung_clusters_qc.xlsx")
write_xlsx(d60_lung_clusters_am, "d60_lung_clusters_am.xlsx")
write_xlsx(d60_lung_clusters_full, "d60_lung_clusters_full.xlsx")
# amount cells in d8 bal
d8_bal_clusters = as.data.frame(table(seurat_obj_experiment_2_bal_workflowed$seurat_clusters))
 d8_bal_clusters_qc = as.data.frame(table(seurat_obj_experiment_2_bal_afterQCdashboard$seurat_clusters))
d8_bal_clusters_am = as.data.frame(table(seurat_d8_bal$seurat_clusters))
d8_bal_clusters_full = as.data.frame(table(seurat_d8_bal_full$seurat_clusters))
write_xlsx(d8_bal_clusters, "d8_bal_clusters.xlsx")
 write_xlsx(d8_bal_clusters_qc, "d8_bal_clusters_qc.xlsx")
write_xlsx(d8_bal_clusters_am, "d8_bal_clusters_am.xlsx")
write_xlsx(d8_bal_clusters_full, "d8_bal_clusters_full.xlsx")
# amount cells in d8 lung
d8_lung_clusters = as.data.frame(table(seurat_obj_experiment_2_lung_workflowed$seurat_clusters))
 d8_lung_clusters_qc = as.data.frame(table(seurat_obj_experiment_2_lung_afterQCdashboard$seurat_clusters))
d8_lung_clusters_am = as.data.frame(table(seurat_d8_lung$seurat_clusters))
d8_lung_clusters_full = as.data.frame(table(seurat_d8_lung_full$seurat_clusters))
write_xlsx(d8_lung_clusters, "d8_lung_clusters.xlsx")
 write_xlsx(d8_lung_clusters_qc, "d8_lung_clusters_qc.xlsx")
write_xlsx(d8_lung_clusters_am, "d8_lung_clusters_am.xlsx")
write_xlsx(d8_lung_clusters_full, "d8_lung_clusters_full.xlsx")
