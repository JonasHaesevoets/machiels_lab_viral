easypackages::libraries("Seurat", "tidyverse", "tidyseurat")
obj.v5 <- read_rds("../../Documents/machiels_lab_viral/intermediate_data/seurat_obj_central.rds")
DefaultAssay(obj.v5) <- "RNA"

# Filter metadata columns in obj.v5 starting with "RNA_snn"
obj.v5$RNA_SNN_res_0.18 = obj.v5$RNA_snn_res.0.18
obj.v5$RNA_snn_res.0.18 = NULL
filtered_metadata <- obj.v5@meta.data[, !grepl("^RNA_snn", colnames(obj.v5@meta.data))]
obj.v5@meta.data <- filtered_metadata
print(colnames(obj.v5@meta.data))


# Filter metadata columns in obj.v5 starting with "harmony" except "harmony_clusters_0.18"
filtered_metadata <- obj.v5@meta.data[, !grepl("^harmony", colnames(obj.v5@meta.data))]
obj.v5@meta.data <- filtered_metadata
print(colnames(obj.v5@meta.data))

obj.v5 |>  as_tibble() |>
  separate(orig.ident, sep = "__", into = c("day", "sample_type"), remove = FALSE) |> 
  mutate(day=str_replace_all(day,c("viral.experiment.1"="d60",
                                   "viral.experiment.2"="d8") )) |> 
  mutate(day=factor(day, levels=c("Mock", "d8","d60")),
         sample_type=as_factor(sample_type),
         condition=as_factor(condition),
         sampletag_Ms4a3=as_factor(sampletag_Ms4a3),
         seurat_clusters=as_factor(seurat_clusters)
  )





##case_study_adjustment_of_Ms4a3_sorting_by_gabbr2
new_meta_data <- obj.v5 |>join_features("Gabbr2", slot="counts") |> 
  dplyr::select(orig.ident, condition, sampletag_Ms4a3,.abundance_RNA) |>
  separate(orig.ident, c("exp", "sample_type"),sep="__") |>
  mutate(day=ifelse(exp=="viral.experiment.1", "d60", "d8")) |>
  mutate(day_sample_type=paste(day, sample_type, sep="_"),
         day_mock=ifelse(condition=="Mock","Mock", day),
         day_mock_sample_type=ifelse(condition=="Mock","Mock", day_sample_type),
         day_mock_sample_type_mockTissueIncluded = ifelse(condition == "Mock", paste("Mock", sample_type, sep="_"), paste(day, sample_type, sep="_")),
         day_sample_type_cond=paste(day, sample_type, condition, sep="_"),
         day_mock_sample_type_cond=paste(day_mock_sample_type, condition, sep="_"),
         day_condition = paste(day, condition),
         day_sample_type_cond_ms4a3=paste(day, sample_type, condition,sampletag_Ms4a3, sep="_")) |> 
  #filter(.abundance_RNA==0) |> filter(sampletag_Ms4a3=="Ms4a3_neg")
  ####
  mutate(sample_tag_ms4a3_pos_gabbr2= case_when(
    .abundance_RNA>0 & sampletag_Ms4a3=="Ms4a3_neg" ~ "Gabbr2_pos_Ms4a3_neg",
    TRUE~sampletag_Ms4a3 )) |> 
  ####
  mutate(day_sample_type_cond_ms4a3_pos_gabbr2=
           paste(day, sample_type, condition,sample_tag_ms4a3_pos_gabbr2, sep="_")) |> 
  filter(!str_detect(day_sample_type_cond_ms4a3_pos_gabbr2, "Gabbr2_pos_Ms4a3_neg")) |> 
  separate(day_sample_type, c("day", "sample_type"), remove = FALSE)



obj.v5$day_sample_type <- pull(new_meta_data, day_sample_type)
obj.v5$day_mock <- pull(new_meta_data, day_mock)
obj.v5$day_condition = pull(new_meta_data, day_condition)
obj.v5$day_mock_sample_type <- pull(new_meta_data, day_mock_sample_type)
obj.v5$day_mock_sample_type_mockTissueIncluded = pull(new_meta_data, day_mock_sample_type_mockTissueIncluded)
obj.v5$day <- pull(new_meta_data, day)
obj.v5$day_sample_type_cond <- pull(new_meta_data, day_sample_type_cond)
obj.v5$day_sample_type_cond_ms4a3 <- pull(new_meta_data,day_sample_type_cond_ms4a3)
obj.v5$sample_type <- pull(new_meta_data,sample_type)
obj.v5$day_factor <- obj.v5 |>  mutate(day_factor=factor(day_mock, levels=c("Mock", "d8","d60"))) |> pull(day_factor)
obj.v5$day_mock_sample_type_cond <-  new_meta_data |>  pull(day_mock_sample_type_cond)
#add ms4a3-based calculation corrected with  Gabbr2 
obj.v5$day_sample_type_cond_ms4a3_pos_gabbr2 <- pull(new_meta_data, day_sample_type_cond_ms4a3_pos_gabbr2)
obj.v5$sample_tag_ms4a3_pos_gabbr2 <- pull(new_meta_data,sample_tag_ms4a3_pos_gabbr2)
obj.v5@meta.data$sampletag_name = factor(obj.v5@meta.data$sampletag_name,levels = c("Multiplet", "Undetermined", "Mock_Ms4a3_neg", 
                                                                                                 "Mock_Ms4a3_pos",
                                                                                                 "MAV1_Ms4a3_neg","MAV1_Ms4a3_pos","MuHV4_Ms4a3_neg", "MuHV4_Ms4a3_pos", "PR8_Ms4a3_neg",
                                                                                                 "PR8_Ms4a3_pos", "PVM_Ms4a3_neg", 
                                                                                                 "PVM_Ms4a3_pos"))
obj.v5$condition <- obj.v5 |> pull(sampletag_name) |> 
  str_split_i(i=1, pattern = "_") |> as_factor()
obj.v5$condition = factor(obj.v5$condition, levels = c("Mock", "MAV1", "MuHV4", "PR8", "PVM"))
obj.v5@meta.data$sampletag_Ms4a3 = factor(obj.v5@meta.data$sampletag_Ms4a3, levels = c("Ms4a3_pos", "Ms4a3_neg"))
obj.v5 |>  write_rds("../../Documents/machiels_lab_viral/intermediate_data/seurat_obj_central.rds")
