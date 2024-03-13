easypackages::libraries("Seurat", "tidyverse", "tidyseurat", "scooter")
obj.v5 <- read_rds("intermediate_data/seurat_obj_central.rds")
DefaultAssay(obj.v5) <- "RNA"


obj.v5 |>  as_tibble() |>
        separate(orig.ident, sep = "__", into = c("day", "sample_type"), remove = FALSE) |> 
        mutate(day=str_replace_all(day,c("viral.experiment.1"="d60",
                                         "viral.experiment.2"="d8") )) |> 
        mutate(day=factor(day, levels=c("Mock", "d8","d60")),
               sample_type=as_factor(sample_type),
               condition=as_factor(condition),
               sampletag_Ms4a3=as_factor(sampletag_Ms4a3),
               harmony_cluster_8dims_rough=as_factor(harmony_cluster_8dims_rough)
        )

## from case_study_DeSeq2_exp1_sex_demultiplex.qmd

y_chromosome_genes<-   c("Ddx3y", "Eif2s3y", "Kdm5d")

obj.v5$sex_classifier <- obj.v5 |> 
        join_features(features=c(y_chromosome_genes, "Xist"), slot="data", shape="wide") |>
        #as_tibble() |> 
        rowwise() |> 
        mutate(max_male_genes=max(c_across(Ddx3y:Kdm5d))) |> 
        mutate(sex_classifier=case_when(
                (Xist==0)&(max_male_genes==0)~"no_expression_of_specific_genes",
                (Xist>0)&(max_male_genes==0)~"female",
                (Xist==0)&(max_male_genes>0)~"male",
                (Xist>0)&(max_male_genes>0)~"detection_of_male_and_female_genes",
                TRUE~"something is wrong"))


##case_study_adjustment_of_Ms4a3_sorting_by_gabbr2
new_meta_data <- obj.v5 |>join_features("Gabbr2", slot="counts") |> 
        select(orig.ident, condition, sampletag_Ms4a3,.abundance_RNA, sex_classifier) |>
        separate(orig.ident, c("exp", "sample_type"),sep="__") |>
        mutate(day=ifelse(exp=="viral.experiment.1", "d60", "d8")) |>
        filter(sex_classifier %in% c("male","female")) |> 
        rename(gender = sex_classifier) |> 
        mutate(day_sample_type=paste(day, sample_type, sep="_"),
               day_mock=ifelse(condition=="Mock","Mock", day),
               day_mock_sample_type=ifelse(condition=="Mock","Mock", day_sample_type),
               day_sample_type_cond=paste(day, sample_type, condition, sep="_"),
               day_mock_sample_type_cond=paste(day_mock_sample_type, condition, sep="_"),
               day_sample_type_cond_ms4a3=paste(day, sample_type, condition,sampletag_Ms4a3, sep="_")) |> 
        #filter(.abundance_RNA==0) |> filter(sampletag_Ms4a3=="Ms4a3_neg")
        ####
        mutate(sample_tag_ms4a3_pos_gabbr2= case_when(
                .abundance_RNA>0 & sampletag_Ms4a3=="Ms4a3_neg" ~ "Gabbr2_pos_Ms4a3_neg",
                TRUE~sampletag_Ms4a3 )) |> 
        ####
        mutate(day_sample_type_cond_ms4a3_pos_gabbr2=
                       paste(day, sample_type, condition,sample_tag_ms4a3_pos_gabbr2, sep="_")) |> 
        separate(day_sample_type, c("day", "sample_type"), remove = FALSE) |> 
        filter(!str_detect(day_sample_type_cond_ms4a3_pos_gabbr2, "Gabbr2_pos_Ms4a3_neg")) |> 
        mutate(day_sample_type_cond_ms4a3_pos_gabbr2_gender=
                       paste(day, sample_type, condition,sample_tag_ms4a3_pos_gabbr2, gender, sep="_"))
        
obj.v5$day_sample_type <- pull(new_meta_data, day_sample_type)
obj.v5$day_mock <- pull(new_meta_data, day_mock)
obj.v5$day_mock_sample_type <- pull(new_meta_data, day_mock_sample_type)
obj.v5$day <- pull(new_meta_data, day)
obj.v5$day_sample_type_cond <- pull(new_meta_data, day_sample_type_cond)
obj.v5$day_sample_type_cond_ms4a3 <- pull(new_meta_data,day_sample_type_cond_ms4a3)
obj.v5$sample_type <- pull(new_meta_data,sample_type)
obj.v5$day_factor <- obj.v5 |>  mutate(day_factor=factor(day, levels=c("Mock", "d8","d60"))) |> pull(day_factor)
obj.v5$day_mock_sample_type_cond <-  new_meta_data |>  pull(day_mock_sample_type_cond)
#add ms4a3-based calculation corrected with  Gabbr2 
obj.v5$day_sample_type_cond_ms4a3_pos_gabbr2 <- pull(new_meta_data, day_sample_type_cond_ms4a3_pos_gabbr2)
obj.v5$sample_tag_ms4a3_pos_gabbr2 <- pull(new_meta_data,sample_tag_ms4a3_pos_gabbr2)
obj.v5$day_tissue_virus_ms4a3_gabbr2corrected_gender <- pull(new_meta_data, day_sample_type_cond_ms4a3_pos_gabbr2_gender)
obj.v5$gender =  pull(new_meta_data, gender)

obj.v5 |>  write_rds("intermediate_data/seurat_obj_central.rds")
