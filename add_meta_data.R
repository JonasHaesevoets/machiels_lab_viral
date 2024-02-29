
easypackages::libraries("Seurat", "tidyverse", "tidyseurat")
obj.v5 <- read_rds("C:\\Users\\danne\\R_projects\\machiels_lab_viral\\intermediate_data\\seurat_obj_central.rds")
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

obj.v5$sex_classifier_2 <- obj.v5 |> 
        join_features(features=c(y_chromosome_genes, "Xist"), slot="data", shape="wide") |>
        #as_tibble() |> 
        rowwise() |> 
        mutate(max_male_genes=max(c_across(Ddx3y:Kdm5d))
        ) |> 
        mutate(sex_classifier_2=case_when(
                (Xist==0)&(max_male_genes==0)~"no_expression_of_specific_genes",
                (Xist>0)&(max_male_genes==0)~"female",
                (Xist==0)&(max_male_genes>0)~"male",
                (Xist>0)&(max_male_genes>0)~"detection_of_male_and_female_genes",
                TRUE~"something is wrong")) |> 
                pull(sex_classifier_2)

obj.v5$sex_classifier_1 <- obj.v5 |>
        join_features(features=c(y_chromosome_genes, "Xist"), slot="data", shape="wide") |>
       # as_tibble() |> 
        rowwise() |> 
        mutate(max_male_genes=max(c_across(Ddx3y:Kdm5d))
        ) |> 
        #mutate(avg_male_gens=avg_male_gns) |>
        mutate(female_sex_ratio=Xist-max_male_genes) |>
        mutate(sex_classifier=case_when(
                female_sex_ratio==0~"no_expression_of_specific_genes",
                female_sex_ratio>0~"classified_as_female",
                female_sex_ratio<0~"classified_as_male",
                TRUE~"something is wrong")) |> 
        pull(sex_classifier)


##case_study_adjustment_of_Ms4a3_sorting_by_gabbr2
new_meta_data <- obj.v5 |>join_features("Gabbr2", slot="counts") |> 
        select(orig.ident, condition, sampletag_Ms4a3,.abundance_RNA) |>
        separate(orig.ident, c("exp", "sample_type"),sep="__") |>
        mutate(day=ifelse(exp=="viral.experiment.1", "d60", "d8")) |>
        mutate(day_sample_type=paste(day, sample_type, sep="_"),
               day_mock=ifelse(condition=="Mock","Mock", day),
               day_mock_sample_type=ifelse(condition=="Mock","Mock", day_sample_type),
               day_sample_type_cond=paste(day, sample_type, condition, sep="_"),
               day_sample_type_cond_ms4a3=paste(day, sample_type, condition,sampletag_Ms4a3, sep="_")) |> 
        #filter(.abundance_RNA==0) |> filter(sampletag_Ms4a3=="Ms4a3_neg")
        ####
        mutate(sample_tag_ms4a3_pos_gabbr2= case_when(
                .abundance_RNA>0 & sampletag_Ms4a3=="Ms4a3_neg" ~ "Gabbr2_pos_Ms4a3_neg",
                TRUE~sampletag_Ms4a3 )) |> 
        ####
        mutate(day_sample_type_cond_ms4a3_pos_gabbr2=
                       paste(day, sample_type, condition,sample_tag_ms4a3_pos_gabbr2, sep="_")) |> 
        separate(day_sample_type, c("day", "sample_type"), remove = FALSE)
        


obj.v5$day_sample_type <- pull(new_meta_data, day_sample_type)
obj.v5$day_mock <- pull(new_meta_data, day_mock)
obj.v5$day_mock_sample_type <- pull(new_meta_data, day_mock_sample_type)
obj.v5$day <- pull(new_meta_data, day)
obj.v5$day_sample_type_cond <- pull(new_meta_data, day_sample_type_cond)
obj.v5$day_sample_type_cond_ms4a3 <- pull(new_meta_data,day_sample_type_cond_ms4a3)
obj.v5$sample_type <- pull(new_meta_data,sample_type)


#add ms4a3-based calculation corrected with  Gabbr2 
obj.v5$day_sample_type_cond_ms4a3_pos_gabbr2 <- pull(new_meta_data, day_sample_type_cond_ms4a3_pos_gabbr2)
obj.v5$sample_tag_ms4a3_pos_gabbr2 <- pull(new_meta_data,sample_tag_ms4a3_pos_gabbr2)