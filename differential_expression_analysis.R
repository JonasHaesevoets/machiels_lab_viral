## Differential expression analysis using the central seurat object containing all the clean data ready for analysis


################################################### 1) Load in the required packages and the central seurat object
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("scran")
easypackages::libraries(
  "Seurat","tidyseurat","dplyr","ggplot2","readr","forcats","ragg","ggpubr","stringr","Matrix","textTinyR","tidyverse","ggrepel","scran","pheatmap","DESeq2","patchwork","ggridges", "scran", "textTinyR", "withr", "htmltools","pheatmap", "viridis", "ggExtra")

obj = read_rds("../Desktop/Analysis_Lucia/R_Projects/machiels_lab_viral-main/intermediate_data/seurat_obj_central.rds")






DefaultAssay(obj) <- "RNA" ## make sure that the right assay is used
#################################################### 2) control for sex related batch effects

all_conditions_set <- c("all", "only_viruses")
cond_select <- "all"# change this to change the comparisons of the analysis


############## this function extracts+analyzes expression data from specific cells 
make_idents_counts <- function(
    ds,
    idents_choose="seurat_clusters",
    assay_="RNA",
    min_cell_num=200,
    keep=NULL,
    perc_expressed=1){
  slot_ <- "counts"
  Idents(ds) <- idents_choose
  group_names <- unique(Idents(ds))
  if (!is.null(keep)) {group_names <- keep
  
  }
  reads <- LayerData(obj,"counts")
  #get rid of false genes / . in there name
  reads <- reads[rownames(reads) %>% str_subset(pattern = "[.]",negate = T),]
  #cut out genes with less then perc_expressed percent of cells that eypress the gene
  i_cells <- WhichCells(ds, idents = keep)
  reads_group <- reads[,i_cells]
  
  reads <- reads[(((reads_group > 1) %>% rowSums())> (length(i_cells)/100*perc_expressed)),]
  print(paste0(as.character(nrow(reads))," genes kept"))  # Extracts the expression matrix from the specified assay and filters genes based on perc_expressed
  sum_counts_ds <- tibble(.rows =  nrow(reads))
  #z=1
  for (i in seq_along(group_names)) { #Loops through the unique cell identities (groups) in ds
    print(i)
    #print(colnames(sum_counts_ds))
    print(group_names[i])
    i_cells <- WhichCells(object =ds,idents = group_names[i]) #Selects cells belonging to the current group using WhichCells.
    
    if (length(i_cells)<min_cell_num) { # the group has fewer cells than min_cell_num, add a column with zero counts and a descriptive name to the results frame and skip to the next group.
      sum_counts_ds <- sum_counts_ds %>%
        add_column(dings=0)
      colnames(sum_counts_ds)[i] <- paste0(as.character(group_names[i]),"_smaller" , as.character(min_cell_num)) #Adds a column to the results frame containing the counts for the current group and assigns the group name as the column name.
      
      next}
    #i_cells <- WhichCells(object =ds,idents = "hla_d_rlo_s100a_monocytes_C19-CB-0003")
    reads_group <- reads[,i_cells] %>%  sparse_Sums(rowSums=T) # selects part of the expression matrix (reads), keeping only those genes expressed in cells beloning to the current group
    
    sum_counts_ds <- sum_counts_ds %>%
      add_column(dings=reads_group) ## ads a column called dings whose values are filled with the calculated group-wise gene expression sums stored in reads_group
    colnames(sum_counts_ds)[i] <- as.character(group_names[i]) ## renames the newly added "dings" column to the actual group name corresponding to the current iteration (i) in the loo
    #print(as.character(colnames(sum_counts_ds)[i]))
    
  }
  
  sum_counts_ds <- sum_counts_ds %>%  mutate(gene=rownames(reads)) %>% relocate(gene) #: Adds a column with gene names to the results frame 
  return(sum_counts_ds)
}





################## add sex specific genes to the object which is done to check for batch effects due to inequal genders and create plots for downstream correction of this batch effects

y_chromosome_genes<-   c("Ddx3y", "Eif2s3y", "Kdm5d")

obj_sex <- obj |>
  join_features(features=c(y_chromosome_genes, "Xist"), slot="data", shape="wide") |>
  rowwise() |> 
  mutate(max_male_genes=max(c_across(Ddx3y:Kdm5d))
  ) |> 
  mutate(female_sex_ratio=Xist-max_male_genes) |>
  mutate(sex_classifier=case_when(
    (Xist==0)&(max_male_genes==0)~"no_expression_of_specific_genes",
    (Xist>0)&(max_male_genes==0)~"female",
    (Xist==0)&(max_male_genes>0)~"male",
    (Xist>0)&(max_male_genes>0)~"detection_of_male_and_female_genes",
    TRUE~"something is wrong")) 

### create a plot showing the distribution of gender associated genes among the ms4a3 positive and negative controls for the 4 different viral conditions, stored as plot_distribution_sex_genes
plot_distribution_sex_genes <-
  ggplot(data = obj_sex,aes(day_sample_type_cond_ms4a3, fill=sex_classifier))+
  geom_bar()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+ ggtitle("sex classification", "female: (Xist>0)AND(max_male_genes==0);male= Xist== AND(max_male_genes>0)")

plot_distribution_sex_genes


### create a histogram for the sex genes in question for all conditions (8, 4 viruses with ms4a3 + and -)
theme_x1 <- theme_minimal()+
  theme(strip.text.x = element_text(size = 15))


p1 <- obj_sex |> ggplot(aes(Ddx3y)) + geom_histogram(binwidth = 0.1)+facet_wrap( ~day_sample_type_cond_ms4a3, ncol = 1 )+theme_x1
p2 <- obj_sex |> ggplot(aes(Eif2s3y))+ geom_histogram(binwidth = 0.1)+facet_wrap( ~day_sample_type_cond_ms4a3, ncol = 1)+theme_x1
p3 <- obj_sex |> ggplot(aes(Kdm5d))+ geom_histogram(binwidth = 0.1)+facet_wrap( ~day_sample_type_cond_ms4a3, ncol = 1)+theme_x1
p4 <- obj_sex |> ggplot(aes(max_male_genes))+ geom_histogram(binwidth = 0.1)+facet_wrap( ~day_sample_type_cond_ms4a3, ncol = 1)+theme_x1

obj_sex |>
  select(.cell,day_sample_type_cond_ms4a3,sex_classifier ) |>
  mutate(paste0(day_sample_type_cond_ms4a3,sex_classifier, sep="_"))

histogram_sexGenes_allconditions = p1+p2+p3+p4
#histogram_sexGenes_allconditions


### create a plot showing the difference between the Xist gene expression and the max male associated genes
plot_difference_xist_maxmalegenes =
  ggplot(data = obj_sex,aes(max_male_genes+1, Xist+1)) + geom_point() +ylim(0,3)     +xlim(0,3) +facet_wrap( ~day_sample_type_cond_ms4a3) + ggtitle("Xist vs. max. male gene")
plot_difference_xist_maxmalegenes

### create a plots showing the densities of Eif2s3y/Kdm5d/maxmalegene/ across all samples
plot_Eif2s3y <- ggplot(data = obj_sex,aes(Eif2s3y, day_sample_type_cond_ms4a3, fill=day_sample_type_cond_ms4a3)) + geom_density_ridges()+ggtitle("Eif2s3y per sample")
plot_Eif2s3y
plot_Kdm5d <- ggplot(data = obj_sex,aes(Kdm5d, day_sample_type_cond_ms4a3, fill=day_sample_type_cond_ms4a3)) + geom_density_ridges()+ggtitle("Kdm5d per sample")
plot_Kdm5d
plot_maxmale <- ggplot(data = obj_sex,aes(max_male_genes, day_sample_type_cond_ms4a3, fill=day_sample_type_cond_ms4a3)) + geom_density_ridges()+ggtitle("maxmale per sample")
plot_maxmale


################## based on the subsequent plots it was decided to filter by this sex classifier and classify them either as m or f while removing those who are androgenous
keep <- obj_sex |> filter(sex_classifier %in% c("male","female")) |> select(.cell,sex_classifier)

obj_sex_corrected <- obj |> filter(.cell %in% pull(keep,.cell)) |> left_join(keep) # left join so only those from the left dataframe where a match in the right dataframe are retained
obj_sex_corrected <- obj_sex_corrected |> mutate(virus_ms4a3_sex=paste(day_sample_type_cond_ms4a3_pos_gabbr2,sex_classifier, sep="-"))

obj_sex_corrected <- obj_sex_corrected %>%
  filter(!str_detect(day_sample_type_cond_ms4a3_pos_gabbr2, "Gabbr2_pos_Ms4a3_neg")) ## remove those positive for Gabbr2 but negative for Ms4a3

obj_sex_corrected_sums <- make_idents_counts(ds=obj_sex_corrected,idents_choose = "virus_ms4a3_sex",perc_expressed = 0.05) ## call the functions specified on the ds with the specified idents and a min percentage that a gene should be expressed

## bar chart agg counts per sample
### from now on obj_sex_corrected will be used for further analysis



#################################################### 3) check for the effects of possible wrong sorting on ms4a3 => to be validated

# ##################### creata a new dataframe for metadata to be used later on based on the sex corrected seurat object
# # new_meta_data <- obj_sex_corrected_sums |>join_features("Gabbr2", slot="counts") |> ## adds a new column Gabbr2 containing expression counts for Gabbr2
# #   select(orig.ident, condition, sampletag_Ms4a3,.abundance_RNA) |> ## select only specific columns
# #   separate(orig.ident, c("exp", "sample_type"),sep="__") |> ## separate the orig.ident column into 2 based on the experiment type and the amount of days
# #   mutate(day=ifelse(exp=="viral.experiment.1", "d60", "d8")) |>
# #   mutate(day_sample_type=paste(day, sample_type, sep="_"),
# #          day_mock=ifelse(condition=="Mock","Mock", day),
# #          day_mock_sample_type=ifelse(condition=="Mock","Mock", day_sample_type),
# #          day_sample_type_cond=paste(day, sample_type, condition, sep="_"),
# #          day_sample_type_cond_ms4a3=paste(day, sample_type, condition,sampletag_Ms4a3, sep="_")) |> 
# #   #filter(.abundance_RNA==0) |> filter(sampletag_Ms4a3=="Ms4a3_neg")
# #   ####
#   mutate(sample_tag_ms4a3_pos_gabbr2= case_when(
#     .abundance_RNA>0 & sampletag_Ms4a3=="Ms4a3_neg" ~ "Gabbr2_pos_Ms4a3_neg",
#     TRUE~sampletag_Ms4a3 )) |> 
#   ####Assigns "Gabbr2_pos_Ms4a3_neg" if .abundance_RNA is greater than 0 and sampletag_Ms4a3 is "Ms4a3_neg" (potentially identifying cells with Gabbr2 expression but negative for Ms4a3 sorting)
#   #### Otherwise, keeps the original value from sampletag_Ms4a3.
#   mutate(day_sample_type_cond_ms4a3_pos_gabbr2=
#            paste(day, sample_type, condition,sample_tag_ms4a3_pos_gabbr2, sep="_"))
#   
# 
# #add ms4a3-based calculation corrected with  Gabbr2 
# obj_sex_corrected$day_sample_type_cond_ms4a3_pos_gabbr2 <- pull(new_meta_data, day_sample_type_cond_ms4a3_pos_gabbr2)
# obj_sex_corrected$sample_tag_ms4a3_pos_gabbr2 <- pull(new_meta_data,sample_tag_ms4a3_pos_gabbr2)
# 
# ##################### identifies differentially expressed markers distinguishing cells based on their Ms4a3 expression status and prepares for further analysis potentially involving different experimental conditions.
# 
# Idents(obj_sex_corrected) <- "sampletag_Ms4a3"
# DefaultAssay(obj_sex_corrected) <- "RNA"
# 
# #calculate top markers for Ms4a3 pos and neg- use these markers later to calculate auc for individual sammples
# markers_ms4a3_neg <- FindMarkers(obj_sex_corrected ,
#                                  ident.1 =  "Ms4a3_neg", ident.2 =  "Ms4a3_pos", test.use = "roc", max.cells.per.ident = 1000,min.pct = 0.2,logfc.threshold = 0.5 ,only.pos = T) |> as_tibble(rownames="gene")|> mutate("sampletag_Ms4a3"="Ms4a3_neg")
# 
# markers_ms4a3_pos <- FindMarkers(obj_sex_corrected , ident.1 =  "Ms4a3_pos", ident.2 =  "Ms4a3_neg", test.use = "roc", max.cells.per.ident = 1000,min.pct = 0.2,logfc.threshold = 0.5 ,only.pos = T)|> as_tibble(rownames="gene") |> mutate("sampletag_Ms4a3"="Ms4a3_pos")
# 
# markers_ms4a3 <- bind_rows(head(markers_ms4a3_neg),head(markers_ms4a3_pos)) ## combine the top markers from Ms4a3 + and - conditions into a single dataframe
# 
# conditions <- unique(obj_sex_corrected$condition)
# 
# 
# ################## Identifies markers specifically for Ms4a3_pos vs. Ms4a3_neg cells within each unique combination of day, sample type, and condition.
# 
# day_sample_type_cond_v <- unique(obj_sex_corrected$day_sample_type_cond) ## extracts unique values from the day_sample_type_condition
# Idents(obj_sex_corrected) <- "day_sample_type_cond_ms4a3" ## sets the idendity of cells in the seurat object to be based on the specified column talking about Ms4a3 status en day sample type and condition
# all_neg_pos_comps_tbl <- tibble() ## creates an empty tibble to be filled up using the subsequent for loop
# for (i in day_sample_type_cond_v) {
#   neg <- paste(i, "Ms4a3_neg", sep="_")
#   pos <- paste(i, "Ms4a3_pos" ,sep="_")
#   ## identify differentially expressed markers between the Ms4a3 + and Ms4a3- cells within the current condition
#   x1 <- FindMarkers(obj_sex_corrected, ident.1 = pos, ident.2 = neg, features = unique(pull(markers_ms4a3, gene)), test.use = "roc") |> as_tibble(rownames="gene") |> mutate(comparison=paste0(neg,"__vs.__",pos))
#   all_neg_pos_comps_tbl<- bind_rows(all_neg_pos_comps_tbl,x1)
#   
# }
# write_rds(all_neg_pos_comps_tbl,"../Desktop/Analysis_Lucia/R_Projects/machiels_lab_viral-main/intermediate_data//ms4a3_roc_specific_markers.rds")
# all_neg_pos_comps_tbl <- read_rds("../Desktop/Analysis_Lucia/R_Projects/machiels_lab_viral-main/intermediate_data/ms4a3_roc_specific_markers.rds")
# 
# 
# obj_sex_corrected$day_sample_type_cond_ms4a3_pos_gabbr2 |> table()
# 
# # calculate 
# Idents(obj_sex_corrected) <- "day_sample_type_cond_ms4a3_pos_gabbr2"
# 
# gabbr2_cor_neg_pos_comps_tbl <- tibble()
# for (i in day_sample_type_cond_v) {
#   neg <- paste(i, "Ms4a3_neg", sep="_")
#   pos <- paste(i, "Ms4a3_pos" ,sep="_")
#   x1 <- FindMarkers(obj_sex_corrected, ident.1 = pos, ident.2 = neg, features = unique(pull(markers_ms4a3, gene)), test.use = "roc") |> as_tibble(rownames="gene") |> mutate(comparison=paste0(neg,"__vs.__",pos))
#   gabbr2_cor_neg_pos_comps_tbl<- bind_rows(gabbr2_cor_neg_pos_comps_tbl,x1)
#   
# }
# 
# 
# 
# write_rds(gabbr2_cor_neg_pos_comps_tbl,"../Desktop/Analysis_Lucia/R_Projects/machiels_lab_viral-main/intermediate_data/ms4a3_roc_specific_markers_pos_gabbr2.rds")
# gabbr2_cor_neg_pos_comps_tbl <- read_rds("../Desktop/Analysis_Lucia/R_Projects/machiels_lab_viral-main/intermediate_data/ms4a3_roc_specific_markers_pos_gabbr2.rds")
# gabbr2_cor_neg_pos_comps_tbl
# all_neg_pos_comps_tbl
# 
# output_path <- "../Desktop/Analysis_Lucia/R_Projects/machiels_lab_viral-main/intermediate_data/"
# 
# gabbr2_cor_neg_pos_comps_tbl <- gabbr2_cor_neg_pos_comps_tbl |> mutate(correction="gabbr2_cut") ##  Adds a new column named "correction" with the value "gabbr2_cut", meaning it is corrected for GabaR2
# auc_top_markers_ms4a3_with_and_without_gabbr2_correction <-  all_neg_pos_comps_tbl |> mutate(correction="no_correction") |> bind_rows(gabbr2_cor_neg_pos_comps_tbl) |>
#   filter(gene %in% c("Gabbr2", "Siglecf", "Csf1r", "Actn1", "Cd2")) |>
#   ggplot(aes(myAUC, comparison,color=correction))+geom_point()+ xlim(0,1)+ facet_wrap(~gene, ncol = 1)
# 
# ggsave(auc_top_markers_ms4a3_with_and_without_gabbr2_correction,
#        filename="auc_top_markers_ms4a3_with_and_without_gabbr2_correction.svg",
#        path=output_path,
#        width=10 , height = 14)
# ## this generates a plot showing the top differentiating genes between the Ms4a3 condition
# 
# plot_difference_xist_maxmalegenes
# ## to note from this graph, there is a clear difference between the distance between them between day 8 and day 60
# 
# ################## Check for the functions of the top markers both negative and positive
# ## Siglec f in ms4a3 negative =>
# 


#################################################### 4) Differential Expression Analysis

################# create a design matrix 

expression_matrix <- as.matrix(obj_sex_corrected_sums %>% select(-gene))   #|> #select(-contains("smaller"))) 
## converts to expression matrix storing expression values for each sample
rownames(expression_matrix) <- obj_sex_corrected_sums %>% pull(gene)
dim(expression_matrix)
## each row represents a specific gene and its corresponding values in the matrix represent its expression across different samples.
## creates a new df with a single column sampletag containing the sample id's from the expression matrix
design_x <- tibble(sampletag=colnames(expression_matrix)) |> separate(col=sampletag, into =c("days","x"),remove = F, extra="merge") |> 
    separate(col=x, into=c("tissue","virus", "x"), remove = F, extra = "merge") |>
    separate(col=x, into = c("Ms4a3", "gender"), sep = "-") |>
    separate(col=gender, into = c("gender", "small_size"), sep = "_")

design_no_small = design_x |> filter(is.na(small_size))
design_no_small

## plots showing the remaining amounts in the design for the deseq experiment
ggplot(design_x, aes(x = gender)) + 
  geom_histogram(stat = "count") +
  labs(x = "Gender", y = "Count", title = "Gender Distribution") 
ggsave("histogram_deseq_exploratory_designx_gender.png")
ggplot(design_x, aes(x = virus)) + 
  geom_histogram(stat = "count") +
  labs(x = "virus", y = "Count", title = "virus Distribution") 
ggsave("histogram_deseq_exploratory_designx_virus.png")

ggplot(design_x, aes(x = Ms4a3)) + 
  geom_histogram(stat = "count") +
  labs(x = "Ms4a3", y = "Count", title = "Ms4a3 Distribution") 
ggsave("histogram_deseq_exploratory_designx_Ms4a3.png")

ggplot(design_x, aes(x = tissue)) + 
  geom_histogram(stat = "count") +
  labs(x = "tissue", y = "Count", title = "tissue Distribution") 
ggsave("histogram_deseq_exploratory_designx_tissue.png")
ggplot(design_x, aes(x = c(day_size))) + 
  geom_histogram(stat = "count") +
  labs(x = "days", y = "Count", title = "days Distribution") 
ggsave("histogram_deseq_exploratory_design_days.png")

ggplot(design_no_small, aes(x = gender)) + 
  geom_histogram(stat = "count") +
  labs(x = "Gender", y = "Count", title = "Gender Distribution") 
ggsave("histogram_deseq_exploratory_designnosmall_gender.png")
ggplot(design_no_small, aes(x = virus)) + 
  geom_histogram(stat = "count") +
  labs(x = "virus", y = "Count", title = "virus Distribution") 
ggsave("histogram_deseq_exploratory_designnosmall_virus.png")

ggplot(design_no_small, aes(x = Ms4a3)) + 
  geom_histogram(stat = "count") +
  labs(x = "Ms4a3", y = "Count", title = "Ms4a3 Distribution") 
ggsave("histogram_deseq_exploratory_designnosmall_Ms4a3.png")

ggplot(design_no_small, aes(x = tissue)) + 
  geom_histogram(stat = "count") +
  labs(x = "tissue", y = "Count", title = "tissue Distribution") 
ggsave("histogram_deseq_exploratory_designnosmall_tissue.png")
ggplot(design_no_small, aes(x = c(day_size))) + 
  geom_histogram(stat = "count") +
  labs(x = "days", y = "Count", title = "days Distribution") 
ggsave("histogram_deseq_exploratory_designnosmall_days.png")


# this matrix contains information about each sample based on the original sample IDs,
# which are further split to capture experimental factors like "virus", "Ms4a3" status, and "sex" information. 
# The design matrix will be used later in the DESeq2 workflow to model and analyze the relationship between these factors and gene expression levels.
design_x
dim(design_x)


#################################################### 4A) Differential Expression Analysis: First try check the effect if MS4A3 alone

dds_Ms4a3 = DESeqDataSetFromMatrix(expression_matrix,
                                   colData=design_x,
                                   design= ~ 
                                     Ms4a3)
## r log transformation of the deseq2 object to do a variance stabilizing transformation
## improves normality (deseq2 assumes it)
## reduces variance => variance across genes become more stable => imporves stat reliability of tests downstream
## blind parameter: set it to false if interested in DEA because allows the transformation to incorporate the design information, potentially leading to more accurate results
dds_ms4a3_rlogtransformation <- vst(dds_Ms4a3, blind=FALSE)



## create a principal components plot of the deseq2 object after r log transformation
## used to visualize variation between expression analysis samples
##  the PCA plot shows samples subjected to the same condition cluster together
##                     the clusters should be reasonably well-separated along the X-axis (“PC1”: the ﬁrst principal component)
## creates three PCA plots, two grouped by specific factors and one potentially related to virus presence
##

pcaplot_Ms4a3 <- plotPCA(object=dds_ms4a3_rlogtransformation, intgroup=c("sampletag"))

p1_ms4a3 <- pcaplot_Ms4a3+theme_classic()+geom_text(aes(label=group),nudge_y = 0.5)

pcaplot_Ms4a3 <- plotPCA(object=dds_ms4a3_rlogtransformation, intgroup=c("Ms4a3"))

p2_Ms4a3 <- pcaplot_Ms4a3+theme_classic()+geom_text(aes(label=group),nudge_y = 0.5)

pcaplot_Ms4a3 <- plotPCA(object=dds_ms4a3_rlogtransformation, intgroup=c("virus"))

p3_Ms4a3 <- pcaplot_Ms4a3+theme_classic()+geom_text(aes(label=group),nudge_y = 0.5)

pcaplot_Ms4a3 <- plotPCA(object=dds_ms4a3_rlogtransformation, intgroup=c("sex"))

p4_Ms4a3 <- pcaplot_Ms4a3+theme_classic()+geom_text(aes(label=group),nudge_y = 0.5)

pcaplot_Ms4a3 <- plotPCA(object=dds_ms4a3_rlogtransformation, intgroup=c("sex"))

p5_Ms4a3 <- pcaplot_Ms4a3+theme_classic()+geom_text(aes(label=group),nudge_y = 0.5)

p1_ms4a3 ## shows overall plot with full experimental setup
p2_Ms4a3 ## shows only info on ms4a3 status ## clear separation based on ms4a3 status
p3_Ms4a3 ## shows only unfo on type of experiment, 8 or 60 days
p4_Ms4a3 ## shows only info on gender => no clear separation because of sex



### create a heatmap to visualize correlation between expression profiles

dds_ms4a3_rlogtransformation_expressionmatrix <- assay(dds_ms4a3_rlogtransformation) ## extracts the expression matrix
dds_expressionMatrix = assay(dds_Ms4a3)
dds_ms4a3_corrmatrix <- cor(dds_expressionMatrix) ## calculates the corr matrix containing Pearson corr coefficients between all pairs of genes

dds_ms4a3_rlogtransformation_expressionmatrix_corrmatrix <- cor(dds_ms4a3_rlogtransformation_expressionmatrix) ## calculates the corr matrix containing Pearson corr coefficients between all pairs of genes
conds <- (names(obj_sex_corrected_sums)[-1]) ## extracts column names from the starting sums matrix except the first one

paste(design_x$virus,design_x$Ms4a3, sep="_")

design_x <- design_x |> mutate(name=paste(virus,Ms4a3, sep="_")) 
# adds a new column named "name" by applying the previously created string  creates a combined factor based on the existing "virus" and "Ms4a3" factors.
# Plot heatmap
corrHeatMap_Ms4a3 = pheatmap(dds_ms4a3_rlogtransformation_expressionmatrix_corrmatrix)
corrHeatMap_Ms4a3_noTransformation = pheatmap(dds_ms4a3_corrmatrix)
corrHeatMap_Ms4a3 


