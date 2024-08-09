## Differential expression analysis using the central seurat object containing all the clean data ready for analysis
## this script is based on differential_expression_analysis.R
## all the parts that aren't needed are removed (commented out in the original one)
## the main difference is that the central object will be split back up into the two experiments
## we will use seurat objects for both experiments separately and investigate the differential expression



################################################### 1) Load in the required packages and the central seurat object
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("scran")

easypackages::libraries(
  "Seurat","tidyseurat","dplyr","ggplot2","readr","forcats","ragg","ggpubr","stringr","Matrix","textTinyR","tidyverse","ggrepel","scran","pheatmap","DESeq2","patchwork","ggridges", "scran", "textTinyR", "withr", "htmltools","pheatmap", "viridis", "ggExtra", "scales")

obj = read_rds("../../Documents/machiels_lab_viral/intermediate_data/seurat_obj_central.rds")
output = "../../Documents/machiels_lab_viral/"


### split the object in 2 based on the experiment i.e. the amount of days
obj_experiment_list = SplitObject(obj, split.by = "day")

obj_experiment_d60 = obj_experiment_list$d60
obj_experiment_1 = obj_experiment_d60 %>%  filter(str_detect(sample_type, "bal"))
obj_experiment_2 = obj_experiment_d60 %>% filter(str_detect(sample_type, "lung"))

DefaultAssay(obj_experiment_1) <- "RNA"
DefaultAssay(obj_experiment_2) <- "RNA"## make sure that the right assay is used
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

obj_sex_experiment1 <- obj_experiment_1 |>
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


obj_sex_experiment2 <- obj_experiment_2 |>
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

### plots showing the distribution of gender associated genes among the ms4a3 positive and negative controls for the 4 different viral conditions, for experiment 1
plot_distribution_sex_genes_experiment1 <-
  ggplot(data = obj_sex_experiment1,aes(day_sample_type_cond_ms4a3, fill=sex_classifier))+
  geom_bar()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+ ggtitle("sex classification", "female: (Xist>0)AND(max_male_genes==0);male= Xist== AND(max_male_genes>0)")

plot_distribution_sex_genes_experiment1


### create a histogram for the sex genes in question for all conditions (8, 4 viruses with ms4a3 + and -)
theme_x1 <- theme_minimal()+
  theme(strip.text.x = element_text(size = 15))


p1_experiment1 <- obj_sex_experiment1 |> ggplot(aes(Ddx3y)) + geom_histogram(binwidth = 0.1)+facet_wrap( ~day_sample_type_cond_ms4a3, ncol = 1 )+theme_x1
p2_experiment1 <- obj_sex_experiment1 |> ggplot(aes(Eif2s3y))+ geom_histogram(binwidth = 0.1)+facet_wrap( ~day_sample_type_cond_ms4a3, ncol = 1)+theme_x1
p3_experiment1 <- obj_sex_experiment1 |> ggplot(aes(Kdm5d))+ geom_histogram(binwidth = 0.1)+facet_wrap( ~day_sample_type_cond_ms4a3, ncol = 1)+theme_x1
p4_experiment1 <- obj_sex_experiment1 |> ggplot(aes(max_male_genes))+ geom_histogram(binwidth = 0.1)+facet_wrap( ~day_sample_type_cond_ms4a3, ncol = 1)+theme_x1

library(dplyr)
obj_sex_experiment1 |>
  dplyr::select(.cell,day_sample_type_cond_ms4a3,sex_classifier ) |>
  dplyr::mutate(paste0(day_sample_type_cond_ms4a3,sex_classifier, sep="_"))

histogram_sexGenes_allconditions_experiment1 = p1_experiment1+p2_experiment1+p3_experiment1+p4_experiment1
#histogram_sexGenes_allconditions


### create a plot showing the difference between the Xist gene expression and the max male associated genes
plot_difference_xist_maxmalegenes_experiment1 =
  ggplot(data = obj_sex_experiment1,aes(max_male_genes+1, Xist+1)) + geom_point() +ylim(0,3)     +xlim(0,3) +facet_wrap( ~day_sample_type_cond_ms4a3) + ggtitle("Xist vs. max. male gene")
plot_difference_xist_maxmalegenes_experiment1

### create a plots showing the densities of Eif2s3y/Kdm5d/maxmalegene/ across all samples
plot_Eif2s3y_experiment1 <- ggplot(data = obj_sex_experiment1,aes(Eif2s3y, day_sample_type_cond_ms4a3, fill=day_sample_type_cond_ms4a3)) + geom_density_ridges()+ggtitle("Eif2s3y per sample")
plot_Eif2s3y_experiment1
plot_Kdm5d_experiment1 <- ggplot(data = obj_sex_experiment1,aes(Kdm5d, day_sample_type_cond_ms4a3, fill=day_sample_type_cond_ms4a3)) + geom_density_ridges()+ggtitle("Kdm5d per sample")
plot_Kdm5d_experiment1
plot_maxmale_experiment1 <- ggplot(data = obj_sex_experiment1,aes(max_male_genes, day_sample_type_cond_ms4a3, fill=day_sample_type_cond_ms4a3)) + geom_density_ridges()+ggtitle("maxmale per sample")
plot_maxmale_experiment1


### plots showing the distribution of gender associated genes among the ms4a3 positive and negative controls for the 4 different viral conditions, for experiment 2
plot_distribution_sex_genes_experiment2 <-
  ggplot(data = obj_sex_experiment2,aes(day_sample_type_cond_ms4a3, fill=sex_classifier))+
  geom_bar()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+ ggtitle("sex classification", "female: (Xist>0)AND(max_male_genes==0);male= Xist== AND(max_male_genes>0)")

plot_distribution_sex_genes_experiment2


### create a histogram for the sex genes in question for all conditions (8, 4 viruses with ms4a3 + and -)
theme_x1 <- theme_minimal()+
  theme(strip.text.x = element_text(size = 15))


p1_experiment2 <- obj_sex_experiment2 |> ggplot(aes(Ddx3y)) + geom_histogram(binwidth = 0.1)+facet_wrap( ~day_sample_type_cond_ms4a3, ncol = 1 )+theme_x1
p2_experiment2 <- obj_sex_experiment2 |> ggplot(aes(Eif2s3y))+ geom_histogram(binwidth = 0.1)+facet_wrap( ~day_sample_type_cond_ms4a3, ncol = 1)+theme_x1
p3_experiment2 <- obj_sex_experiment2 |> ggplot(aes(Kdm5d))+ geom_histogram(binwidth = 0.1)+facet_wrap( ~day_sample_type_cond_ms4a3, ncol = 1)+theme_x1
p4_experiment2 <- obj_sex_experiment2 |> ggplot(aes(max_male_genes))+ geom_histogram(binwidth = 0.1)+facet_wrap( ~day_sample_type_cond_ms4a3, ncol = 1)+theme_x1

obj_sex_experiment2 |>
  dplyr::select(.cell,day_sample_type_cond_ms4a3,sex_classifier ) |>
  dplyr::mutate(paste0(day_sample_type_cond_ms4a3,sex_classifier, sep="_"))

histogram_sexGenes_allconditions_experiment2 = p1_experiment2+p2_experiment2+p3_experiment2+p4_experiment2
#histogram_sexGenes_allconditions


### create a plot showing the difference between the Xist gene expression and the max male associated genes
plot_difference_xist_maxmalegenes_experiment2 =
  ggplot(data = obj_sex_experiment2,aes(max_male_genes+1, Xist+1)) + geom_point() +ylim(0,3)     +xlim(0,3) +facet_wrap( ~day_sample_type_cond_ms4a3) + ggtitle("Xist vs. max. male gene")
plot_difference_xist_maxmalegenes_experiment2

### create a plots showing the densities of Eif2s3y/Kdm5d/maxmalegene/ across all samples
plot_Eif2s3y_experiment2 <- ggplot(data = obj_sex_experiment2,aes(Eif2s3y, day_sample_type_cond_ms4a3, fill=day_sample_type_cond_ms4a3)) + geom_density_ridges()+ggtitle("Eif2s3y per sample")
plot_Eif2s3y_experiment2
plot_Kdm5d_experiment2 <- ggplot(data = obj_sex_experiment2,aes(Kdm5d, day_sample_type_cond_ms4a3, fill=day_sample_type_cond_ms4a3)) + geom_density_ridges()+ggtitle("Kdm5d per sample")
plot_Kdm5d_experiment2
plot_maxmale_experiment2 <- ggplot(data = obj_sex_experiment2,aes(max_male_genes, day_sample_type_cond_ms4a3, fill=day_sample_type_cond_ms4a3)) + geom_density_ridges()+ggtitle("maxmale per sample")
plot_maxmale_experiment2


################## based on the subsequent plots it was decided to filter by this sex classifier and classify them either as m or f while removing those who are androgenous
keep_experiment1 <- obj_sex_experiment1 |> filter(sex_classifier %in% c("male","female")) |> dplyr::select(.cell,sex_classifier)

obj_sex_corrected_experiment1 <- obj_experiment_1 |> filter(.cell %in% pull(keep_experiment1,.cell)) |> left_join(keep_experiment1) # left join so only those from the left dataframe where a match in the right dataframe are retained
obj_sex_corrected_experiment1 <- obj_sex_corrected_experiment1 |> mutate(virus_ms4a3_sex=paste(day_sample_type_cond_ms4a3_pos_gabbr2,sex_classifier, sep="-"))

obj_sex_corrected_experiment1 <- obj_sex_corrected_experiment1 %>%
  filter(!str_detect(day_sample_type_cond_ms4a3_pos_gabbr2, "Gabbr2_pos_Ms4a3_neg")) ## remove those positive for Gabbr2 but negative for Ms4a3

obj_sex_corrected_experiment1_sums <- make_idents_counts(ds=obj_sex_corrected_experiment1,idents_choose = "virus_ms4a3_sex",perc_expressed = 0.05) ## call the functions specified on the ds with the specified idents and a min percentage that a gene should be expressed


keep_experiment2 <- obj_sex_experiment2 |> filter(sex_classifier %in% c("male","female")) |> dplyr::select(.cell,sex_classifier)

obj_sex_corrected_experiment2 <- obj_experiment_2|> filter(.cell %in% pull(keep_experiment2,.cell)) |> left_join(keep_experiment2) # left join so only those from the left dataframe where a match in the right dataframe are retained
obj_sex_corrected_experiment2 <- obj_sex_corrected_experiment2 |> mutate(virus_ms4a3_sex=paste(day_sample_type_cond_ms4a3_pos_gabbr2,sex_classifier, sep="-"))

obj_sex_corrected_experiment2 <- obj_sex_corrected_experiment2 %>%
  filter(!str_detect(day_sample_type_cond_ms4a3_pos_gabbr2, "Gabbr2_pos_Ms4a3_neg")) ## remove those positive for Gabbr2 but negative for Ms4a3

obj_sex_corrected_experiment2_sums <- make_idents_counts(ds=obj_sex_corrected_experiment2,idents_choose = "virus_ms4a3_sex",perc_expressed = 0.05) ## call the functions specified on the ds with the specified idents and a min percentage that a gene should be expressed




#################################################### 4) Differential Expression Analysis for experiment 1

################# create a design matrix 

expression_matrix_experiment1 <- as.matrix(obj_sex_corrected_experiment1_sums %>% dplyr::select(-gene))   #|> #select(-contains("smaller"))
## converts to expression matrix storing expression values for each sample
rownames(expression_matrix_experiment1) <- obj_sex_corrected_experiment1_sums %>% dplyr::pull(gene)
dim(expression_matrix_experiment1)
## each row represents a specific gene and its corresponding values in the matrix represent its expression across different samples.
## creates a new df with a single column sampletag containing the sample id's from the expression matrix
design_x_experiment1 <- tibble(sampletag=colnames(expression_matrix_experiment1)) |> separate(col=sampletag, into =c("days","x"),remove = F, extra="merge") |> 
  separate(col=x, into=c("tissue","virus", "x"), remove = F, extra = "merge") |>
  separate(col=x, into = c("Ms4a3", "gender"), sep = "-") |>
  separate(col=gender, into = c("gender", "small_size"), sep = "_")


## design no small will be used downstreal
## the small size parameter corresponds to those samples who don't exist or those where the sequencing was not deep enough
# this matrix contains information about each sample based on the original sample IDs,
# which are further split to capture experimental factors like "virus", "Ms4a3" status, and "sex" information. 
# The design matrix will be used later in the DESeq2 workflow to model and analyze the relationship between these factors and gene expression levels.

design_no_small_experiment1 = design_x_experiment1 |> filter(is.na(small_size))
design_no_small_experiment1

## because we the amount of columns in the expression matrix has to be equal to the amount of rows in the design matrix 
## need to filter the columns where the rowSums are 0 in the expression matrix to have it correctly
keep_cols_experiment1 <- colSums(abs(expression_matrix_experiment1)) > 0
expression_matrix_experiment1 = expression_matrix_experiment1[, keep_cols_experiment1]

## plots showing the remaining amounts in the design for the deseq experiment
histogram_experiment1_exploration_deseq_gender = ggplot(design_no_small_experiment1, aes(x = gender)) + 
  geom_histogram(stat = "count") +
  labs(x = "Gender", y = "Count", title = "Gender Distribution") 
histogram_experiment1_exploration_deseq_gender

histogram_experiment1_exploration_deseq_virus =ggplot(design_no_small_experiment1, aes(x = virus)) + 
  geom_histogram(stat = "count") +
  labs(x = "virus", y = "Count", title = "virus Distribution") 
histogram_experiment1_exploration_deseq_virus

histogram_experiment1_exploration_deseq_Ms4a3 = ggplot(design_no_small_experiment1, aes(x = Ms4a3)) + 
  geom_histogram(stat = "count") +
  labs(x = "Ms4a3", y = "Count", title = "Ms4a3 Distribution") 
histogram_experiment1_exploration_deseq_Ms4a3

histogram_experiment1_exploration_deseq_tissue =ggplot(design_no_small_experiment1, aes(x = tissue)) + 
  geom_histogram(stat = "count") +
  labs(x = "tissue", y = "Count", title = "tissue Distribution") 
histogram_experiment1_exploration_deseq_tissue




#################################################### 4A) Differential Expression Analysis: based on the full model with tissue

dds_full_experiment1_tissue = DESeqDataSetFromMatrix(expression_matrix_experiment1,
                                   colData=design_no_small_experiment1,
                                   design= ~ virus + Ms4a3)
## in this model we estimate the effect the effect of tissue after adjusting for gender virus and ms4a3



## r log transformation of the deseq2 object to do a variance stabilizing transformation
## improves normality (deseq2 assumes it)
## reduces variance => variance across genes become more stable => imporves stat reliability of tests downstream
## blind parameter: set it to false if interested in DEA because allows the transformation to incorporate the design information, potentially leading to more accurate results
dds_full_experiment1_tissue_rlogtransformation <- rlog(dds_full_experiment1_tissue, blind=FALSE)





## create a principal components plot of the deseq2 object after r log transformation
## used to visualize variation between expression analysis samples
##  the PCA plot shows samples subjected to the same condition cluster together
##                     the clusters should be reasonably well-separated along the X-axis (“PC1”: the ﬁrst principal component)
## creates three PCA plots, two grouped by specific factors and one potentially related to virus presence
##

pcaplot_full_experiment1_tissue <- plotPCA(object=dds_full_experiment1_tissue_rlogtransformation, intgroup=c("sampletag"))
plot_pca_sampletag_full_experiment1_tissue <- pcaplot_full_experiment1_tissue+theme_classic()+geom_text(aes(label=group),nudge_y = 0.5)
ggsave("pca_d60_bal_tissue.png", width = 10, height = 10)

pcaplot_full_experiment1_tissue <- plotPCA(object=dds_full_experiment1_tissue_rlogtransformation, intgroup=c("Ms4a3"))
plot_pca_Ms4a3_full_experiment1_tissue <- pcaplot_full_experiment1_tissue+theme_classic()+geom_text(aes(label=group),nudge_y = 0.5)
ggsave("pca_d60_bal_Ms4a3.png", width = 10, height = 10)

pcaplot_full_experiment1_tissue <- plotPCA(object=dds_full_experiment1_tissue_rlogtransformation, intgroup=c("virus"))
plot_pca_virus_full_experiment1_tissue <- pcaplot_full_experiment1_tissue+theme_classic()+geom_text(aes(label=group),nudge_y = 0.5)
ggsave("pca_d60_bal_condition.png", width = 10, height = 10)

# pcaplot_full_experiment1_tissue <- plotPCA(object=dds_full_experiment1_tissue_rlogtransformation, intgroup=c("gender"))
# plot_pca_gender_full_experiment1_tissue <- pcaplot_full_experiment1_tissue+theme_classic()+geom_text(aes(label=group),nudge_y = 0.5)
# ggsave("pca_d60_bal_gender.png")

# pcaplot_full_experiment1_tissue <- plotPCA(object=dds_full_experiment1_tissue_rlogtransformation, intgroup=c("tissue"))
# plot_pca_tissue_full_experiment1_tissue <- pcaplot_full_experiment1_tissue+theme_classic()+geom_text(aes(label=group),nudge_y = 0.5)
# ggsave("pca_d60_bal_tissuee.png")
plot_pca_sampletag_full_experiment1_tissue ## shows overall plot with full experimental setup
### 3 clusters 2 for the lung
plot_pca_Ms4a3_full_experiment1_tissue ## shows only info on ms4a3 status ## no clear separation
plot_pca_virus_full_experiment1_tissue ## shows only unfo on type of virus => no clear separation
# plot_pca_gender_full_experiment1_tissue ## shows only info on gender => no clear separation because of sex
# plot_pca_tissue_full_experiment1_tissue  ## shows only info on tissue => clear separation between bal and lung











#################################################### 4) Differential Expression Analysis for experiment 1

################# create a design matrix 

expression_matrix_experiment2 <- as.matrix(obj_sex_corrected_experiment2_sums %>% dplyr::select(-gene))   #|> #select(-contains("smaller"))
## converts to expression matrix storing expression values for each sample
rownames(expression_matrix_experiment2) <- obj_sex_corrected_experiment2_sums %>% dplyr::pull(gene)
dim(expression_matrix_experiment2)
## each row represents a specific gene and its corresponding values in the matrix represent its expression across different samples.
## creates a new df with a single column sampletag containing the sample id's from the expression matrix
design_x_experiment2 <- tibble(sampletag=colnames(expression_matrix_experiment2)) |> separate(col=sampletag, into =c("days","x"),remove = F, extra="merge") |> 
  separate(col=x, into=c("tissue","virus", "x"), remove = F, extra = "merge") |>
  separate(col=x, into = c("Ms4a3", "gender"), sep = "-") |>
  separate(col=gender, into = c("gender", "small_size"), sep = "_")


## design no small will be used downstreal
## the small size parameter corresponds to those samples who don't exist or those where the sequencing was not deep enough
# this matrix contains information about each sample based on the original sample IDs,
# which are further split to capture experimental factors like "virus", "Ms4a3" status, and "sex" information. 
# The design matrix will be used later in the DESeq2 workflow to model and analyze the relationship between these factors and gene expression levels.

design_no_small_experiment2 = design_x_experiment2 |> filter(is.na(small_size))
design_no_small_experiment2

## because we the amount of columns in the expression matrix has to be equal to the amount of rows in the design matrix 
## need to filter the columns where the rowSums are 0 in the expression matrix to have it correctly
keep_cols_experiment2 <- colSums(abs(expression_matrix_experiment2)) > 0
expression_matrix_experiment2 = expression_matrix_experiment2[, keep_cols_experiment2]

## plots showing the remaining amounts in the design for the deseq experiment
histogram_experiment2_exploration_deseq_gender = ggplot(design_no_small_experiment2, aes(x = gender)) + 
  geom_histogram(stat = "count") +
  labs(x = "Gender", y = "Count", title = "Gender Distribution") 
histogram_experiment2_exploration_deseq_gender

histogram_experiment2_exploration_deseq_virus =ggplot(design_no_small_experiment2, aes(x = virus)) + 
  geom_histogram(stat = "count") +
  labs(x = "virus", y = "Count", title = "virus Distribution") 
histogram_experiment2_exploration_deseq_virus

histogram_experiment2_exploration_deseq_Ms4a3 = ggplot(design_no_small_experiment2, aes(x = Ms4a3)) + 
  geom_histogram(stat = "count") +
  labs(x = "Ms4a3", y = "Count", title = "Ms4a3 Distribution") 
histogram_experiment2_exploration_deseq_Ms4a3

histogram_experiment2_exploration_deseq_tissue =ggplot(design_no_small_experiment2, aes(x = tissue)) + 
  geom_histogram(stat = "count") +
  labs(x = "tissue", y = "Count", title = "tissue Distribution") 
histogram_experiment2_exploration_deseq_tissue




#################################################### 4A) Differential Expression Analysis: based on the full model with tissue

dds_full_experiment2_tissue = DESeqDataSetFromMatrix(expression_matrix_experiment2,
                                                     colData=design_no_small_experiment2,
                                                     design= ~ virus + Ms4a3)
## in this model we estimate the effect the effect of tissue after adjusting for gender virus and ms4a3



## r log transformation of the deseq2 object to do a variance stabilizing transformation
## improves normality (deseq2 assumes it)
## reduces variance => variance across genes become more stable => imporves stat reliability of tests downstream
## blind parameter: set it to false if interested in DEA because allows the transformation to incorporate the design information, potentially leading to more accurate results
dds_full_experiment2_tissue_rlogtransformation <- rlog(dds_full_experiment2_tissue, blind=FALSE)





## create a principal components plot of the deseq2 object after r log transformation
## used to visualize variation between expression analysis samples
##  the PCA plot shows samples subjected to the same condition cluster together
##                     the clusters should be reasonably well-separated along the X-axis (“PC1”: the ﬁrst principal component)
## creates three PCA plots, two grouped by specific factors and one potentially related to virus presence
##

pcaplot_full_experiment2_tissue <- plotPCA(object=dds_full_experiment2_tissue_rlogtransformation, intgroup=c("sampletag"))

plot_pca_sampletag_full_experiment2_tissue <- pcaplot_full_experiment2_tissue+theme_classic()+geom_text(aes(label=group),nudge_y = 0.5)
ggsave("pca_d60_lung_tissue.png",  width = 10, height = 10)
pcaplot_full_experiment2_tissue <- plotPCA(object=dds_full_experiment2_tissue_rlogtransformation, intgroup=c("Ms4a3"))

plot_pca_Ms4a3_full_experiment2_tissue <- pcaplot_full_experiment2_tissue+theme_classic()+geom_text(aes(label=group),nudge_y = 0.5)
ggsave("pca_d60_lung_Ms4a3.png", width = 10, height = 10)
pcaplot_full_experiment2_tissue <- plotPCA(object=dds_full_experiment2_tissue_rlogtransformation, intgroup=c("virus"))

plot_pca_virus_full_experiment2_tissue <- pcaplot_full_experiment2_tissue+theme_classic()+geom_text(aes(label=group),nudge_y = 0.5)
ggsave("pca_d60_lung_condition.png",  width = 10, height = 10)
# 
# pcaplot_full_experiment2_tissue <- plotPCA(object=dds_full_experiment2_tissue_rlogtransformation, intgroup=c("gender"))
# 
# plot_pca_gender_full_experiment2_tissue <- pcaplot_full_experiment2_tissue+theme_classic()+geom_text(aes(label=group),nudge_y = 0.5)
# ggsave("pca_d60_lung_gender.png")

# pcaplot_full_experiment2_tissue <- plotPCA(object=dds_full_experiment2_tissue_rlogtransformation, intgroup=c("tissue"))
# 
# plot_pca_tissue_full_experiment2_tissue <- pcaplot_full_experiment2_tissue+theme_classic()+geom_text(aes(label=group),nudge_y = 0.5)
# ggsave("pca_d60_lung_tisssue.png")

plot_pca_sampletag_full_experiment2_tissue ## shows overall plot with full experimental setup
### 3 clusters 2 for the lung
plot_pca_Ms4a3_full_experiment2_tissue ## shows only info on ms4a3 status ## no clear separation
plot_pca_virus_full_experiment2_tissue ## shows only unfo on type of virus => no clear separation
# plot_pca_tissue_full_experiment2_tissue 




















# ### create a heatmap to visualize correlation between expression profiles
# 
# dds_full_rlogtransformation_expressionmatrix_experiment1 <- assay(dds_full_experiment1_tissue_rlogtransformation) ## extracts the expression matrix
# dds_expressionMatrix_experiment1 = assay(dds_full_experiment1_tissue)
# dds_full_corrmatrix_experiment1 <- cor(dds_expressionMatrix_experiment1) ## calculates the corr matrix containing Pearson corr coefficients between all pairs of genes
# 
# dds_full_rlogtransformation_expressionmatrix_corrmatrix_experiment1 <- cor(dds_full_rlogtransformation_expressionmatrix_experiment1) ## calculates the corr matrix containing Pearson corr coefficients between all pairs of genes
# conds <- (names(obj_sex_corrected_experiment1)[-1]) ## extracts column names from the starting sums matrix except the first one
# 
# 
# design_no_small_experiment1 <- design_no_small_experiment1 |> mutate(name=paste(virus,Ms4a3, sep="_")) 
# # adds a new column named "name" by applying the previously created string  creates a combined factor based on the existing "virus" and "Ms4a3" factors.
# #+ Plot heatmap
# plot_corrHeatMap_full_experiment1 = pheatmap(dds_full_rlogtransformation_expressionmatrix_corrmatrix_experiment1)
# print(plot_corrHeatMap_full_experiment1)
# 
# 
# ##### run the deseq analysis
# dds_full_experiment1_tissue_deseq = DESeq(dds_full_experiment1_tissue)
# # return the available comparisons
# resultsNames(dds_full_experiment1_tissue_deseq)
# # test_no_gender <- as_tibble(results(DESeq(dds_full_experiment1_tissue_deseq, test="LRT", reduced= ~tissue + virus + Ms4a3)))
# # test_no_gender <- as_tibble(results(DESeq(dds_full_experiment1_tissue_deseq, test="LRT", reduced= ~tissue + virus + Ms4a3)))
# 
# 
# 
# #+
# #### export the normalized counts
# dds_full_experiment1_tissue_deseq_normalized = counts(dds_full_experiment1_tissue_deseq, normalized = T)
# write.csv(dds_full_experiment1_tissue_deseq_normalized, "C:/Users/Jonas/Desktop/Analysis_Lucia/R_Projects/machiels_lab_viral-main/output/deseq/normalizedCountsExperiment1.csv")
# 
# 
# #+
# #### size factors
# #### higher numbers tend to correlate with higher sequencing depths (depth = average number of times a given position in the genetic sequence is sequenced)
# sizeFactors(dds_full_experiment1_tissue_deseq)
# 
# ### IMPORTANT
# ### Genes that are highly expressed will have a more consistent level of variations, but it will be higher than the mean.
# ### Lowly expressed genes will exhibit variation that hovers around the mean (but with a higher amount of variability).
# 
# 
# #### check the model fit by looking at the plot of dispersion estimates obtained from running DEseq2() on the deseqdataset object
# #### dispersion models within group variability by describing how much variance deviates from mean
# #### dispersion estimates reflecr variance in gene expression for a given mean value
# #### red line is drawn using info across all genes to get better estimates => this principle is called shrinking
# plot_dispersion_full_experiment1 = plotDispEsts(dds_full_experiment1_tissue_deseq)
# # evaluation of the plot: data scatters well around the curve                       LOOKS FINE
# #                         not too strong shrinkage
# #                         decrease in dispersion with increasing mean expression
# 
# 
# 
# 
# #+
# ##### Deseq results
# 
# # returns results from output showing the differentially expressed genes
# # alpha was set to 0.05 meaning the adjusted p value will be 0.05
# # only those genes with an adjusted p value below 0.05 are termed as differentially expressed genes
# # contains results of every stat test for every single row (gene)
# # wald test performed to obtain the output comparing bal vs lung tissue 
# dds_full_experiment1_tissue_results_contrastTissue = results(dds_full_experiment1_tissue_deseq, alpha = 0.01, contrast = c("tissue", "lung", "bal"))
# # baseMean is mean expression of a gene across all samples
# # p value is obtained using the wald test using H0 saying there is no differential between the groups 
# # if the p value is below .01 the H0 gets rejected because only a 1% chance that H0 is true
# # however in cases like this when 13801 genes are tested
# # by chance there is a 1% chance that genes are not differentially expressed but obtained significant p values => false positives
# # hence why the adjusted p value is used over the ordinary p value
# # log2foldchange (log2FC) refers to the magnitude and direction of gene expression differences between groups being compared
# # A positive log2FC indicates that a gene is upregulated in one group compared to the other. The higher the value, the greater the upregulation.
# # Conversely, a negative value signifies downregulation, with a more negative value indicating stronger downregulation.
# # in this case a positive log2fc score indicates biased expression towards the lung
# # in this case a negative log2fc score indicates biased expression towards the bal
# summary(dds_full_experiment1_tissue_results_contrastTissue)
# # we obtain a very high amount of genes with significant p adjusted values for both cases
# # 5059 genes are upregulated in the lung 
# # 3860 genes are upregulated in the bal
# # this is a very high amount and will be filtered upon to further refine the results to obtain the differentially expressed genes in lung and bal
# 
# 
# #+
# # #### create a MA plot
# # #### displays a log ratio (M) vs an average (A) in order to visualize the differences between two groups
# # #### Visualize the relationship between gene expression and fold change between two conditions or groups being compared.
# # #### Identify potential DEGs based on their location on the plot
# # #### x-axis(A): Represents the average expression level (mean) of a gene across samples in both conditions being compared. 
# # ####            This indicates the overall expression level of the gene.
# # #### y-axis(M): Represents the log2 fold change (M) in gene expression between the two conditions. + => upreg -=> downreg
# # #### ideally majority genes near origin, genes far from origin on y axis are potential deg's
# # 
# plot_MA_full_experiment1 = plotMA(dds_full_experiment1_tissue_results_contrastTissue)
# plot_MA_full_experiment1
# 
# #+
# #### make the MA better
# # Coerce to a data frame
# dds_full_experiment1_tissue_deseq_contrastTissue_df <- as.data.frame(dds_full_experiment1_tissue_results_contrastTissue)
# 
# #Examine this data frame
# head(dds_full_experiment1_tissue_deseq_contrastTissue_df)
# #
# # Set a boolean column for significance based on the adjusted p value and log2foldChange
# dds_full_experiment1_tissue_deseq_contrastTissue_df$significant <- ifelse(dds_full_experiment1_tissue_deseq_contrastTissue_df$padj < 0.1 & abs(dds_full_experiment1_tissue_deseq_contrastTissue_df$log2FoldChange) > 1,"Significant", NA)
# # Plot the results similar to DEseq2
# ggplot(dds_full_experiment1_tissue_deseq_contrastTissue_df, aes(baseMean, log2FoldChange, colour=significant)) +
#     geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10()+
#     geom_hline(yintercept = 0, colour="tomato1", size=2) + labs(x="mean of normalized counts", y="log fold change") +
#     scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + theme_bw()
# 
# # Let's add some more detail
# plot_MA_full_experiment1_detail = ggplot(dds_full_experiment1_tissue_deseq_contrastTissue_df, aes(baseMean, log2FoldChange, colour=padj)) +
#     geom_point(size=1) + scale_y_continuous(limits=c(-3, 5), oob=squish) + scale_x_log10() +
#     geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash")+
#     labs(x="mean of normalized counts", y="log fold change")
# plot_MA_full_experiment1_detail
# 
# #+ create a volcano plot
# library(EnhancedVolcano)
# 
# plot_volcano_full_experiment1_tissue_contrastTissue = EnhancedVolcano(dds_full_experiment1_tissue_results_contrastTissue,
#                                                                       lab = rownames(dds_full_experiment1_tissue_results_contrastTissue),
#                                                                       x = 'log2FoldChange',
#                                                                       y = 'pvalue',
#                                                                       title = 'lung vs ball using full model')
# plot_volcano_full_experiment1_tissue_contrastTissue
# 





