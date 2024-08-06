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

obj = read_rds("../../Documents/machiels_lab_viral/intermediate_data/seurat_obj_central_am.rds")
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



################## based on the subsequent plots it was decided to filter by this sex classifier and classify them either as m or f while removing those who are androgenous


obj_corrected_experiment1 <- obj_experiment_1 %>%
  filter(!str_detect(day_sample_type_cond_ms4a3_pos_gabbr2, "Gabbr2_pos_Ms4a3_neg")) ## remove those positive for Gabbr2 but negative for Ms4a3

obj_corrected_experiment1_sums <- make_idents_counts(ds=obj_corrected_experiment1,idents_choose = "sampletag_name",perc_expressed = 0.05) ## call the functions specified on the ds with the specified idents and a min percentage that a gene should be expressed




obj_corrected_experiment2 <- obj_experiment_2 %>%
  filter(!str_detect(day_sample_type_cond_ms4a3_pos_gabbr2, "Gabbr2_pos_Ms4a3_neg")) ## remove those positive for Gabbr2 but negative for Ms4a3

obj_corrected_experiment2_sums <- make_idents_counts(ds=obj_corrected_experiment2,idents_choose = "sampletag_name",perc_expressed = 0.05) ## call the functions specified on the ds with the specified idents and a min percentage that a gene should be expressed




#################################################### 4) Differential Expression Analysis for experiment 1

################# create a design matrix 

expression_matrix_experiment1 <- as.matrix(obj_corrected_experiment1_sums %>% dplyr::select(-gene))   #|> #select(-contains("smaller"))
## converts to expression matrix storing expression values for each sample
rownames(expression_matrix_experiment1) <- obj_corrected_experiment1_sums %>% dplyr::pull(gene)
dim(expression_matrix_experiment1)
## each row represents a specific gene and its corresponding values in the matrix represent its expression across different samples.
## creates a new df with a single column sampletag containing the sample id's from the expression matrix
design_x_experiment1 <- tibble(sampletag=colnames(expression_matrix_experiment1)) |> separate(col=sampletag, into =c("condition","Ms4a3"),remove = F, extra="merge")
  





keep_cols_experiment1 <- colSums(abs(expression_matrix_experiment1)) > 0
expression_matrix_experiment1 = expression_matrix_experiment1[, keep_cols_experiment1]



histogram_experiment1_exploration_deseq_virus =ggplot(design_x_experiment1, aes(x = condition)) + 
  geom_histogram(stat = "count") +
  labs(x = "virus", y = "Count", title = "virus Distribution") 
histogram_experiment1_exploration_deseq_virus

histogram_experiment1_exploration_deseq_Ms4a3 = ggplot(design_x_experiment1, aes(x = Ms4a3)) + 
  geom_histogram(stat = "count") +
  labs(x = "Ms4a3", y = "Count", title = "Ms4a3 Distribution") 
histogram_experiment1_exploration_deseq_Ms4a3





#################################################### 4A) Differential Expression Analysis: based on the full model with tissue

dds_full_experiment1_tissue = DESeqDataSetFromMatrix(expression_matrix_experiment1,
                                                     colData=design_x_experiment1,
                                                     design= ~ condition + Ms4a3)

dds_full_experiment1_tissue_rlogtransformation <- rlog(dds_full_experiment1_tissue, blind=FALSE)







pcaplot_full_experiment1_tissue <- plotPCA(object=dds_full_experiment1_tissue_rlogtransformation, intgroup=c("sampletag"))
plot_pca_sampletag_full_experiment1_tissue <- pcaplot_full_experiment1_tissue+theme_classic()+geom_text(aes(label=group),nudge_y = 0.5)
ggsave("pca_d60_am_bal_sampletag.png", width = 10, height = 10)

pcaplot_full_experiment1_tissue <- plotPCA(object=dds_full_experiment1_tissue_rlogtransformation, intgroup=c("Ms4a3"))
plot_pca_Ms4a3_full_experiment1_tissue <- pcaplot_full_experiment1_tissue+theme_classic()+geom_text(aes(label=group),nudge_y = 0.5)
ggsave("pca_d60_am_bal_Ms4a3.png", width = 10, height = 10)

pcaplot_full_experiment1_tissue <- plotPCA(object=dds_full_experiment1_tissue_rlogtransformation, intgroup=c("condition"))
plot_pca_virus_full_experiment1_tissue <- pcaplot_full_experiment1_tissue+theme_classic()+geom_text(aes(label=group),nudge_y = 0.5)
ggsave("pca_d60_am_bal_condition.png", width = 10, height = 10)


plot_pca_sampletag_full_experiment1_tissue
plot_pca_Ms4a3_full_experiment1_tissue 
plot_pca_virus_full_experiment1_tissue 











#################################################### 4) Differential Expression Analysis for experiment 1

################# create a design matrix 

expression_matrix_experiment2 <- as.matrix(obj_corrected_experiment2_sums %>% dplyr::select(-gene))   #|> #select(-contains("smaller"))
## converts to expression matrix storing expression values for each sample
rownames(expression_matrix_experiment2) <- obj_corrected_experiment2_sums %>% dplyr::pull(gene)
dim(expression_matrix_experiment2)
## each row represents a specific gene and its corresponding values in the matrix represent its expression across different samples.
## creates a new df with a single column sampletag containing the sample id's from the expression matrix
design_x_experiment2 <- tibble(sampletag=colnames(expression_matrix_experiment2)) |> separate(col=sampletag, into =c("condition","Ms4a3"),remove = F, extra="merge")



## design no small will be used downstreal
## the small size parameter corresponds to those samples who don't exist or those where the sequencing was not deep enough
# this matrix contains information about each sample based on the original sample IDs,
# which are further split to capture experimental factors like "virus", "Ms4a3" status, and "sex" information. 
# The design matrix will be used later in the DESeq2 workflow to model and analyze the relationship between these factors and gene expression levels.



## because we the amount of columns in the expression matrix has to be equal to the amount of rows in the design matrix 
## need to filter the columns where the rowSums are 0 in the expression matrix to have it correctly
keep_cols_experiment2 <- colSums(abs(expression_matrix_experiment2)) > 0
expression_matrix_experiment2 = expression_matrix_experiment2[, keep_cols_experiment2]

## plots showing the remaining amounts in the design for the deseq experiment

histogram_experiment2_exploration_deseq_virus =ggplot(design_x_experiment2, aes(x = condition)) + 
  geom_histogram(stat = "count") +
  labs(x = "virus", y = "Count", title = "virus Distribution") 
histogram_experiment2_exploration_deseq_virus

histogram_experiment2_exploration_deseq_Ms4a3 = ggplot(design_x_experiment2, aes(x = Ms4a3)) + 
  geom_histogram(stat = "count") +
  labs(x = "Ms4a3", y = "Count", title = "Ms4a3 Distribution") 
histogram_experiment2_exploration_deseq_Ms4a3





#################################################### 4A) Differential Expression Analysis: based on the full model with tissue

dds_full_experiment2_tissue = DESeqDataSetFromMatrix(expression_matrix_experiment2,
                                                     colData=design_x_experiment2,
                                                     design= ~ condition + Ms4a3)




dds_full_experiment2_tissue_rlogtransformation <- rlog(dds_full_experiment2_tissue, blind=FALSE)





## create a principal components plot of the deseq2 object after r log transformation
## used to visualize variation between expression analysis samples
##  the PCA plot shows samples subjected to the same condition cluster together
##                     the clusters should be reasonably well-separated along the X-axis (“PC1”: the ﬁrst principal component)
## creates three PCA plots, two grouped by specific factors and one potentially related to virus presence
##

pcaplot_full_experiment2_tissue <- plotPCA(object=dds_full_experiment2_tissue_rlogtransformation, intgroup=c("sampletag"))

plot_pca_sampletag_full_experiment2_tissue <- pcaplot_full_experiment2_tissue+theme_classic()+geom_text(aes(label=group),nudge_y = 0.5)
ggsave("pca_d60_am_lung_sampletag.png",  width = 10, height = 10)
pcaplot_full_experiment2_tissue <- plotPCA(object=dds_full_experiment2_tissue_rlogtransformation, intgroup=c("Ms4a3"))

plot_pca_Ms4a3_full_experiment2_tissue <- pcaplot_full_experiment2_tissue+theme_classic()+geom_text(aes(label=group),nudge_y = 0.5)
ggsave("pca_d60_am_lung_Ms4a3.png", width = 10, height = 10)
pcaplot_full_experiment2_tissue <- plotPCA(object=dds_full_experiment2_tissue_rlogtransformation, intgroup=c("condition"))

plot_pca_virus_full_experiment2_tissue <- pcaplot_full_experiment2_tissue+theme_classic()+geom_text(aes(label=group),nudge_y = 0.5)
ggsave("pca_d60_am_lung_condition.png",  width = 10, height = 10)


plot_pca_sampletag_full_experiment2_tissue 
plot_pca_Ms4a3_full_experiment2_tissue 
plot_pca_virus_full_experiment2_tissue 



















