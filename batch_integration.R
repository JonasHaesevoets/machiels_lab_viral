
library(SingleR)
library(clustree)
library(tidyseurat)
library(readr)
library(Seurat)
set.seed(2023)

# load and merge data files
exp_2_lung <- "intermediate_data/seurat_obj_experiment_2_lung_afterQCdashboard.rds" |> read_rds()
exp_2_bal <- "intermediate_data/seurat_obj_experiment_2_bal_afterQCdashboard.rds" |> read_rds()

exp_2 <- merge(exp_2_lung, exp_2_bal,add.cell.ids = c("lung", "bal")) 
misc_list_2 <- list(exp_2_lung@misc,exp_2_bal@misc)
exp_2_bal <- NULL
exp_2_lung <- NULL
gc()

exp_1_lung <- "intermediate_data/seurat_obj_experiment_1_lung_afterQCdashboard.rds" |> read_rds()
exp_1_bal <- "intermediate_data/seurat_obj_experiment_1_bal_afterQCdashboard.rds" |> read_rds()

exp_1 <- merge(exp_1_lung, exp_1_bal,add.cell.ids = c("lung", "bal")) 
misc_list_1 <- list(exp_1_lung@misc, exp_1_bal@misc)
exp_1_bal <- NULL
exp_1_lung <- NULL
gc()

obj <- merge(exp_1, exp_2,add.cell.ids = c("exp_1", "exp_2"))
exp_2 <- NULL
exp_1<- NULL
gc()
obj@misc <- list(misc_list_1,misc_list_2)
obj <- JoinLayers(obj, assay = "RNA")
# obj <- JoinLayers(obj, assay = "protein")
# obj <- JoinLayers(obj, assay = "sampletags")
#write_rds(obj," intermediate_data/seurat_obj_experiment_1_2_merged.rds")
#obj <- read_rds(" intermediate_data/seurat_obj_experiment_1_2_merged.rds")

### convert object to v5 object

options(Seurat.object.assay.version = "v5")
meta_data <- obj[[]]


obj.v5 <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = meta_data)
#obj.v5[[]] <- meta_data
#obj.v5$orig.ident <- meta_data$orig.ident
protein <- CreateAssay5Object(counts = obj[["protein"]]$counts)
sampletags <- CreateAssay5Object(counts = obj[["sampletags"]]$counts)
obj.v5[["protein"]] <- protein
obj.v5[["sampletags"]] <- sampletags

obj <- NULL
gc()

obj.v5[["RNA"]] <- split(obj.v5[["RNA"]], f = obj.v5$orig.ident)
gc()

obj.v5 <- NormalizeData(obj.v5) # individual size factors accoriding to Ahlmann-Eltze et al (2023) could be added here
obj.v5 <- FindVariableFeatures(obj.v5,
                               selection.method = "vst",
                               nfeatures = 2000, 
                               verbose = FALSE)
obj.v5 <- ScaleData(obj.v5)
obj.v5 <- RunPCA(obj.v5)



obj.v5 <- IntegrateLayers(
  object = obj.v5, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated.harmony",
  verbose = FALSE
)
obj.v5 <- JoinLayers(obj.v5, assay= "RNA")
#obj <- JoinLayers(obj, assay = "protein")
#obj <- JoinLayers(obj, assay = "sampletags")

obj.v5 <- FindNeighbors(obj.v5, reduction = "integrated.harmony", dims = 1:40)


# Find clusters using a range of resolutions
resolution.range <- seq(from = 0, to = 1, by = 0.2)

obj.v5 <- Seurat::FindClusters(object = obj.v5, resolution = resolution.range)
resolution.range <- seq(from = 0, to = 0.2, by = 0.02)

obj.v5 <- Seurat::FindClusters(object = obj.v5, resolution = resolution.range)
clustree(obj.v5)
ggsave("clustree_integration.png", path = "../../Documents/machiels_lab_viral/CaseStudy_QC_integration/integration", width = 14, height = 14)

# based on this we decided to put the resolution at 0.18

obj.v5 <- FindClusters(obj.v5, cluster.name = "harmony_clusters_0.18", resolution = 0.18)
obj.v5 <- RunUMAP(obj.v5, reduction = "integrated.harmony", reduction.name = "umap_harmony", dims = 1:40)
obj.v5 |> DimPlot(group.by = "harmony_clusters_0.18", label=T)
ggsave("umap_0.18.png", path = "../../Documents/machiels_lab_viral/CaseStudy_QC_integration/integration", width = 12, height = 10)

obj.v5 |>  write_rds("../../Documents/machiels_lab_viral/intermediate_data/seurat_obj_experiment_1_2_integrated.rds")

