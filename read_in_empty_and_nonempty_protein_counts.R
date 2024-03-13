

library(data.table)
library(R.utils)
library(readr)
unfiltered_csv_paths <- c( "../../Documents/machiels_lab_viral/raw_data\\machiels_lab\\viral\\output_lung_d8\\BD-Analysis-BMachiels-Lung_DBEC_MolsPerCell_Unfiltered.csv.gz",
                           "../../Documents/machiels_lab_viral/raw_data\\machiels_lab\\viral\\output_bal_d8\\BD-Analysis-BMachiels-BAL_DBEC_MolsPerCell_Unfiltered.csv.gz",
                           "../../Documents/machiels_lab_viral/raw_data\\machiels_lab\\viral\\2023-10-02_output_lung\\Output_Lung\\BD-Analysis-BMachiels_DBEC_MolsPerCell_Unfiltered.csv.gz",
                           "../../Documents/machiels_lab_viral/raw_data\\machiels_lab\\viral\\2023-10-02_output_bal\\Output_BAL\\BD-Analysis-BMachiels_DBEC_MolsPerCell_Unfiltered.csv.gz")



exp_name <- c("lung_d8","bal_d8",
              "lung_d60", "bal_d60")

for (i in seq_along(unfiltered_csv_paths)) {
        counts_fread = fread(unfiltered_csv_paths[i], select = c("Cell_Index",
                                                                 "Siglec-F|Siglecf|AMM2013|pAbO",
                                                                 "I-A_I-E|H2-Ab_Ad_Aq_Ed_Ek|AMM2019|pAbO",
                                                                 "CD274|Cd274|AMM2038|pAbO",
                                                                 "CD11c:HL3|Itgax|AMM2008|pAbO",
                                                                 "Ly-6G|Ly6g|AMM2009|pAbO",
                                                                 "Ly-6A_Ly-6E|Ly6a_Ly6e|AMM2026|pAbO"), 
                             showProgress = T)
        
        colnames(counts_fread) <- c("cell_index","Siglecf_AbSeq","H2_ia_ie_AbSeq","Cd274_AbSeq","Cd11c","Ly_6g_AbSeq","Ly_6a_AbSeq")
        write_csv(counts_fread,paste0("../../Documents/machiels_lab_viral/intermediate_data/",exp_name[i],"unfiltered_prot_counts_fread.csv") )
        
}


