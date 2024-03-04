
library(shiny)
library(dplyr)
library(ggplot2)
library(readr)
library(forcats)
#library(tidyseurat)
# bal_alv_mac <- read_rds( "C:\\Users\\danne\\R_projects\\machiels_lab_viral\\intermediate_data\\seurat_obj_experiment_1_bal_alv_macs_with_merged_groups.rds")
# write_rds(bal_alv_mac, "app_data\\")
# bal_alv_mac$virus <- tibble(sampletag_name=bal_alv_mac |> pull(sampletag_name)) |> separate(sampletag_name,into = c("virus")) |> pull("virus")

# bootstrap_mean_tbl <- read_rds("C:\\Users\\danne\\R_projects\\machiels_lab_viral\\intermediate_data\\bootstrap_mean_ddfc_all_cond_downsample_size_50_genes.rds")
# bootstrap_mean_tbl |> write_rds("C:\\Users\\danne\\R_projects\\machiels_lab_viral\\boostrap_fc_fc\\app_data\\bootstrap_mean_tbl.rds")
bootstrap_mean_tbl <- read_rds("bootstrap_mean_tbl.rds")
gene_vec <- bootstrap_mean_tbl$gene |> unique()
condition_vec <- bootstrap_mean_tbl$condition |> unique()
sample_size_vec <- bootstrap_mean_tbl$group_size |> unique()



 fc_plot <-function(dat, x,y,contrast)
    
{         ggplot(data = dat,aes({{x}},{{y}},color={{contrast}}))+
        geom_density_2d(contour_var = "ndensity") +  
        geom_abline(intercept = 0, slope = 1, color="black", linetype="dashed",linewidth=1)+
        geom_vline(xintercept =  0,  color="black")+
        geom_hline(yintercept =  0,  color="black")+
        xlim(c(-2,2))+ylim(c(-2,2))+
        theme_bw()
}


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Bootstrap foldchange comparison of \n
               differentially expressed genes Msa4a3+ vs Msa4a3"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectizeInput(
                'genes_selected', 'gene multi-select', choices = gene_vec, multiple = TRUE
            ),
            selectizeInput(
                'condition_x', 'condition x vs. Mock', choices = condition_vec, multiple = FALSE
            ),
            selectizeInput(
                'sample_size', 'sample size', choices = sample_size_vec, multiple = FALSE
            )
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("bootstrap_fc_fc_plot"),
           plotOutput("bar_plot")
        )
    )
)



# 
server <- function(input, output) {

    output$bootstrap_fc_fc_plot <- renderPlot({
        
        bootstrap_mean_tbl|> filter(gene%in%c(input[["genes_selected"]])) |>
            filter(condition%in%c(input[["condition_x"]])) |> 
            filter(group_size%in%c(input[["sample_size"]])) |>
            mutate(mean_diff_a=a_x-a_y,
                   mean_diff_b=b_x-b_y ) |>
            fc_plot(mean_diff_a,mean_diff_b,gene)
        
})
    
    
    # output$bar_plot <- renderPlot({
    # 
    # bootstrap_mean_tbl |>
    #     filter(condition!="Mock") |> 
    #     mutate(mean_diff_a = a_x - a_y,
    #            mean_diff_b = b_x - b_y) |>
    #     mutate(
    #         quadrant = case_when(
    #             mean_diff_a > 0 & mean_diff_b > 0 ~ "Q1",
    #             mean_diff_a < 0 & mean_diff_b > 0 ~ "Q2",
    #             mean_diff_a < 0 & mean_diff_b < 0 ~ "Q3",
    #             mean_diff_a > 0 & mean_diff_b < 0 ~ "Q4",
    #             TRUE ~ "On Axis"
    #         )
    #     ) |> 
    #     mutate(relation=case_when(
    #         mean_diff_a>mean_diff_b ~"a_larger_b",
    #         mean_diff_a<mean_diff_b ~"b_larger_a",
    #         TRUE ~ "equal")) |> 
    #     mutate(quadrant_ralation=paste(quadrant,relation,sep = "_")) |> 
    #     group_by(gene,condition,group_size,quadrant_ralation#,relation
    #     ) |>
    #     summarise(n=n()) |> 
    #     group_by(gene,condition,group_size) |> 
    #     
    #     mutate(fraction=n/(sum(n))) |>
    #     filter(group_size==700) |> 
    #     ungroup() |> 
    #     #arrange(quadrant, fraction) |> 
    #     mutate(gene=fct_reorder2(gene,quadrant_ralation,fraction)) |> 
    #     ungroup() |> 
    #     ggplot(aes(quadrant_ralation,fraction,fill=quadrant_ralation)) +theme_bw()+ geom_col()+facet_grid(vars(condition),vars(gene))
    # 
    
    
    # output$bar_plot <- renderPlot({
    #     bal_alv_mac |>  dplyr::select(virus,sampletag_Ms4a3) |> ggplot(aes(virus,fill=sampletag_Ms4a3))+geom_bar()+theme_minimal()
    #     
    
         })
}

# Run the application 
shinyApp(ui = ui, server = server)
