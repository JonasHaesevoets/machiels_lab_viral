bootstrap_mean_tbl
bootstrap_mean_tbl <- read_rds("intermediate_data\\bootstrap_mean_tbl.rds")
library(forcats)

bootstrap_mean_tbl <- bootstrap_mean_tbl #|> head(100)
bootstrap_mean_tbl |>
        filter(condition!="Mock") |> 
        mutate(mean_diff_a = a_x - a_y,
               mean_diff_b = b_x - b_y) |>
        mutate(
                quadrant = case_when(
                        mean_diff_a > 0 & mean_diff_b > 0 ~ "Q1",
                        mean_diff_a < 0 & mean_diff_b > 0 ~ "Q2",
                        mean_diff_a < 0 & mean_diff_b < 0 ~ "Q3",
                        mean_diff_a > 0 & mean_diff_b < 0 ~ "Q4",
                        TRUE ~ "On Axis"
                )
        ) |> 
        mutate(relation=case_when(
                mean_diff_a>mean_diff_b ~"a_larger_b",
                mean_diff_a<mean_diff_b ~"b_larger_a",
                TRUE ~ "equal")) |> 
        mutate(quadrant_ralation=paste(quadrant,relation,sep = "_")) |> 
        group_by(gene,condition,group_size,quadrant_ralation#,relation
                 ) |>
        summarise(n=n()) |> 
        group_by(gene,condition,group_size) |> 
                
        mutate(fraction=n/(sum(n))) |>
        filter(group_size==700) |> 
        ungroup() |> 
        #arrange(quadrant, fraction) |> 
        mutate(gene=fct_reorder2(gene,quadrant_ralation,fraction)) |> 
        ungroup() |> 
        ggplot(aes(quadrant_ralation,fraction,fill=quadrant_ralation)) +theme_bw()+ geom_col()+facet_grid(vars(condition),vars(gene))
