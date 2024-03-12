# User interface functions

#' Normalize data and run GSEA 
#' 
#' Function to perform count normalization using the DESeq2 package and GSEA analysis using the GSEA command line interface
#' @param expr_data Gene expression dataset (dataframe or TSV). The first column must be named "hgnc_symbol".
#' @param groups Experimental group metadata for each sample (dataframe or TSV). The columns must be named "sample_id", "group". Group names may not contain special characters or spaces.
#' @param control Name of the control group.
#' @param permute_mode Whether to perform GSEA permutations by "phenotype" (random assignment of sample groups) or "gene_set" (randomly select genes).
#' @param out_name Name under which to save output files.
#' @param out_dir Directory in which to save output files.
#' @param GSEAforR Location of the GSEAforR directory.
#' @param load_results Boolean (TRUE/FALSE) for whether to load the GSEA results files
#' @return No return value. Produces a GSEA results directory.
normalize_and_run_GSEA <- function(expr_data, groups, control, permute_mode, out_name, out_dir, GSEAforR) {
  # if inputs are TSVs, leave them as is; if inputs are dataframes, save as TSVs first
  if (is.data.frame(expr_data)) {
    temp_expr_out_name <- str_glue("{out_dir}/{out_name}_temp_input_data.tsv")
    write_tsv(expr_data, temp_expr_out_name)
    expr_data <- temp_expr_out_name
  }
  
  if (is.data.frame(groups)) {
    temp_groups_out_name <- str_glue("{out_dir}/{out_name}_temp_groups_data.tsv")
    write_tsv(groups, temp_groups_out_name)
    groups <- temp_groups_out_name
  }
  
  system(str_glue("sh {GSEAforR}/src/gsea_wrapper_script.sh -e {expr_data} -g {groups} -c '{control}' -r '{permute_mode}' -n '{out_name}' -o {out_dir} -p {GSEAforR}"))

  if (exists("temp_expr_out_name")) {
    file.remove(temp_expr_out_name)
  }
  
  if (exists("temp_groups_out_name")) {
    file.remove(temp_groups_out_name)
  }
  
  # the GSEA command line interface code downloaded from Broad Institute creates an folder named today's date that must be removed
  date_mmmdd <- tolower(format(Sys.Date(), "%b%d"))
  dated_empty_dir <- str_glue("{dirname(out_dir)}/{date_mmmdd}")
  if (dir.exists(dated_empty_dir)) {
    unlink(dated_empty_dir, recursive = T)
  }
}

#' Load GSEA results
#' 
#' Loads the GSEA Hallmark Gene Set Analysis results summary files
#' @param gsea_results_dir GSEA results directory filepath
#' @return Dataframe of GSEA Hallmark Gene Set Analysis results
load_GSEA_results <- function(gsea_results_dir) {
  gsea_results_files <- list.files(gsea_results_dir, pattern = "gsea_report_for_.*\\.tsv", full.names = T)
  
  suppressMessages(result1 <- read_tsv(gsea_results_files[1], show_col_types = F) %>%
    janitor::clean_names())
  suppressMessages(result2 <- read_tsv(gsea_results_files[2], show_col_types = F) %>%
    janitor::clean_names())
  
  result <- bind_rows(result1, result2) %>%
    arrange(desc(nes))
  
  return(result)
}

#' Plot GSEA results
#' 
#' Creates a waterfall plot of the GSEA Hallmark Gene Set Analysis results summary
#' @param result GSEA result, loaded using `load_GSEA_results()`
#' @param fdr_q_cutoff Desired significance cutoff for FDR Q-value
#' @param control_group Name of the control group
#' @param experimental_group Name of the experimental group
#' @param order Order the pathways from largest to smallest NES value (order=1) or smallest to largest NES value (order=-1). Default value: 1
#' @param save_file (Optional) Save the plot to this filepath
#' @return Waterfall plot of GSEA results, with significant results indicated by colour
plot_GSEA_results <- function(result, fdr_q_cutoff, control_group, experimental_group, order=1, save_file=NULL) {
  if(!order %in% c(-1, 1)) {
    stop("order must be 1 or -1")
  }
  
  p <- result %>%
    ggplot(aes(reorder(name, order*nes), nes, fill = fdr_q_val < fdr_q_cutoff)) +
    geom_col() +
    scale_fill_scico_d(palette = "oslo", begin = 0.4, end = 0.9, direction = -1) +
    labs(x = "", y = "NES",
         fill = str_glue("FDR q-value < {fdr_q_cutoff}"),
         title = "GSEA Hallmark Analysis",
         subtitle = str_glue("Comparison: {experimental_group} vs. {control_group}")) +
    coord_flip() +
    theme(axis.text = element_text(colour = "black"),
          legend.key.size = unit(0.5, "cm"),
          legend.margin = margin(0,0,0,-2))
  
  if (!is.null(save_file)) {
    ggsave(save_file, plot = p, width = 10, height = 10)
  }
  
  print(p)
}

