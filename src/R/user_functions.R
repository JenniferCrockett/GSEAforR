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
normalize_and_run_GSEA <- function(expr_data, groups, control, permute_mode, out_name, out_dir, GSEAforR) {
  # if inputs are dataframes, leave them as is; if inputs are TSVs, load them
  if (!is.data.frame(expr_data) & file.exists(expr_data)) {
    expr_data <- read_tsv(expr_data, show_col_types = FALSE)
  }
  
  if (!is.data.frame(groups) & file.exists(groups)) {
    groups <- read_tsv(groups, show_col_types = FALSE)
  }
  
  system(str_glue("sh {GSEAforR}/src/gsea_automated.sh -e {expr_data} -g {groups} -c '{control}' -r '{permute_mode}' -n '{out_name}' -o {out_dir} -p {GSEAforR}"))
}

#' Load GSEA results
#' 
#' @param name description
load_GSEA_results <- function() {
  
}

#' Plot GSEA results
#' 
#' @param
plot_GSEA_results <- function() {
  
}