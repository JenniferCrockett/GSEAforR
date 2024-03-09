# User interface functions

#' Normalize data and run GSEA 
#' 
#' Function to perform count normalization using the DESeq2 package and GSEA analysis using the GSEA command line interface
#' @param expr_data Gene expression dataset. The first column must be named "hgnc_symbol"
#' @param groups Experimental group metadata for each sample. The columns must be named "sample_id", "group". The samples must be in the same order as in expr_data.
#' @param control The name of the control group.
#' @param permute_mode Whether to perform GSEA permutations by "phenotype" (random assignment of sample groups) or "gene_set" (randomly select genes).
#' @param out_name The name under which to save output files.
#' @param out_dir The directory in which to save output files.
#' @param GSEAforR The location of the GSEAforR directory.
normalize_and_run_GSEA <- function(expr_data, groups, control, permute_mode, out_name, out_dir, GSEAforR) {
  system(str_glue("{GSEAforR}/src/gsea_automated.sh -e {expr_data} -g {groups} -c '{control}' -r '{permute_mode}' -n '{out_name}' -o {out_dir} -p {GSEAforR}"))
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