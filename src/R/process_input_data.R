#!/usr/bin/env Rscript

'process_input_data

Usage:
  process_input_data.R [-h] (--control <group_name>) <expr_data> <groups> <out_name> <out_dir> <pipeline>

Options:
  -h --help                          Show this help message and exit.
  --control <control>                Control group name. (Required) GSEA comparisons will be reported 
                                     as experimental group vs. control group.

Arguments:
  <expr_data>                        TSV file of gene symbols (rows) by patient IDs (columns) with 
                                     RNA-seq non-normalized counts as values. DESeq2 normalization 
                                     will be performed. Naming requirements: The first column of the 
                                     expr_data TSV file must be "hgnc_symbol". Values requirements:
                                     Values must be RNA-seq non-normalized counts.
  <groups>                           TSV file of patient IDs (column 1) and group labels (column 2).
                                     Naming requirements: The columns of the groups TSV file must be 
                                     "sample_id", "group".
  <out_name>                         Name of GSEA run, to be used for naming output files. 
  <out_dir>                          Path to output directory.
  <pipeline>                         Path to pipeline directory.
' -> doc

suppressPackageStartupMessages(library(docopt))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))

arguments <- docopt(doc)

#####
# temporary: Load test inputs, for development only
# setwd("/projects/karsanlab/jgrants/gsea_automated/src/R")
# source("../../testing/load_test_inputs.R")
#####
pipeline <- arguments$pipeline
source(str_glue("{pipeline}/src/R/functions_gsea_automated.R"))

# Data processing script begins

# Load arguments and check
## Expression data
expr_data <- read_tsv(arguments$expr_data, show_col_types = F)

### check that hgnc_symbol is first column
if (colnames(expr_data[,1]) != "hgnc_symbol") {
  stop("The first column of the expr_data TSV file must be 'hgnc_symbol'.")
}

### check no duplicate sample IDs
if (any(duplicated(colnames(expr_data)))) {
  stop("Sample IDs in expr_data may not contain duplicates.")
}

## Groups data
groups <- read_tsv(arguments$groups, show_col_types = F)

### check column names
if (paste(colnames(groups), collapse = ", ") != "sample_id, group") {
  stop("The columns of the groups TSV file must be 'sample_id', 'group'.")
}

### check "group" is character (categorical)
### TODO: Later add option for continuous phenotype
if (!is.character(groups$group)) {
  stop("The group variable must be categorical (character vector)")
}

### check that group names do not contain spaces
if (any(grepl(" ", unique(groups$group)))) {
  stop("Group names may not contain spaces.")
}

### check that sample names match to expr_data
### non-matching could indicate user error in input data, better for them to fix it themselves
if (!identical(groups$sample_id, colnames(select(expr_data, -hgnc_symbol)), attrib.as.set = F)) {
  stop("Sample ID order does not match between expr_data and groups inputs.")
}

### check that the provided 'control' group matches what is written in the groups data
if (!any(arguments$control %in% groups$group)) {
  stop("The 'control' group argument (-c) does not match any of the group names in the groups.tsv file.")
}

## Use control group to relevel the factor group variable
groups <- mutate(groups, group = relevel(as.factor(group), ref = arguments$control))

## Reorder samples in 'groups' and 'expr_data' so they appear as: {all controls}, {all experimentals}

groups <- arrange(groups, group)

expr_data <- expr_data[,c("hgnc_symbol", groups$sample_id)]

## Name of this run
out_name <- arguments$out_name

## Output directory for this run
out_dir <- arguments$out_dir

## Output file name basename
out_basename <- str_glue("{out_dir}/{out_name}")


# Normalization
## https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

## Formatting for DESeq
print("Normalizing count data using DESeq2")

counts <- column_to_rownames(expr_data, "hgnc_symbol") %>%
  as.matrix()

design_mat <- column_to_rownames(groups, "sample_id")

## Calculate DESeq2 normalized counts
dds <- DESeqDataSetFromMatrix(countData = round(counts, 0), colData = design_mat, design = ~ group)
dds <- estimateSizeFactors(dds)

normalized_counts <- counts(dds, normalized=TRUE)

expr_norm <- as.data.frame(normalized_counts) %>%
  rownames_to_column("hgnc_symbol")

# Generate standard GSEA inputs

## Expression dataset
## format: https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#TXT:_Text_file_format_for_expression_dataset_.28.2A.txt.29
expr_format <- mutate(expr_norm, DESCRIPTION = "na") %>%
  select(NAME = hgnc_symbol, DESCRIPTION, everything())

expr_file <- str_glue("{out_basename}_exprset.txt")

write_delim(expr_format, expr_file, delim = "\t")
print(paste("Saved expression dataset to:", expr_file))

## Phenotype labels
## categorical format: https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#CLS:_Categorical_.28e.g_tumor_vs_normal.29_class_file_format_.28.2A.cls.29
# start here: make a function instead of this ugly mess
write_cls(groups, out_basename)