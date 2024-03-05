#!/bin/sh
# Script to run GSEA from the command line

set -e
set -o pipefail

usage()
{
    echo "
gsea_automated.sh : Processes input data and runs GSEA from the command line.

Usage:
-h Show this help message and exit.
-e expr_data.tsv (Required)
 	TSV file of gene symbols (rows) by patient IDs (columns) with 
	RNA-seq non-normalized counts as values. DESeq2 normalization 
	will be performed. Naming requirements: The first column of the 
	expr_data TSV file must be \"hgnc_symbol\". Values requirements:
	Values must be RNA-seq non-normalized counts.
-g groups.tsv (Required)
	TSV file of patient IDs (column 1) and group labels (column 2).
 	Naming requirements: The columns of the groups TSV file must be 
	sample_id, group.
-c control (Required)
    Control group name. GSEA comparisons will be reported as
    experimental group vs. control group.
-n out_name (Required)
	Name of GSEA run, to be used for naming output files.
-r permute mode (Required)
    Set GSEA permute mode to 'phenotype' (permutation on phenotype)
    or 'gene_set' (permutation on gene set).
-x gmx selection (Required)
    Select GSEA gmx gene set collection.
    Valid options are: 'h' (Hallmark gene sets), 'c2' (C2 curated 
    gene sets), 'c3' (C3 regulatory target gene sets: microRNA and 
    TF targets), 'c5' (C5 ontology gene sets: GO terms), or a path 
    to a custom .gmt file. 
-v gmx version (Optional)
    Name of desired gmx version directory. 
    Default: latest (/projects/karsanlab/RESOURCES/msigdb/human/latest)
    Example: v2023.2 (for: /projects/karsanlab/RESOURCES/msigdb/human/v2023.2)
-o out_dir (Required)
    Path to output directory.
-p pipeline (Required)
    Path to pipeline directory.
"
}


# Step 1: Processing arguments
echo "Running GSEA automated pipeline"

while getopts ":e:g:c:n:o:p:hr:" opt; do
case $opt in
    h)
        usage
        exit 0
        ;;
    e)
        expr_data="$OPTARG"
        ;;
    g)
        groups="$OPTARG"
        ;;
    c)
        control="$OPTARG"
        ;;
    n)
        out_name="$OPTARG"
        ;;
    o)
        out_dir="$OPTARG"
        ;;
    p)
        PIPELINE="$OPTARG"
        ;;
    r)
        permute="$OPTARG"
        ;;
    \?)
        echo "invalid option"
        exit 1
esac
done
shift $((OPTIND-1)) # this allows you to parse position arguments at the end if you have any

## check validity of input files
if [ -z "$expr_data" ] || [ ! -f "$expr_data" ]
then
    usage
    echo "Please provide an existing expr_data file."
    exit 1
fi

if [ -z "$groups" ] || [ ! -f "$groups" ]
then
    usage
    echo "Please provide an existing groups file."
    exit 1
fi


## check validity of permute mode selection
if [ "$permute" != "phenotype" ] && [ "$permute" != "gene_set" ]
then
    echo "Valid options for -r (permute) are 'phenotype' or 'gene_set'"
    exit 1
fi


