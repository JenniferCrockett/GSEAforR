# GSEAforR

## About

### What is GSEA?

[Gene Set Enrichment Analysis (GSEA, Broad Institute)](https://www.gsea-msigdb.org/gsea/index.jsp) is a widely-used software by biology researchers. It allows researchers to input gene expression data, which is a measure of how 'active' or 'inactive' genes are, for cells that have undergone different drug treatments, genetic treatments, or other experimental conditions. GSEA analyzes which pathways, represented by 'gene sets', are activated or inactivated by a particular condition. Determining which pathways respond to a condition can help researchers to understand the biological mechanism that is happening when cells experience that condition.  

Examples of why this information is useful to researchers:  

* Researchers want to know how new potential pharmaceutical agents affect cells - which cellular pathways are turned on or off when the drug works on cells?
* Researchers want to know why a particular gene mutation causes a genetic disease or cancer - which cellular pathways are turned on or off when the gene mutation is introduced into cells in the lab?

### Why use _GSEAforR_?

Normally researchers run GSEA as follows:  

1. Load gene expression counts data _(typically performed in R)_.
2. Perform data normalization _(typically performed in R)_.
3. Format the data according to the [GSEA User Guide](https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html) and [Data Formats](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats) specifications _(typically performed in R)_.
4. Save the formatted data _(typically performed in R)_.
5. Load the data into the GSEA app.
6. Run GSEA to produce GSEA results files.
7. Analyze the results files and plot the results _(typically performed in R)_.  

Notice how the above workflow typically is performed in R, **except for the GSEA analysis steps**. _GSEAforR_ fills the gap to create a seamless workflow in R, by providing functions to execute GSEA in an R environment.  

Here is the GSEA workflow using _GSEAforR_, **all in R:** 

1. Load gene expression counts data.
2. Format the data, in a much simpler format than what is required by the GSEA Data Formats guide.
3. Run `normalize_and_run_GSEA()` to automatically perform: data normalization, data formatting for GSEA input, run the GSEA source code (open source code provided by the GSEA project), produce GSEA results files.   
4. Load the GSEA results summary files using `load_GSEA_results()`.
5. Plot the results using `plot_GSEA_results()`.  


## Installation

1. Clone this repository:

```
git clone https://github.com/JenniferCrockett/GSEAforR.git
```

2. In the command line terminal, run the `setup/installation.sh` script to set up a GSEAforR conda environment and install the necessary software.

```
# from the GSEAforR directory, run the following
sh ./setup/installation.sh
```

3. In RStudio, run the `setup/installation.R` script to install the required R packages.

4. In the `src/gsea_wrapper_script.sh` script's _Magic Strings_ section, update the following:  

* Update the `CONDA_ENV=` variable to the path to your GSEAforR conda environment
* Update the `RPATH=` variable to your R path

**Tips:**  

To see the path to the conda environment, run the following in the terminal:
```
conda info --env
```

To see your R path, run the following in the terminal:  
```
which R
```

## Usage

**An example analysis is included in `GSEAforR/testing/Analysis_using_GSEAforR.Rmd`.**

In an Rmarkdown document, load the GSEAforR functions:  

```
devtools::load_all("/absolute/path/to/GSEAforR/src/R")
```

Formatting input data:  

* Expression data: Gene symbols (rows) by patient IDs (columns) with RNA-seq non-normalized counts as values. The first column name must be *hgnc_symbol*.
* Patient groups data: Patient IDs (column 1) and group labels (column 2). The column names must be *sample_id*, *group*. The patient IDs included in the *sample_id* column must be the same as the patient IDs included in the expressio data column names. 
* Accepted formats: data.frame or filepath to a saved TSV file

Normalize data and run GSEA analysis:  

```
normalize_and_run_GSEA(expr_data = expression_data, 
                       groups = groups_data, 
                       control = "my_control_group_name", 
                       permute_mode = "gene_set", 
                       out_name = "my_experiment_name", 
                       out_dir = "/absolute/path/to/my/output/directory", 
                       GSEAforR = "/absolute/path/to/GSEAforR/")
```

Load GSEA results:  

```
result <- load_GSEA_results(gsea_results_dir = "./path/to/GSEA_results_directory")
```

Plot GSEA results:

```
plot_GSEA_results(result = result, 
                  fdr_q_cutoff = 0.05, 
                  control_group = "My Control Group Name", 
                  experimental_group = "My Experimental Group Name", 
                  order = -1, 
                  save_file = "output_filename.png")
```

For detailed information about these functions, see their help documentation by running:  

```
?normalize_and_run_GSEA
?load_GSEA_results
?plot_GSEA_results
```

## Updates

Gene set updates are released by MSigDB annually or biannually. The current gene set provided in this repository is the Hallmark Gene Set _h.all.v2023.2_. If you wish to manually update to a new version:  

1. Create a new versioned directory in the `GSEAforR/resources` directory.
2. To the versioned directory, download the new gene set (gmt) file by clicking on the "Gene Symbols" link under the H: hallmark gene set section on the following web page: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#H
3. Create a symlink for this file in the `GSEAforR/resources/MSigDB_latest` directory, and remove the symlink to the old file.

## Sources

1. GSEA v4.3.3 for the command line: Downloaded from https://www.gsea-msigdb.org/gsea/downloads.jsp
2. GSEA run command: Accessed using the "command" feature of the GSEA app (https://www.gsea-msigdb.org/gsea/downloads.jsp)
3. GSEA Hallmark Gene Set: Downloaded from https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#H (Gene Symbols)
4. Test data: Downloaded from the GEO data repository at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72737
