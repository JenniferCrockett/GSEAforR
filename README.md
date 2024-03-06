# GSEA functions for R

## About

### What is GSEA?

[Gene Set Enrichment Analysis (GSEA, Broad Institute)](https://www.gsea-msigdb.org/gsea/index.jsp) is a widely-used software by biology researchers. It allows researchers to input gene expression data, which is a measure of how 'active' or 'inactive' genes are, for cells that have undergone different drug treatments, genetic treatments, or other experimental conditions (let's call this 'conditions' from now on). GSEA analyzes which pathways, represented by 'gene sets', are activated or inactivated by a particular condition. Determining which pathways respond to a condition can help researchers to understand the biological mechanism that is happening when cells experience that condition.  

Examples of why this information is useful to researchers:  

* Researchers want to know how newly discovered pharmaceutical agents affect cells - which cellular pathways are turned on or off when the drug works on cells?
* Researchers want to know why a particular gene mutation causes a genetic disease or cancer - which cellular pathways are turned on or off when the gene mutation is introduced into cells in the lab?

### Why use _GSEA functions for R_

Normally researchers run GSEA as follows:  

1. Load gene expression counts data _(typically performed in R)_.
2. Perform data normalization _(typically performed in R)_.
3. Format the data according to the [GSEA User Guide](https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html) and [Data Formats](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats) specifications _(typically performed in R)_.
4. Save the formatted data _(typically performed in R)_.
5. Load the data into the GSEA app.
6. Run GSEA to produce GSEA results files.
7. Analyze the results files and plot the results _(typically performed in R)_.  

Notice how the above workflow typically is performed in R, **except for the GSEA analysis steps**. _GSEA functions for R_ fills the gap to create a seamless workflow in R, by providing functions to execute GSEA in an R environment.  

Here is the GSEA workflow using _GSEA functions for R_, **all in R:** 

1. Load gene expression counts data.
2. Format the data, in a much simpler format than what is required by the GSEA Data Formats guide.
3. Run `normalize_and_run_GSEA()` to automatically perform: data normalization, data formatting for GSEA input, run the GSEA source code (open source code provided by the GSEA project), produce GSEA results files.
4. Analyze the results files and plot the results using `load_GSEA_results()` and `plot_GSEA_results()`.  

## Installation

1. Clone this repository:

```
git clone https://github.com/JenniferCrockett/GSEAforR.git
```

2. In the command line terminal, run the `setup/installation.sh` script to set up a GSEAforR conda environment and install the necessary software.
3. In the `src/gsea_wrapper_script.sh` script's _Magic Strings_ section: update the `CONDA_ENV=` variable to the path to your GSEAforR conda environment, and update the `PIPELINE=` variable to the path to your GSEAforR directory.  

To see the path to the conda environment, run the following in the terminal:
```
conda info --env
```

4. In RStudio, run the `setup/installation.R` script to install the necessary R packages to run the GSEAforR functions in your R environment.

## Usage

## Sources

1. GSEA v4.3.3 for the command line: Downloaded from https://www.gsea-msigdb.org/gsea/downloads.jsp
2. GSEA run command: Accessed using the "command" feature of the GSEA app (https://www.gsea-msigdb.org/gsea/downloads.jsp)
3. GSEA Hallmark Gene Set: Downloaded from https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#H (Gene Symbols)
