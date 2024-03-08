# Installation script #2: R package installs

install.packages(c("tidyverse", "here", "janitor", "lemon", "scico", "docopt", "survival"))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")