#!/bin/sh

# Installation script #1: command line installs

## Create conda environment
conda create --name "GSEAforR"
source activate GSEAforR

## Install required packages in conda environment
conda install conda-forge::openjdk=11.0.1

conda deactivate