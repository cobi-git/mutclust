# Mutclust

## Overview
This repository contains the codes and data used in our study to identify mutation hotspots of SARS-CoV-2. Our research focuses on identifying hotspot regions, which are areas densely populated with nucleotides where mutations occur frequently, based on the H-score calculated from 224,318 sequences in the GISAID database. These identified hotspots were further applied to the COVID-19 dataset to identify severity-related marker hotspots.

A brief tutorial with example data can be found in the Tutorial.ipynb jupyter script.

## Repository Contents
- example_data/input: Input data for performing Mutclust algorithm.
- src/: Scripts of Mutclust algorithm.

## Data
The dataset used in this study consists of 224,318 SARS-CoV-2 sequences from the GISAID database, spanning the period from January 2020 to November 2022. Sequence alignment was performed using MAFFT v7.490 with hCoV-19/Wuhan/WIV04/2019 (GISAID Access ID: EPI_ISL_402124, GenBank Access: MN908947) as the reference genome.

## Code Availability
The code used for performing Mutclust algorithm and visualize mutation hotspots is available in this repository. To reproduce our results, please follow the instructions provided in the respective directories.

## Getting Started
To get started with the code in this repository, clone the repository to your local machine and install the required dependencies listed in requirements.txt.

