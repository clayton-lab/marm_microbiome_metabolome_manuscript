# Repository for the marm_microbiome_metabolome_manuscript

## Overview
This repository hosts the software used to analyze and generate data for the marm_microbiome_metabolome manuscript. The repo consists of 3 directories that are described below:
- **src** — The code used to analyze 16S and metabolomics data. The majority of the code is contained within a Jupyter notebook (marm-microbiome-metabolome-manuscript.ipynb), with some helper functions in the file called q2Handler.py. An additional Cytoscape file (marm-micro-metab-associations.cys) contains the portion of analysis that took place within Cytoscape.
- **env** — The Conda environment that was used to run the code in src. Contains a .yaml file that can be used to recreate the environment. While Dokdo was installed via pip, the process for doing so was slightly unconventional. Further details can be found in Dokdo's [documentation](https://dokdo.readthedocs.io/en/latest/readme.html).
- **results** — Processed 16S and metabolomics data, as well as the output files from several statistical analysis tools. The files in this folder should be sufficient to recreate any figures present in the published manuscript, with the caveat that additional steps in Cytoscape may be necessary for some figures.
