# Robust enhancer-gene association prediction using single cell transcriptomes and epigenomes


This is the repository that hosts customized scripts for the analysis of enhancer-gene associations in neurons from mouse primary motor cortex.

Reference:
- [Xie, Armand et al. 2021; Robust enhancer-gene regulation identified by single-cell transcriptomes and epigenomes](https://www.biorxiv.org/content/10.1101/2021.10.25.465795v1)

Correspondence: [Eran A. Mukamel](mailto:emukamel@ucsd.edu) and [Fangming Xie](mailto:f7xie@ucsd.edu)

# Getting started
**System requirements**

This package is tested using a Ubuntu 18.04.6 LTS (Bionic Beaver) server. However, we expect it can be operated under a wide range of systems.
We recommend users to use a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html) to install dependencies.

**Installation**

After pre-installing Anaconda, run the following command to clone this repo and install dependencies.
```bash

# clone this repo
git clone https://github.com/FangmingXie/scf_enhancer_paper.git

# create an conda env and install dependancies.
cd SingleCellFusion_EnhancerPaper
conda env create -f env.yml
```
The installation of the conda environment takes less than 20 minutes. After installation, activate the environment using
```bash
conda activate conda_enhancer
```

**Demo**

`demo/data`
`demo/script`
The whole demo takes about xx minutes to run through.


**Instructions for use**

- how to run your data
