# MAGPIE: Multimodal alignment of genes and peaks for integrative exploration
<p align="center">
<img src="figures/magpie_logo.png" width="200">
</p>

Recent developments in spatially resolved -omics have enabled the joint study of gene expression, metabolite levels and tissue morphology, offering greater insights into biological pathways. 
Integrating these modalities from matched tissue sections to probe spatially-coordinated processes, however, remains challenging. 
Here we introduce _MAGPIE_, a framework for co-registering spatially resolved transcriptomics, metabolomics, and tissue morphology from the same or consecutive sections. 
We show _MAGPIE_’s generalisability and scalability on spatial multi-omics data from multiple tissues, combining Visium with MALDI and DESI mass spectrometry imaging. 
_MAGPIE_ was also applied to new multi-modal datasets generated with a specialised sampling strategy to characterise the metabolic and transcriptomic landscape in an in vivo model of drug-induced pulmonary fibrosis and to link small-molecule co-detection with endogenous lung responses. 
_MAGPIE_ demonstrates the refined resolution and enhanced interpretability that spatial multi-modal analyses provide for studying tissue injury especially in pharmacological contexts, and delivers a modular, accessible workflow for data integration.

_MAGPIE_ has now been published in Nature Communications, please cite the following if you use our tool:

> **Spatially resolved integrative analysis of transcriptomic and metabolomic changes in tissue injury studies.**  
> Eleanor C. Williams, Lovisa Franzén, Martina Olsson Lindvall, Gregory Hamm, Steven Oag,
> Muntasir Mamun Majumder, James Denholm, Azam Hamidinekoo, Javier Escudero Morlanes,
> Marco Vicari, Joakim Lundeberg, Laura Setyo, Trevor M. Godfrey, Livia S. Eberlin,
> Aleksandr Zakirov, Jorrit J. Hornberg, Marianna Stamou, Patrik L. Ståhl, Anna Ollerstam,
> Jennifer Y. Tan, Irina Mohorianu.
> _Nature Communications_, 17, Article 205. (2026)  
> https://doi.org/10.1038/s41467-025-68003-w

## Installation

The MAGPIE pipeline requires a Python installation and the following package dependencies:
* snakemake
* shiny
* matplotlib
* pandas
* numpy
* scikit-image
* pathlib
* scikit-learn
* scipy
* json
* collections
* shutil
* gzip
* h5py
* scanpy

We recommend to create a conda environment with from which the whole pipeline can be run. You can install all required dependencies using the magpie_environment.yml file within the snakemake folder in this repository using the following command:
```
conda env create -f magpie_environment.yml
```

The pipeline has been previously tested on the following systems:
* macOS: Sequoia (15.3.2)
* Windows: 11 (22H2)

Installation should take up to ~10 minutes on a normal desktop computer.

## Input structure

The MAGPIE pipeline automatically detects the files in your input folder and makes decisions accordingly so you must ensure your files follow the following structure:

    [sample name]
    ├── visium                               # Spaceranger outputs
    │   ├── filtered_feature_bc_matrix.h5
    │   ├── spatial
    │   │   ├── aligned_fiducials.jpg
    │   │   ├── detected_tissue_image.jpg
    │   │   ├── scalefactors_json.json
    │   │   ├── tissue_hires_image.png
    │   │   ├── tissue_lores_image.png
    │   │   ├── tissue_positions_list.csv
    ├── msi                    
    │   ├── MSI_intensities.csv              # Table of intensities with MSI peaks on columns and pixels on rows
    │   ├── MSI_metadata.csv                 # Table of metadata about MSI pixels, including x and y coordinate columns
    │   │── MSI_HE.[jpg,png,tiff]            # (OPTIONAL) intermediate MSI image to assist with coregistration
    ├── landmarks_MSI2HE.csv                 # (OPTIONAL) Table of identified landmarks between MSI image and MSI H&E image (added by shiny app or identified externally)
    ├── landmarks_MSI2HE.csv                 # (OPTIONAL) Table of identified landmarks between MSI H&E and Visium H&E image (added by shiny app or identified externally)
    └── landmarks_noHE.csv                   # (OPTIONAL) Table of identified landmarks between MSI image and Visium H&E (added by shiny app or identified externally). 
                                               landmarks_noHE.csv or landmarks_MSI2HE.csv and landmarks_MSI2HE.csv are required for coregistration.
    


## Running the shiny app

To run the pipeline, you need to be in the folder with all files in the _snakemake_ folder in this repository as well as an _input_ folder as described in the previous section.

To start the shiny app for manual landmark selection, run ``` shiny run magpie_shiny_app.py ```

For each sample you will be prompted to select some manual landmarks then download. At the point you download them they will be saved into the file structure described above. If you would prefer to use your own landmarks please save them into that structure instead and you can skip the shiny app step.

## Running the snakemake pipeline

Once landmarks have been selected for each sample, you can switch to the snakemake pipeline to perform the coregistration. Again you must be in the folder with all files in the _snakemake_ folder in this repository as well as an _input_ folder as described in the previous section with your newly selected landmarks. You can then run the pipeline using ``` snakemake --cores [n] ``` where _n_ is the number of cores you would like to use. You can explicitly state which samples you would like to use by listing them in a *selected.txt* file within the *input* folder and equivalently specify some files you would like to exclude using a *exclude.txt* file.

## Tutorial

We provide extensive [documentation describing the pipeline](https://core-bioinformatics.github.io/magpie/) and a tutorial with example data described [here](https://core-bioinformatics.github.io/magpie/tutorial/SMA_tutorial.html). The tutorial should take around 5-10 minutes to run. 
