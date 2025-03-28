# MAGPIE: Multimodal alignment of genes and peaks for integrative exploration
<p align="center">
<img src="figures/magpie_logo.png" width="200">
</p>

Recent developments in spatially resolved -omics have enabled studies linking gene expression and metabolite levels to tissue morphology, offering new insights into biological pathways. By capturing multiple modalities on matched tissue sections, one can better probe how different biological entities interact in a spatially coordinated manner. However, such cross-modality integration presents experimental and computational challenges. 
To align multimodal datasets into a shared coordinate system and facilitate enhanced integration and analysis, we propose *MAGPIE* (Multi-modal Alignment of Genes and Peaks for Integrated Exploration), a framework for co-registering spatially resolved transcriptomics, metabolomics, and tissue morphology from the same or consecutive sections. 
We illustrate the generalisability and scalability of *MAGPIE* on spatial multi-omics data from multiple tissues, combining Visium with both MALDI and DESI mass spectrometry imaging. *MAGPIE* was also applied to newly generated multimodal datasets created using specialised experimental sampling strategy to characterise the metabolic and transcriptomic landscape in an in vivo model of drug-induced pulmonary fibrosis, to showcase the linking of small-molecule co-detection with endogenous responses in lung tissue.
*MAGPIE* highlights the refined resolution and increased interpretability of spatial multimodal analyses in studying tissue injury, particularly in pharmacological contexts, and offers a modular, accessible computational workflow for data integration.

Preprint: https://www.biorxiv.org/content/10.1101/2025.02.26.640381v1

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
