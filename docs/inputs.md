(input-structure)=
# MAGPIE directory structure

The directory from which you run the pipeline should have the following general structure:

    ├── Snakefile
    ├── magpie_shiny_app.py
    ├── figures  
    │   ├── magpie_logo.png
    ├── scripts
    │   ├── alter_data.py
    │   ├── create_mock_spaceranger.py
    │   ├── create_perbarcode_matrix.py
    ├── input
    │   ├── (optional) exclude.txt
    │   ├── (optional) selected.txt
    │   ├── ... 

The exclude.txt and selected.txt files are optional and allow the user to list a number of samples to include (rather than using all sub-folders in the input directory by default) or to exclude from analysis. 

The Snakefile, magpie_shiny_app.py and scripts and figures folders should all be copied from the GitHub repository.

(inputstructure)=
## Input folder structure

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
    │   │── MSI_HE.[jpg,png,tiff]                       # (OPTIONAL) intermediate MSI image to assist with coregistration
    ├── landmarks_MSI2HE.csv                 # (OPTIONAL) Table of identified landmarks between MSI image and MSI H&E image (added by shiny app or identified externally)
    ├── landmarks_MSI2HE.csv                 # (OPTIONAL) Table of identified landmarks between MSI H&E and Visium H&E image (added by shiny app or identified externally)
    └── landmarks_noHE.csv                   # (OPTIONAL) Table of identified landmarks between MSI image and Visium H&E (added by shiny app or identified externally). 
                                               Either landmarks_noHE.csv or both landmarks_MSI2HE.csv and landmarks_MSI2HE.csv are required for coregistration.

## MSI data

The MSI data should be split into an intensity table and a metadata table. The intensity table should look like this, with peaks on the columns and pixels on the rows. 

|spot_id|mz-101.2345   | mz-102.85754 | mz-303.35855   | mz-344.48575  | mz-321.38583  | mz-112.28485 |
|:---:|:---:|:---:|:---:|:---:|:---:|:---: |
|pixel_0| 5774.812  | 675.361  | 23.555  | 8444.958  | 777.234  | 20.332  |
|pixel_1| 8794.013  | 444.523  | 81.294  | 6775.393  | 899.284  | 10.275  |
|pixel_2| 6777.358  | 857.585  | 14.326  | 9468.367  | 747.385  | 24.521  |
|| ...  | ...  | ...  | ...  | ... | ... |

The metadata table must contains columns spot_id, x and y (and can additionally include others) and the spot_id columns must exactly match the spot_id column in the intensity table.

| spot_id           | x   | y   | Sample   | Treatment | Fibrotic |
| :---:        | :---:       | :---:       | :---:       | :---:   | :---: |
| pixel_0 | 1        | 1        | sample_1        | treatment_1 | normal |
| pixel_1 | 2        | 1        | sample_1        | treatment_1 | normal |
| pixel_2 | 2        | 2        | sample_1       | treatment_1 | fibrotic |
|| ...  | ...  | ...  | 

The first stage of the snakemake pipeline will check all these criteria and save a summary table called summary.txt. 
