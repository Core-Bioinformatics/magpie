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
    │   │── MSI_HE.jpg                       # (OPTIONAL) intermediate MSI image to assist with coregistration
    ├── landmarks_MSI2HE.csv                 # (OPTIONAL) Table of identified landmarks between MSI image and MSI H&E image (added by shiny app or identified externally)
    ├── landmarks_MSI2HE.csv                 # (OPTIONAL) Table of identified landmarks between MSI H&E and Visium H&E image (added by shiny app or identified externally)
    └── landmarks_noHE.csv                   # (OPTIONAL) Table of identified landmarks between MSI image and Visium H&E (added by shiny app or identified externally). 
                                               landmarks_noHE.csv or landmarks_MSI2HE.csv and landmarks_MSI2HE.csv are required for coregistration.
