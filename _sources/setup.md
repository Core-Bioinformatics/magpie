# Setup Instructions

To run this pipeline, you need to have `Snakemake` installed along with several bioinformatics tools. We recommend using `conda` for package management.

## 1. Install Conda

If you haven't already, install Miniconda using instructions at https://docs.anaconda.com/miniconda/install/, for example using the below.

For Mac:
```bash
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
bash ~/Miniconda3-latest-MacOSX-arm64.sh
```

For Linux:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash ~/Miniconda3-latest-Linux-x86_64.sh
```

For Windows:
```bash
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe --output .\Downloads\Miniconda3-latest-Windows-x86_64.exe
```
Then launch the installer from Downloads and follow the steps


## 2. Create Conda Environment 

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

We recommend to create a conda environment with from which the whole pipeline can be run. You can install all required dependencies using the magpie_environment.yml file within the snakemake folder in the GitHub repository using the following command:

```bash
conda env create -f magpie_environment.yml
```

The pipeline has been primary developed and tested on Mac so there may be some difficulties with installation on Windows and Linux and packages may need to be installed manually.

## 3. Activate Environment

The MAGPIE Conda environment can then be activated using

```bash
conda activate magpie
```

## 4. Copy the MAGPIE pipeline locally

In order to run the MAGPIE pipeline, you need to have several files from the GitHub repository https://github.com/Core-Bioinformatics/magpie locally. Specifically you need the structure:

    ├── Snakefile
    ├── magpie_shiny_app.py
    ├── figures  
    │   ├── magpie_logo.png
    ├── scripts
    │   ├── alter_data.py
    │   ├── create_mock_spaceranger.py
    │   └── create_perbarcode_matrix.py

You will then create a folder in the same directory called 'input' where all your inputs will be placed, following the structure specified in [the Input folder structure section](inputstructure). 

In order to get these files locally, you can clone the reposity, either using git e.g. `git clone https://github.com/Core-Bioinformatics/magpie.git`, using the GitHub Desktop app by going to `File > Clone Repository...` or by downloading the whole repository as a zip file by going to the repository at https://github.com/Core-Bioinformatics/magpie then going to `Code > Download ZIP` then extracting the resulting folder locally.
