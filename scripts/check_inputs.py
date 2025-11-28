import os
import pandas as pd
import glob

def summarize_folders(parent_folder,
                      verbose=True):
    """
    Summarise all MSI datasets contained in subfolders of a parent directory.

    This function:
    - Identifies sample folders automatically, or reads them from `input/selected.txt`
    - Excludes samples listed in `input/exclude.txt`
    - For each sample, checks the presence and integrity of:
        * MSI metadata table (must contain spot_id, x, y)
        * MSI intensity table (must begin with spot_id)
        * Optional intermediate MSI H&E image (MSI_HE.[tiff/png/jpg])
    - Confirms metadata and intensity tables align on spot_id
    - Extracts basic dataset dimensions and example entries
    - Returns a summary DataFrame listing per-sample checks and properties

    Parameters
    ----------
    parent_folder : str
        Path to the folder that contains all sample subdirectories.
    verbose : bool, optional
        Whether to print progress updates during processing.

    Returns
    -------
    pd.DataFrame
        Summary table with columns:
        ["Subfolder", "MSI image", "Num pixels", "Name of first pixel",
         "Num peaks", "Name of first peak"]
    """

    summary = []

    # --------------------------------------------------------------
    # Determine which samples to process
    # --------------------------------------------------------------

    # Priority 1: explicitly selected samples
    if os.path.isfile('input/selected.txt'):
        if verbose:
            print("Using sample list from input/selected.txt")
        with open("input/selected.txt", "r") as file:
            samples = [line.rstrip() for line in file]
    else:
        # Priority 2: infer sample folders by scanning for */msi directories
        samples = [
            x.replace('/msi', '').replace("\\", "")
            for x in glob.glob('*/msi', root_dir='input')
        ]

    # Remove excluded samples if exclude.txt exists
    if os.path.isfile('input/exclude.txt'):
        if verbose:
            print("Excluding samples from input/exclude.txt")
        with open("input/exclude.txt", "r") as file:
            excluded = set(line.rstrip() for line in file)
        samples = list(set(samples).difference(excluded))

    samples.sort()

    # --------------------------------------------------------------
    # Iterate over samples and perform validity checks
    # --------------------------------------------------------------
    for subdir in samples:
        subdir_path = os.path.join(parent_folder, subdir)

        if os.path.isdir(subdir_path):

            if verbose:
                print(f"\n===== Working on sample: {subdir} =====")

            # Paths to core MSI files
            msi_metadata_path = os.path.join(subdir_path, "msi", "MSI_metadata.csv")
            msi_intensities_path = os.path.join(subdir_path, "msi", "MSI_intensities.csv")

            # ----------------------------------------------------------
            # Detect presence of an intermediate MSI H&E image (optional)
            # ----------------------------------------------------------
            if verbose:
                print("Checking for intermediate MSI image...")
            if glob.glob(subdir_path + '/msi/MSI_HE.*') != []:
                msi_he_img = glob.glob(subdir_path + '/msi/MSI_HE.*')
                msi_he_img = [x for x in msi_he_img if any(ext in x for ext in ['tiff','png','jpg'])]
                if msi_he_img == []:
                    msi_he_img = None
                else :
                    msi_he_img = msi_he_img[0]
                    if verbose:
                        print(f"MSI image {msi_he_img} found")
            else :
                msi_he_img = None
            
            # ----------------------------------------------------------
            # Read metadata and intensity tables
            # ----------------------------------------------------------
            if verbose:
                print("Checking metadata and intensity tables exist and can be read...")
            # Check if metadata file exists
            try:
                msi_metadata = pd.read_csv(msi_metadata_path)
            except Exception as e:
                print(f"Error reading {msi_metadata_path}: {e}")

        # Check if intensity file exists
            try:
                msi_intensities = pd.read_csv(msi_intensities_path)
            except Exception as e:
                print(f"Error reading {msi_metadata_path}: {e}")

            # ----------------------------------------------------------
            # Column validation for metadata
            # ----------------------------------------------------------
            if verbose:
                print("Checking metadata contains spot_id, x and y columns...")
            if not set(['spot_id','x','y']).issubset(msi_metadata.columns):
                missing_cols = set(['spot_id','x','y']) - set(msi_metadata.columns)
                raise ValueError(f"Missing columns in MSI_metadata.csv for sample '{subdir}': {missing_cols}")
            

            # ----------------------------------------------------------
            # Validate intensity table structure
            # ----------------------------------------------------------
            if verbose:
                print("Checking intensity table format...")

            if msi_intensities.columns[0] != 'spot_id':
                raise ValueError(
                    f"First column of MSI_intensities.csv must be 'spot_id' (sample '{subdir}')"
                )

            # ----------------------------------------------------------
            # Validate sample alignment
            # ----------------------------------------------------------
            if verbose:
                print("Checking spot_id alignment between tables...")

            if not msi_metadata['spot_id'].equals(msi_intensities['spot_id']):
                raise ValueError(
                    f"spot_id mismatch between MSI_metadata.csv and MSI_intensities.csv (sample '{subdir}')"
                )

            # ----------------------------------------------------------
            # Summarize dataset properties
            # ----------------------------------------------------------
            if verbose:
                print("Extracting dataset dimensions and examples...")

            msi_intensities.set_index('spot_id', inplace=True)
            num_pixels, num_peaks = msi_intensities.shape
            first_pixel = msi_intensities.index[0]
            first_peak = msi_intensities.columns[0]

            summary.append([
                subdir,
                msi_he_img,
                num_pixels,
                first_pixel,
                num_peaks,
                first_peak
            ])

    # Convert summary to a DataFrame
    summary_df = pd.DataFrame(
        summary,
        columns=[
            "Subfolder",
            "MSI image",
            "Num pixels",
            "Name of first pixel",
            "Num peaks",
            "Name of first peak"
        ]
    )

    return summary_df


if __name__ == "__main__":
    parent_directory = "input/"
    result_df = summarize_folders(parent_directory, verbose=snakemake.params['verbose'])

    # Save results
    result_csv_path = os.path.join("output", "summary.csv")
    result_df.to_csv(result_csv_path, index=False)

    if snakemake.params['verbose']:
        print(f"\nSummary saved to {result_csv_path}")
