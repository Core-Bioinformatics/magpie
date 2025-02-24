import os
import pandas as pd
import glob

def summarize_folders(parent_folder,
                      verbose=True):
    summary = []

    if os.path.isfile('input/selected.txt'):
        file = open("input/selected.txt", "r")
        samples = [line.rstrip() for line in file]
    else:
        samples = [x.replace('/msi', '').replace("\\","") for x in glob.glob('*/msi',root_dir='input')]
    if os.path.isfile('input/exclude.txt'):
        file = open("input/exclude.txt", "r")
        samples = list(set(samples).difference(set([line.rstrip() for line in file])))
    samples.sort()
    
    # Iterate through subdirectories
    for subdir in samples:
        subdir_path = os.path.join(parent_folder, subdir)

        if os.path.isdir(subdir_path):
            if verbose:
                print(f"Working on sample {subdir}")

            msi_metadata_path = os.path.join(subdir_path, "msi", "MSI_metadata.csv")
            msi_intensities_path = os.path.join(subdir_path, "msi", "MSI_intensities.csv")

            # find if there is an MSI intermediate image
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

            if verbose:
                print("Checking metadata contains spot_id, x and y columns...")
            if not set(['spot_id','x','y']).issubset(msi_metadata.columns):
                missing_cols = set(['spot_id','x','y']) - set(msi_metadata.columns)
                raise ValueError(f"Missing columns in MSI_metadata.csv for sample '{subdir}': {missing_cols}")
            
            if verbose:
                print("Checking intensity table begins with metadata column...")
            if msi_intensities.columns[0] != 'spot_id':
                raise ValueError(f"First column of MSI_intensities.csv is not called 'spot_id' in sample '{subdir}'")

            if verbose:
                print("Checking spot_id columns in metadata and intensity table match exactly")
            if not (msi_metadata['spot_id'].equals(msi_intensities['spot_id'])):
                raise ValueError(f"Column 'spot_id' does not match between MSI_intensities.csv and MSI_metadata.csv in sample '{subdir}'")

            if verbose:
                print("Storing dimensions and examples of data")

            msi_intensities.set_index('spot_id', drop=True, inplace=True)
            msi_intensities_rows, msi_intensities_cols = msi_intensities.shape

            first_peak = msi_intensities.columns[0]
            first_pixel = msi_intensities.index[0]

            summary.append([subdir, msi_he_img, msi_intensities_rows, first_pixel, msi_intensities_cols, first_peak])
    
    # Convert summary to a DataFrame
    summary_df = pd.DataFrame(summary, columns=["Subfolder", "MSI image", "Num pixels", "Name of first pixel", "Num peaks", "Name of first peak"])
    
    return summary_df

if __name__ == "__main__":
    parent_directory = "input/"  # Change this to your actual folder path
    result_df = summarize_folders(parent_directory,verbose=snakemake.params['verbose'])
    
    # Save results to a CSV file
    result_csv_path = os.path.join("output", "summary.csv")
    result_df.to_csv(result_csv_path, index=False)
    if snakemake.params['verbose']:
        print(f"Summary saved to {result_csv_path}")