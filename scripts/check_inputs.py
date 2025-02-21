import os
import snakemake
import pandas as pd
import glob

if os.path.isfile('input/selected.txt'):
    file = open("input/selected.txt", "r")
    samples = [line.rstrip() for line in file]
else:
    samples = glob_wildcards("input/{sample}/msi").sample
if os.path.isfile('input/exclude.txt'):
    file = open("input/exclude.txt", "r")
    samples = list(set(samples).difference(set([line.rstrip() for line in file])))

def summarize_folders(parent_folder):
    summary = []
    
    # Iterate through subdirectories
    for subdir in samples:
        subdir_path = os.path.join(parent_folder, subdir)
        if os.path.isdir(subdir_path):
            
            msi_metadata_path = os.path.join(subdir_path, "msi", "MSI_metadata.csv")
            msi_intensities_path = os.path.join(subdir_path, "msi", "MSI_intensities.csv")
            if glob.glob(subdir_path + '/msi/MSI_HE.*') != []:
                msi_he_img = glob.glob(subdir_path + '/msi/MSI_HE.*')
                msi_he_img = [x for x in msi_he_img if any(ext in x for ext in ['tiff','png','jpg'])]
                if msi_he_img == []:
                    msi_he_img = None
                else :
                    msi_he_img = msi_he_img[0]
            else :
                msi_he_img = None
            
            # Check if a.csv exists and read its dimensions
            if os.path.exists(msi_metadata_path):
                try:
                    msi_metadata = pd.read_csv(msi_metadata_path)
                    msi_metadata_rows, msi_metadata_cols = msi_metadata.shape
                except Exception as e:
                    print(f"Error reading {msi_metadata_path}: {e}")

            if os.path.exists(msi_intensities_path):
                try:
                    msi_intensities = pd.read_csv(msi_intensities_path)
#                    msi_intensities_rows, msi_intensities_cols = msi_intensities.shape
                except Exception as e:
                    print(f"Error reading {msi_metadata_path}: {e}")

            if not set(['spot_id','x','y']).issubset(msi_metadata.columns):
                missing_cols = set(['spot_id','x','y']) - set(msi_metadata.columns)
                raise ValueError(f"Missing columns in MSI_metadata.csv for sample '{subdir}': {missing_cols}")
            
            if msi_intensities.columns[0] != 'spot_id':
                raise ValueError(f"First column of MSI_intensities.csv is not called 'spot_id' in sample '{subdir}'")

            if not (msi_metadata['spot_id'].equals(msi_intensities['spot_id'])):
                raise ValueError(f"Column 'spot_id' does not match between MSI_intensities.csv and MSI_metadata.csv in sample '{subdir}'")

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
    result_df = summarize_folders(parent_directory)
    
    # Save results to a CSV file
    result_csv_path = os.path.join("output", "summary.csv")
    result_df.to_csv(result_csv_path, index=False)
    print(f"Summary saved to {result_csv_path}")