import os
import json
from collections import namedtuple
from scipy.spatial.distance import pdist
from scipy.sparse import csr_matrix, coo_array
import pandas as pd
import shutil
from scipy.io import mmwrite
import gzip
import h5py


def create_mock_spaceranger(
    sample_name='sample1',
    output_folder_name="mock_spaceranger",
    msi_image_path="msi_he.png",
    msi_coord_fname="msi_coords.csv",
    msi_spot_prefix="MSI", 
    visium_sf_json_path=None,
    msi_peak_data_path=None,
    verbose=True
):
    """
    Creates mock Space Ranger formatted output data for a sample.

    Args:
        sample (str): Name of sample which will be used to identify the relevant existing input/output files. Defaults to 'sample1'.
        output_folder_name (str, optional): Name of the output folder. Defaults to "mock_spaceranger".
        msi_image_path (str, optional): Path to the MSI image file. Defaults to "msi_he.png".
        msi_coord_fname (str, optional): Name of the MSI coordinate file. Defaults to "msi_coords.csv".
        msi_spot_prefix (str, optional): Prefix for spot IDs in the MSI data. Defaults to "MSI".
        visium_sf_json_path (str): Path to the Visium scalefactor JSON file.
        msi_peak_data_path (str, optional): Path to the MSI peak data file.
        verbose (bool, optional): Print informational messages during execution. Defaults to True.

    Returns:
        None
    """

    # Create output folders
    output_path = os.path.join(output_folder_name)
    os.makedirs(output_path, exist_ok=True)
    spatial_path = os.path.join(output_path, "spatial")
    os.makedirs(spatial_path, exist_ok=True)
    filtered_path = os.path.join(output_path, "filtered_feature_bc_matrix")
    os.makedirs(filtered_path, exist_ok=True)

    # Load matching Visium scalefactor JSON file
    with open(visium_sf_json_path, "r") as f:
        st_json = json.load(f)
    scale_factor = st_json['tissue_hires_scalef']
    msi_json = st_json

    # Read MSI coordinates
    if verbose:
        print("Reading MSI coordinate file...")
        
    msi_coords = pd.read_csv(msi_coord_fname)
    msi_coords['x'] = msi_coords['x']/scale_factor
    msi_coords['y'] = msi_coords['y']/scale_factor
    msi_coords['x_2'] = msi_coords['x']
    msi_coords['y_2'] = msi_coords['y']
    msi_coords.set_index('spot_id', drop=False, inplace=True)

    # Add "in_tissue" column
    msi_coords["in_tissue"] = 1
    msi_tissue_pos = msi_coords.rename(
        columns={
            "spot_id":"barcode",
            'x_2': "array_col",
            'y_2': "array_row",
            'x': "pxl_col_in_fullres",
            'y': "pxl_row_in_fullres",
        }
    )
    
    # Rename "barcodes" so that they are characters
    msi_tissue_pos["barcode"] = msi_spot_prefix + "_" + sample_name + "_" + msi_tissue_pos["barcode"].astype(str)

    msi_tissue_pos = msi_tissue_pos[['barcode','in_tissue','array_row','array_col','pxl_row_in_fullres','pxl_col_in_fullres']]

    # Write tissue_positions.csv
    msi_tissue_pos.to_csv(os.path.join(spatial_path, "tissue_positions_list.csv"),index=False,header=False)
    if verbose:
        print(f"The new MSI coordinate file has been saved, containing {len(msi_tissue_pos)} spots/pixels.")

    # Update spot diameter based on MSI center-to-center distance
    # Assuming pxl_row_in_fullres and pxl_col_in_fullres are spatial coordinates
    distances = pdist(msi_tissue_pos[["pxl_row_in_fullres", "pxl_col_in_fullres"]])
    c_c_dist = distances[distances > 0].min()
    msi_json["spot_diameter_fullres"] = c_c_dist
    msi_json["tissue_lowres_scalef"] = msi_json["tissue_hires_scalef"]

    # Write new scale factor JSON (if data available)
    if (msi_peak_data_path and visium_sf_json_path) or msi_json:
        if verbose:
            print("Creating JSON file...")
        with open(os.path.join(spatial_path, "scalefactors_json.json"), "w") as f:
            json.dump(msi_json, f, indent=4)

    # Copy image into spatial folder
    if verbose:
        print("Copying tissue image...")
    shutil.copyfile(msi_image_path, os.path.join(spatial_path, "tissue_hires_image.png"))
    shutil.copyfile(msi_image_path, os.path.join(spatial_path, "tissue_lowres_image.png"))

    if verbose:
        print("Spaceranger 'spatial' folder content ready!")

    # Read MSI peak data file
    if verbose:
        print("Reading MSI peak data file...")
    msi_peaks = pd.read_csv(msi_peak_data_path, header=0)
    msi_peaks.set_index('spot_id', drop=True, inplace = True)

    # Create feature IDs
    msi_peaks_features = pd.DataFrame(
        {"id1": msi_peaks.columns}
    )

    # Create barcode IDs
    msi_peaks_barcodes = msi_spot_prefix + "_" + sample_name + "_" + msi_peaks.index.astype(str)

    # Prepare MSI peak data matrix
    if verbose:
        print("Preparing MSI peak data...")
    
    msi_peaks_matrix = msi_peaks.transpose().copy()  # Select and transpose data
    msi_peaks_matrix.index = msi_peaks_features["id1"].tolist()  # Set row names
    msi_peaks_matrix.columns = msi_peaks_barcodes.tolist()  # Set column names

    # Convert to sparse matrix (optional)
    # You might need to install `scikit-learn` for sparse matrix functionality
    msi_peaks_mtx = csr_matrix(msi_peaks_matrix)

    # Write output files
    if verbose:
        print("Writing new feature-barcode matrix data...")
        
    msi_peaks_features.to_csv(
        os.path.join(filtered_path, "features.tsv.gz"),
        sep="\t",
        index=False,
        compression="gzip",
        header=False
    )
    pd.DataFrame({"barcodes": msi_peaks_barcodes}).to_csv(
        os.path.join(filtered_path, "barcodes.tsv.gz"),
        sep="\t",
        index=False,
        compression="gzip",
        header=False
    )

    # Write matrix data (modify based on chosen sparse matrix format)
    mmwrite(os.path.join(filtered_path, "matrix.mtx"),msi_peaks_mtx)

    with open(os.path.join(filtered_path, "matrix.mtx"), 'rb') as src, gzip.open(os.path.join(filtered_path, "matrix.mtx.gz"), 'wb') as dst:
        dst.writelines(src)

    if verbose:
        print("Spaceranger 'filtered_feature_bc_matrix' folder content ready!")

    if verbose:
        print("Writing filtered_feature_bc_matrix.h5 file...")

    # create h5 object
    msi_peaks_mtx_csr = msi_peaks_mtx.tocoo()
    msi_peaks_mtx_csr = coo_array((msi_peaks_mtx_csr.data, (msi_peaks_mtx_csr.col, msi_peaks_mtx_csr.row)), shape=(msi_peaks_mtx_csr.shape[1],msi_peaks_mtx_csr.shape[0])).tocsr()
    hf = h5py.File(os.path.join(output_path, "filtered_feature_bc_matrix.h5"), 'w')
    g1 = hf.create_group('matrix')
    g1.create_dataset('data', data=msi_peaks_mtx_csr.data)
    g1.create_dataset('indices', data=msi_peaks_mtx_csr.indices)
    g1.create_dataset('indptr', data=msi_peaks_mtx_csr.indptr)
    g1.create_dataset('shape', data=msi_peaks_mtx.shape)
    g1.create_dataset('barcodes', data=msi_peaks_barcodes.tolist(), dtype=h5py.string_dtype('utf-8'))
    g2 = g1.create_group('features')
    g2.create_dataset('_all_tag_keys',data=['genome'],dtype=h5py.string_dtype('utf-8'))
    g2.create_dataset('id', data=msi_peaks_features["id1"].tolist())
    g2.create_dataset('name', data=msi_peaks_features["id1"].tolist())
    g2.create_dataset('feature_type', data=len(msi_peaks_features["id1"].tolist()) * ['Gene Expression'])
    g2.create_dataset('genome', data=len(msi_peaks_features["id1"].tolist()) * ['unknown'])
    hf.close()

    if verbose:
        print("Spaceranger 'filtered_feature_bc_matrix.h5' file ready!")


    if verbose:
        print("Mock Space Ranger output completed for sample at path "+output_path)


if __name__ == "__main__":
    sample = snakemake.params['sample']
    create_mock_spaceranger(
                        sample,
                        'output/'+sample+'/spaceranger/',
                        snakemake.input[1],
                        snakemake.input[0],
                        visium_sf_json_path = "input/"+sample+"/visium/spatial/scalefactors_json.json",
                        msi_peak_data_path = "input/"+sample+"/msi/MSI_intensities.csv",
                        verbose=False)
