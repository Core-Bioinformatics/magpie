import os
import json
from collections import namedtuple
from scipy.spatial.distance import pdist
from scipy.spatial import cKDTree
import numpy as np
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
    Generate a mock Spaceranger output folder structure from MSI data.

    This function creates:
        - spatial/ folder with:
            - tissue image(s)
            - tissue_positions_list.csv
            - scalefactors_json.json
        - filtered_feature_bc_matrix/ folder with:
            - features.tsv.gz
            - barcodes.tsv.gz
            - matrix.mtx.gz
        - filtered_feature_bc_matrix.h5 (10x-style HDF5 matrix)

    It scales MSI coordinates to Visium units using an existing Visium
    scalefactors JSON file, rewrites the MSI peak matrix in 10x format,
    and synthesizes a minimal valid Spaceranger output suitable for
    downstream tools expecting Visium-format input.

    Parameters
    ----------
    sample_name : str
        Name of the sample, used to prefix barcode IDs.
    output_folder_name : str
        Path where the mock Spaceranger output will be written.
    msi_image_path : str
        Path to the MSI image (used as mock H&E image).
    msi_coord_fname : str
        Path to the MSI coordinate table (with spot_id, x, y).
    msi_spot_prefix : str
        Prefix added to barcodes to distinguish MSI-derived spots.
    visium_sf_json_path : str
        Path to the Visium `scalefactors_json.json` file to mimic Visium scaling.
    msi_peak_data_path : str
        Path to MSI peak intensity table, indexed by spot_id.
    verbose : bool
        Whether to print progress statements.
    """

    # --- Create necessary output directories ---
    output_path = os.path.join(output_folder_name)
    os.makedirs(output_path, exist_ok=True)
    spatial_path = os.path.join(output_path, "spatial")
    os.makedirs(spatial_path, exist_ok=True)
    filtered_path = os.path.join(output_path, "filtered_feature_bc_matrix")
    os.makedirs(filtered_path, exist_ok=True)

    # --- Load Visium scalefactor JSON to reuse scaling parameters ---
    with open(visium_sf_json_path, "r") as f:
        st_json = json.load(f)
    scale_factor = st_json['tissue_hires_scalef']
    msi_json = st_json

    # --- Load MSI coordinates and scale them to Visium coordinate units ---
    if verbose:
        print("Reading MSI coordinate file...")

    msi_coords = pd.read_csv(msi_coord_fname)
    # Scale MSI coordinates by Visium scale factor
    msi_coords['x'] = msi_coords['x']/scale_factor
    msi_coords['y'] = msi_coords['y']/scale_factor
    msi_coords['x_2'] = msi_coords['x']
    msi_coords['y_2'] = msi_coords['y']
    msi_coords.set_index('spot_id', drop=False, inplace=True)

    # --- Mark all MSI spots as 'in tissue' for Spaceranger compatibility ---
    msi_coords["in_tissue"] = 1

    # Reformat into Visium-style tissue_positions_list.csv structure
    msi_tissue_pos = msi_coords.rename(
        columns={
            "spot_id":"barcode",
            'x_2': "array_col",
            'y_2': "array_row",
            'x': "pxl_col_in_fullres",
            'y': "pxl_row_in_fullres",
        }
    )

    # --- Create descriptive barcodes consistent with Visium format ---
    msi_tissue_pos["barcode"] = (
        msi_spot_prefix + "_" + sample_name + "_" + msi_tissue_pos["barcode"].astype(str)
    )

    # Keep only columns required by Visium format
    msi_tissue_pos = msi_tissue_pos[['barcode','in_tissue','array_row','array_col','pxl_row_in_fullres','pxl_col_in_fullres']]

    # --- Write tissue positions file ---
    msi_tissue_pos.to_csv(os.path.join(spatial_path, "tissue_positions_list.csv"),
                          index=False, header=False)

    if verbose:
        print(f"The new MSI coordinate file has been saved, containing {len(msi_tissue_pos)} spots/pixels.")

    # --- Compute mock "spot diameter" based on nearest-neighbor spacing ---
    coords = msi_tissue_pos[["pxl_row_in_fullres", "pxl_col_in_fullres"]].values
    tree = cKDTree(coords)
    dists, idxs = tree.query(coords, k=2)  # k=2 to get nearest neighbor
    min_dist = np.min(dists[:, 1])

    # Adjust JSON scalefactor fields
    msi_json["spot_diameter_fullres"] = min_dist
    msi_json["tissue_lowres_scalef"] = msi_json["tissue_hires_scalef"]

    # --- Write updated JSON file ---
    if (msi_peak_data_path and visium_sf_json_path) or msi_json:
        if verbose:
            print("Creating JSON file...")
        with open(os.path.join(spatial_path, "scalefactors_json.json"), "w") as f:
            json.dump(msi_json, f, indent=4)

    # --- Copy MSI image as both high-res and low-res Visium images ---
    if verbose:
        print("Copying tissue image...")
    shutil.copyfile(msi_image_path, os.path.join(spatial_path, "tissue_hires_image.png"))
    shutil.copyfile(msi_image_path, os.path.join(spatial_path, "tissue_lowres_image.png"))

    if verbose:
        print("Spaceranger 'spatial' folder content ready!")

    # --- Load MSI peak intensity table ---
    if verbose:
        print("Reading MSI peak data file...")
    msi_peaks = pd.read_csv(msi_peak_data_path, header=0)
    msi_peaks.set_index('spot_id', drop=True, inplace=True)

    # Feature list = peak names
    msi_peaks_features = pd.DataFrame({"id1": msi_peaks.columns})

    # Barcode list must match tissue_positions_list.csv
    msi_peaks_barcodes = msi_spot_prefix + "_" + sample_name + "_" + msi_peaks.index.astype(str)

    # --- Prepare peak data matrix in 10x (features x barcodes) orientation ---
    if verbose:
        print("Preparing MSI peak data...")

    msi_peaks_matrix = msi_peaks.transpose()
    msi_peaks_matrix.index = msi_peaks_features["id1"].values
    msi_peaks_matrix.columns = msi_peaks_barcodes.values

    # Convert dense data to CSR sparse matrix (10x compatible)
    msi_peaks_mtx = csr_matrix(msi_peaks_matrix)

    # --- Write 10x-format TSVs and MTX ---
    if verbose:
        print("Writing new feature-barcode matrix data...")

    msi_peaks_features.to_csv(
        os.path.join(filtered_path, "features.tsv.gz"),
        sep="\t", index=False, compression="gzip", header=False
    )

    msi_peaks_barcodes.to_frame(index=False).to_csv(
        os.path.join(filtered_path, "barcodes.tsv.gz"),
        sep="\t", index=False, compression="gzip", header=False
    )

    # Write MTX matrix
    with gzip.open(os.path.join(filtered_path, "matrix.mtx.gz"), 'wb') as f:
        mmwrite(f, msi_peaks_mtx)

    if verbose:
        print("Spaceranger 'filtered_feature_bc_matrix' folder content ready!")

    # --- Write HDF5 version of the MSI peak matrix ---
    if verbose:
        print("Writing filtered_feature_bc_matrix.h5 file...")

    # Reorder to CSR in (barcodes x features) orientation as required by 10x spec
    mtx_coo = msi_peaks_mtx.tocoo()
    mtx_swapped = coo_array(
        (mtx_coo.data, (mtx_coo.col, mtx_coo.row)),
        shape=(mtx_coo.shape[1], mtx_coo.shape[0])
    ).tocsr()

    hf = h5py.File(os.path.join(output_path, "filtered_feature_bc_matrix.h5"), 'w')
    g1 = hf.create_group('matrix')
    g1.create_dataset('library_ids', data=["utf-8"], dtype=h5py.special_dtype(vlen=str))
    g1.create_dataset('data', data=mtx_swapped.data)
    g1.create_dataset('indices', data=mtx_swapped.indices)
    g1.create_dataset('indptr', data=mtx_swapped.indptr)
    g1.create_dataset('shape', data=msi_peaks_mtx.shape)
    g1.create_dataset('barcodes', data=msi_peaks_barcodes.tolist(), dtype=h5py.special_dtype(vlen=str))

    # Feature metadata group
    g2 = g1.create_group('features')
    g2.create_dataset('_all_tag_keys', data=['genome'] * len(msi_peaks_features["id1"]), dtype=h5py.special_dtype(vlen=str))
    g2.create_dataset('id', data=msi_peaks_features["id1"].tolist(), dtype=h5py.special_dtype(vlen=str))
    g2.create_dataset('name', data=msi_peaks_features["id1"].tolist(), dtype=h5py.special_dtype(vlen=str))
    g2.create_dataset('feature_type', data=len(msi_peaks_features) * ['Gene Expression'], dtype=h5py.special_dtype(vlen=str))
    g2.create_dataset('genome', data=len(msi_peaks_features) * ['unknown'], dtype=h5py.special_dtype(vlen=str))

    hf.close()

    if verbose:
        print("Spaceranger 'filtered_feature_bc_matrix.h5' file ready!")

    # --- Completion message ---
    if verbose:
        print("Mock Space Ranger output completed for sample at path " + output_path)


if __name__ == "__main__":
    sample = snakemake.params['sample']
    create_mock_spaceranger(
        sample,
        'output/'+sample+'/spaceranger/',
        snakemake.input[1],
        snakemake.input[0],
        visium_sf_json_path="input/"+sample+"/visium/spatial/scalefactors_json.json",
        msi_peak_data_path="input/"+sample+"/msi/MSI_intensities.csv",
        verbose=False
    )
