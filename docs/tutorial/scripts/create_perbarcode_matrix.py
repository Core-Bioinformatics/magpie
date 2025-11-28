import pandas as pd
import scanpy
from sklearn.neighbors import KNeighborsClassifier
from scipy.spatial import distance_matrix
import os
import json
from scipy.sparse import csr_matrix, coo_array
import shutil
from scipy.io import mmwrite
import gzip
import h5py
from pathlib import Path
import numpy as np

# ------------------------------------------------------------
# Function to estimate the minimal inter-spot distance
# This is used to define the extended radius for grouping MSI pixels
# ------------------------------------------------------------
def estimate_interspot_distance(data):
    # define a dummy class vector for kNN
    classes = [1] * data.shape[0]
    
    # KNN with 1 neighbor to get nearest distances
    knn = KNeighborsClassifier(n_neighbors=1, n_jobs=1, metric='euclidean')
    knn.fit(data, classes)
    
    # get distances to the two nearest neighbors (self included)
    dist, ind = knn.kneighbors(data, 2)
    
    # take the distance to the nearest neighbor that is NOT self
    distances = [dist[i][1] for i in range(len(dist))]
    
    # return half of the minimal distance, used as grouping radius
    return min(distances) / 2

# ------------------------------------------------------------
# Function to assign MSI pixels to Visium spots based on distance
# ------------------------------------------------------------
def group_msi_perspot(msi_obj,
                      visium_allcoords,
                      visium_index,
                      scale_factors,
                      sample_name,
                      verbose=True):

    # determine max distance for grouping based on user choice
    if snakemake.params['radius_to_use'] == 'visium_expanded':
        min_distance = estimate_interspot_distance(visium_allcoords)
        print("Using maximum distance between spots and pixels: " + str(min_distance))
    elif snakemake.params['radius_to_use'] == 'visium':
        spot_diameter_px = scale_factors['spot_diameter_fullres']
        min_distance = spot_diameter_px / 2
        print("Using maximum distance between spots and pixels: " + str(min_distance))
    else:
        spot_diameter_px = scale_factors['spot_diameter_fullres']
        # take max of MSI or Visium inter-spot distances plus spot diameter
        min_distance = np.max([
            estimate_interspot_distance(visium_allcoords),
            estimate_interspot_distance(msi_obj.obsm['spatial'])
        ]) + spot_diameter_px / 2
        print("Using maximum distance between spots and pixels: " + str(min_distance))

    # ------------------------------------------------------------
    # Compute distance matrix between all Visium spots and MSI pixels
    # ------------------------------------------------------------
    if verbose:
        print("Calculating distances between MSI pixels and Visium spots...")

    spot_distances = pd.DataFrame(
        distance_matrix(visium_allcoords, msi_obj.obsm['spatial']),
        columns=msi_obj.obs.index
    )
    spot_distances['visium_spot'] = visium_index

    # ------------------------------------------------------------
    # Reshape distance matrix to long format and filter by min_distance
    # ------------------------------------------------------------
    if verbose:
        print("Grouping MSI pixels by Visium spot...")

    spot_distances_long = pd.wide_to_long(
        spot_distances,
        i="visium_spot",
        stubnames='MSI_' + sample_name,
        j='MSI_spot',
        sep='_',
        suffix=r".+"
    )
    spot_distances_long.columns = ['distance']

    # keep only MSI pixels within the extended Visium radius
    close_points = spot_distances_long[spot_distances_long['distance'] < min_distance].reset_index()
    close_points['MSI_spot'] = ['MSI_' + sample_name + "_" + str(spot_id) for spot_id in close_points['MSI_spot']]

    # save matched coordinates
    close_points.to_csv('output/' + sample_name + '/matched_coords.csv')

    # count matched and unmatched MSI/Visium points
    matched_msi = close_points['MSI_spot'].nunique()
    matched_visium = close_points['visium_spot'].nunique()
    unmatched_visium = len(set(visium_index) - set(close_points['visium_spot']))
    unmatched_msi = len(set(msi_obj.obs.index) - set(close_points['MSI_spot']))

    if verbose:
        print("Number of MSI pixels matching to Visium spots: " + str(matched_msi))
        print("Number of Visium spots matching to MSI pixels: " + str(matched_visium))
        print("Number of MSI pixels not matching to Visium spots: " + str(unmatched_msi))
        print("Number of Visium spots not matching to MSI pixels: " + str(unmatched_visium))

    return close_points

# ------------------------------------------------------------
# Helper function to compute mean of a group (fallback method)
# ------------------------------------------------------------
def group_mean(group):
    return np.mean(group, axis=0)

# ------------------------------------------------------------
# Aggregate MSI intensities per Visium spot
# ------------------------------------------------------------
def create_mean_intensity_table(
        msi_obj,
        visium_allcoords,
        scale_factors,
        sample_name='sample1',
        verbose=True):
    
    # group MSI pixels to Visium spots
    close_points = group_msi_perspot(
        msi_obj,
        visium_allcoords,
        visium_allcoords.index,
        scale_factors,
        sample_name,
        verbose=verbose
    )

    # convert sparse MSI matrix to pandas DataFrame
    msi_data = pd.DataFrame.sparse.from_spmatrix(msi_obj.X)
    msi_data.index = msi_obj.obs.index
    msi_data.columns = msi_obj.var['gene_ids']

    # set MSI_spot as index for merging
    close_points.set_index('MSI_spot', inplace=True, drop=False)
    close_points.iloc[:, :2].to_csv('output/' + sample_name + '/matched_Visium_MSI_IDs.csv', index=False)

    # merge MSI intensities with matched pixels
    intensity_matrix = close_points.merge(msi_data, left_index=True, right_index=True)

    # ------------------------------------------------------------
    # Aggregate intensities per Visium spot using chosen method
    # ------------------------------------------------------------
    if verbose:
        print("Aggregating MSI pixels per Visium spot by " + snakemake.params['agg_fn'] + "...")

    if snakemake.params['agg_fn'] == 'sum':
        # simple sum of intensities
        intensity_matrix_mean = intensity_matrix.drop(labels=['distance', 'MSI_spot'], axis=1).groupby(['visium_spot']).sum()

    elif snakemake.params['agg_fn'] == 'weighted_average':
        # weight pixels by inverse distance
        epsilon = 1e-6
        intensity_matrix['raw_weight'] = 1 / (intensity_matrix['distance'] + epsilon)
        intensity_matrix['weight'] = (
            intensity_matrix.groupby('visium_spot')['raw_weight']
            .transform(lambda w: w / w.sum())
        )

        weighted_intensities = intensity_matrix.drop(columns=['distance', 'raw_weight', 'MSI_spot']).copy()
        peak_cols = [col for col in weighted_intensities.columns if col not in ['visium_spot', 'weight']]
        weighted_intensities[peak_cols] = weighted_intensities[peak_cols].multiply(intensity_matrix['weight'], axis=0)

        numerator = weighted_intensities.groupby('visium_spot').sum()
        denominator = intensity_matrix.groupby('visium_spot')['weight'].sum()
        intensity_matrix_mean = numerator.div(denominator, axis=0).drop(columns=['weight'])

    else:
        # fallback: simple mean using group_mean function
        intensity_matrix_mean = intensity_matrix.drop(labels=['distance', 'MSI_spot'], axis=1)\
            .groupby(['visium_spot']).apply(lambda x: group_mean(x.values))

    # ensure columns are in same order as MSI peaks
    intensity_matrix_mean_df = pd.DataFrame.from_dict(
        dict(zip(intensity_matrix_mean.index, intensity_matrix_mean.values))
    ).transpose()
    intensity_matrix_mean_df.columns = msi_obj.var['gene_ids']

    return intensity_matrix_mean_df

# ------------------------------------------------------------
# Create a mock Space Ranger-style output for aggregated MSI
# ------------------------------------------------------------
def create_mock_spaceranger_aggregated_intensity(
    mean_intensity_table,
    output_folder_name="mock_spaceranger_aggregated",
    msi_image_path="MSI_HE.jpg",
    visium_dir=None,
    verbose=True
):
    """
    Creates mock Space Ranger formatted output data aggregated to give matching observations between modalities.

    Args:
        mean_intensity_table: DataFrame with aggregated MSI intensity per Visium spot
        output_folder_name (str, optional): Name of the output folder. Defaults to "mock_spaceranger_aggregated".
        msi_image_path: Path to the MSI image file. Defaults to "MSI_HE.jpg".
        visium_dir: Path to the Visium spaceranger output.
        verbose (bool, optional): Print informational messages during execution. Defaults to True.

    Returns:
        None
    """

    # Create output folders if they don't already exist
    output_path = os.path.join(output_folder_name)
    Path(output_folder_name).mkdir(parents=True, exist_ok=True)
    spatial_path = os.path.join(output_path, "spatial")
    Path(spatial_path).mkdir(parents=True, exist_ok=True)
    filtered_path = os.path.join(output_path, "filtered_feature_bc_matrix")
    Path(filtered_path).mkdir(parents=True, exist_ok=True)

    # read Visium scale factors to copy to MSI
    with open(visium_dir + "/spatial/scalefactors_json.json", "r") as f:
        msi_json = json.load(f)

    # read Visium coordinates
    if verbose:
        print("Copying Visium coordinates which map to some MSI pixels...")

    if os.path.isfile(visium_dir + '/spatial/tissue_positions.csv'):
        visium_coords = pd.read_csv(visium_dir + '/spatial/tissue_positions.csv', index_col=0)
    else:
        visium_coords = pd.read_csv(visium_dir + '/spatial/tissue_positions_list.csv', header=None, index_col=0)

    msi_tissue_pos = visium_coords.loc[list(mean_intensity_table.index), :]

    # write tissue_positions.csv
    msi_tissue_pos.to_csv(os.path.join(spatial_path, "tissue_positions_list.csv"), index=True, header=False)

    if verbose:
        print(f"The new MSI coordinate file has been saved, containing {len(msi_tissue_pos)} spots/pixels.")

    # create updated scale factor JSON
    if verbose:
        print("Creating JSON file...")
    msi_json["tissue_lowres_scalef"] = msi_json["tissue_hires_scalef"]
    with open(os.path.join(spatial_path, "scalefactors_json.json"), "w") as f:
        json.dump(msi_json, f, indent=4)

    # copy tissue image
    if verbose:
        print("Copying tissue image...")
    shutil.copyfile(msi_image_path, os.path.join(spatial_path, "tissue_hires_image.png"))
    shutil.copyfile(msi_image_path, os.path.join(spatial_path, "tissue_lowres_image.png"))

    if verbose:
        print("Spaceranger 'spatial' folder content ready!")

    # prepare feature and barcode files
    if verbose:
        print("Reading MSI peak data file...")
    msi_peaks = mean_intensity_table
    msi_peaks_features = pd.DataFrame({"id1": msi_peaks.columns})
    msi_peaks_features["id2"] = msi_peaks_features["id1"]
    msi_peaks_barcodes = msi_peaks.index.astype(str)

    # prepare peak data matrix
    if verbose:
        print("Preparing MSI peak data...")
    msi_peaks_matrix = msi_peaks.transpose()
    msi_peaks_matrix.index = msi_peaks_features["id1"].tolist()
    msi_peaks_matrix.columns = msi_peaks_barcodes.tolist()
    msi_peaks_mtx = csr_matrix(msi_peaks_matrix)

    # save features, barcodes, and matrix
    if verbose:
        print("Writing new feature-barcode matrix data...")
    msi_peaks_features.to_csv(os.path.join(filtered_path, "features.tsv.gz"),
                              sep="\t", index=False, compression="gzip", header=False)
    msi_peaks_barcodes.to_frame(index=False).to_csv(os.path.join(filtered_path, "barcodes.tsv.gz"),
                                                    sep="\t", index=False, compression="gzip", header=False)
    with gzip.open(os.path.join(filtered_path, "matrix.mtx.gz"), 'wb') as f:
        mmwrite(f, msi_peaks_mtx)

    if verbose:
        print("Spaceranger 'filtered_feature_bc_matrix' folder content ready!")

    # create HDF5 object for matrix
    if verbose:
        print("Writing filtered_feature_bc_matrix.h5 file...")
    msi_peaks_mtx_csr = msi_peaks_mtx.tocoo()
    msi_peaks_mtx_csr = coo_array((msi_peaks_mtx_csr.data,
                                   (msi_peaks_mtx_csr.col, msi_peaks_mtx_csr.row)),
                                  shape=(msi_peaks_mtx_csr.shape[1], msi_peaks_mtx_csr.shape[0])).tocsr()
    hf = h5py.File(os.path.join(output_path, "filtered_feature_bc_matrix.h5"), 'w')
    g1 = hf.create_group('matrix')
    g1.create_dataset('library_ids', data="utf-8", dtype=h5py.special_dtype(vlen=str))
    g1.create_dataset('data', data=msi_peaks_mtx_csr.data)
    g1.create_dataset('indices', data=msi_peaks_mtx_csr.indices)
    g1.create_dataset('indptr', data=msi_peaks_mtx_csr.indptr)
    g1.create_dataset('shape', data=msi_peaks_mtx.shape)
    g1.create_dataset('barcodes', data=msi_peaks_barcodes.tolist(), dtype=h5py.special_dtype(vlen=str))
    g2 = g1.create_group('features')
    g2.create_dataset('_all_tag_keys', data='genome', dtype=h5py.special_dtype(vlen=str))
    g2.create_dataset('id', data=msi_peaks_features["id1"].tolist(), dtype=h5py.special_dtype(vlen=str))
    g2.create_dataset('name', data=msi_peaks_features["id1"].tolist(), dtype=h5py.special_dtype(vlen=str))
    g2.create_dataset('feature_type', data=len(msi_peaks_features["id1"].tolist()) * ['Gene Expression'], dtype=h5py.special_dtype(vlen=str))
    g2.create_dataset('genome', data=len(msi_peaks_features["id1"].tolist()) * ['unknown'], dtype=h5py.special_dtype(vlen=str))
    hf.close()

    if verbose:
        print("Spaceranger 'filtered_feature_bc_matrix.h5' file ready!")
        print("Mock Space Ranger output completed for sample at path " + output_path)

# ------------------------------------------------------------
# Main script execution
# ------------------------------------------------------------
def main():
    sample = snakemake.params['sample']

    # read MSI object using scanpy
    msi_obj = scanpy.read_visium('output/' + sample + '/spaceranger/', library_id='myLib')

    # read Visium coordinates
    if os.path.isfile('input/' + sample + '/visium/spatial/tissue_positions.csv'):
        visium_allcoords = pd.read_csv("input/" + sample + '/visium/spatial/tissue_positions.csv', index_col=0)
    else:
        visium_allcoords = pd.read_csv("input/" + sample + '/visium/spatial/tissue_positions_list.csv', header=None, index_col=0)

    # assign column names for Visium coordinates
    visium_allcoords.columns = ['in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']

    # filter to tissue spots only if requested
    if snakemake.params['only_within_tissue']:
        visium_allcoords = visium_allcoords[visium_allcoords['in_tissue'] == 1]

    # load Visium scale factors
    with open("input/" + sample + "/visium/spatial/scalefactors_json.json", "r") as f:
        visium_json = json.load(f)

    # create aggregated mean intensity table
    mean_intensity = create_mean_intensity_table(
        msi_obj,
        visium_allcoords.iloc[:, [4, 3]],  # pxl_col_in_fullres, pxl_row_in_fullres
        visium_json,
        sample,
        verbose=snakemake.params['verbose']
    )

    # create mock Space Ranger output
    create_mock_spaceranger_aggregated_intensity(
        mean_intensity,
        'output/' + sample + '/spaceranger_aggregated/',
        'output/' + sample + '/spaceranger/spatial/tissue_hires_image.png',
        visium_dir="input/" + sample + "/visium/",
        verbose=snakemake.params['verbose']
    )

if __name__ == "__main__":
    main()
