import pandas as pd
import scanpy
from sklearn.neighbors import KNeighborsClassifier
from scipy.spatial import distance_matrix
import os
import json
from collections import namedtuple
from scipy.spatial.distance import pdist
from scipy.sparse import csr_matrix,coo_array
import shutil
from scipy.io import mmwrite
import gzip
import h5py
from pathlib import Path
import numpy as np

# calculate actual distance between spots to find extended spot radius
def estimate_interspot_distance(data):

    # define kNN model with Euclidean distance
    classes = [1]*data.shape[0]
    knn = KNeighborsClassifier(n_neighbors=1,n_jobs=1,metric='euclidean')
    knn.fit(data, classes)
    # find all minimum distances to neighbours
    dist, ind = knn.kneighbors(data,2)
    distances = [dist[i][1] for i in range(len(dist))]
    return(min(distances)/2)

# group MSI coordinates which fall within a Visium spot extended radius
def group_msi_perspot(msi_obj,
                      visium_allcoords,
                      visium_index,
                      sample_name,
                      verbose = True):

    # get minimum distance between spots to find extended Visium radius
    min_distance = estimate_interspot_distance(visium_allcoords)
    # find distance from all MSI pixels to all Visium spots
    if verbose:
        print("Calculating distances between MSI pixels and Visium spots...")

    spot_distances = pd.DataFrame(distance_matrix(visium_allcoords,msi_obj.obsm['spatial']),columns=msi_obj.obs.index)
    spot_distances['visium_spot'] = visium_index

    # convert to long format and filter to those within the extended Visium radius
    if verbose:
        print("Grouping MSI pixels by Visium spot...")

    spot_distances_long = pd.wide_to_long(spot_distances, i="visium_spot",stubnames='MSI_' + sample_name,j='MSI_spot', sep='_',suffix=r".+")
    spot_distances_long.columns = ['distance']
    close_points = spot_distances_long[spot_distances_long['distance'] < min_distance].reset_index()
    close_points['MSI_spot'] = ['MSI_'+ sample_name + "_" +str(spot_id) for spot_id in close_points['MSI_spot']]
    matched_msi = close_points['MSI_spot'].nunique()
    matched_visium = close_points['visium_spot'].nunique()
    unmatched_visium = len(list(set([x for x in visium_index if x not in list(close_points['visium_spot'])])))
    unmatched_msi = len(list(set([x for x in msi_obj.obs.index if x not in list(close_points['MSI_spot'])])))
    if verbose:
        print("Number of MSI pixels matching to Visium spots: "+str(matched_msi))
        print("Number of Visium spots matching to MSI pixels: "+str(matched_visium))
        print("Number of MSI pixels not matching to Visium spots: "+str(unmatched_msi))
        print("Number of Visium spots not matching to MSI pixels: "+str(unmatched_visium))

    return(close_points)

# helper function to get mean for each group
def group_mean(group):
    return np.mean(group, axis=0)

# combine MSI pixels in one Visium spot and take average
def create_mean_intensity_table(
        msi_obj,
        visium_allcoords,
        sample_name = 'sample1',
        verbose = True):

    """
    Aggregates MSI intensities per Visium spot.

    Args:
        msi_obj: Path to Space Ranger-style output from stage 2 of snakemake pipeline.
        visium_allcoords: DataFrame of Visium coordinates as read from Space Ranger output (pxl_col_in_fullres and pxl_row_in_fullres). The row index should be the spot barcodes.
        sample_name: Name of sample which will be used to rename MSI pixels for uniqueness. Defaults to "sample1".
        verbose (bool, optional): Print informational messages during execution. Defaults to True.

    Returns:
        None
    """
    # get all points grouped by Visium spot
    close_points = group_msi_perspot(msi_obj,
                                     visium_allcoords,
                                     visium_allcoords.index,
                                     sample_name,
                                     verbose = verbose)
    msi_data = msi_obj.X
    # convert to sparse format
    msi_data = pd.DataFrame.sparse.from_spmatrix(msi_data)
    msi_data.index = msi_obj.obs.index
    msi_data.columns = msi_obj.var['gene_ids']
    close_points.set_index('MSI_spot',inplace=True,drop=False)
    close_points.iloc[:,:2].to_csv('output/'+sample_name+'/matched_Visium_MSI_IDs.csv', index=False)
    # get MSI intensities for all grouped points
    intensity_matrix = close_points.merge(msi_data,left_index=True,right_index=True)

    # group by Visium spot and take mean for each peak
    if verbose:
        print("Aggregating MSI pixels per Visium spot by "+snakemake.params['agg_fn']+"...")
    if snakemake.params['agg_fn'] == 'sum':
        intensity_matrix_mean = intensity_matrix.drop(labels=['distance','MSI_spot'],axis=1).groupby(['visium_spot']).sum()
    else :
        intensity_matrix_mean = intensity_matrix.drop(labels=['distance','MSI_spot'],axis=1).groupby(['visium_spot']).apply(lambda x: group_mean(x.values))
    
    intensity_matrix_mean_df = pd.DataFrame.from_dict(dict(zip(intensity_matrix_mean.index, intensity_matrix_mean.values))).transpose()
    intensity_matrix_mean_df.columns = msi_obj.var['gene_ids']

    return(intensity_matrix_mean_df)

# create spaceranger-style output
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

    # get Visium scale factors and copy to MSI
    with open(visium_dir+"/spatial/scalefactors_json.json", "r") as f:
        msi_json = json.load(f)

    # Get Visium coordinates which had an MSI pixel mapped to it
    if verbose:
        print("Copying Visium coordinates which map to some MSI pixels...")
    
    visium_coords = pd.read_csv(visium_dir+"/spatial/tissue_positions_list.csv",header=None,index_col=0)
    msi_tissue_pos = visium_coords.loc[list(mean_intensity_table.index),:]

    # Write tissue_positions.csv
    msi_tissue_pos.to_csv(os.path.join(spatial_path, "tissue_positions_list.csv"),index=True,header=False)

    if verbose:
        print(f"The new MSI coordinate file has been saved, containing {len(msi_tissue_pos)} spots/pixels.")

    # Write new scale factor JSON
    if verbose:
        print("Creating JSON file...")
    msi_json["tissue_lowres_scalef"] = msi_json["tissue_hires_scalef"]
    with open(os.path.join(spatial_path, "scalefactors_json.json"), "w") as f:
        json.dump(msi_json, f, indent=4)

    # Copy image into spatial folder
    if verbose:
        print("Copying tissue image...")
    shutil.copyfile(msi_image_path, os.path.join(spatial_path, "tissue_hires_image.png"))
    shutil.copyfile(msi_image_path, os.path.join(spatial_path, "tissue_lowres_image.png"))

    if verbose:
        print("Spaceranger 'spatial' folder content ready!")

    # Get MSI peak data file and create h5 object

    if verbose:
        print("Reading MSI peak data file...")
    msi_peaks = mean_intensity_table
   
    # Create feature IDs
    msi_peaks_features = pd.DataFrame(
        {"id1": msi_peaks.columns}
    )
    msi_peaks_features["id2"] = msi_peaks_features["id1"]

    # Create barcode IDs
    msi_peaks_barcodes = msi_peaks.index.astype(str)

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
        
    # Save feature and barcode files to filtered_feature_bc_matrix folder
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
    g1.create_dataset('library_ids', data="utf-8",dtype=h5py.special_dtype(vlen=str))
    g1.create_dataset('data', data=msi_peaks_mtx_csr.data)
    g1.create_dataset('indices', data=msi_peaks_mtx_csr.indices)
    g1.create_dataset('indptr', data=msi_peaks_mtx_csr.indptr)
    g1.create_dataset('shape', data=msi_peaks_mtx.shape)
    g1.create_dataset('barcodes', data=msi_peaks_barcodes.tolist(),dtype=h5py.special_dtype(vlen=str))
    g2 = g1.create_group('features')
    g2.create_dataset('_all_tags_keys',data='genome',dtype=h5py.special_dtype(vlen=str))
    g2.create_dataset('id', data=msi_peaks_features["id1"].tolist(),dtype=h5py.special_dtype(vlen=str))
    g2.create_dataset('name', data=msi_peaks_features["id1"].tolist(),dtype=h5py.special_dtype(vlen=str))
    g2.create_dataset('feature_type', data=len(msi_peaks_features["id1"].tolist()) * ['Gene Expression'],dtype=h5py.special_dtype(vlen=str))
    g2.create_dataset('genome', data=len(msi_peaks_features["id1"].tolist()) * ['unknown'],dtype=h5py.special_dtype(vlen=str))
    hf.close()

    if verbose:
        print("Spaceranger 'filtered_feature_bc_matrix.h5' file ready!")

    if verbose:
        print("Mock Space Ranger output completed for sample at path "+output_path)

def main():
    sample = snakemake.params['sample']
    msi_obj = scanpy.read_visium('output/'+sample+'/spaceranger/',library_id='myLib')
    if os.path.isfile('input/'+sample+'/visium/spatial/tissue_positions.csv'):
        visium_allcoords = pd.read_csv("input/"+sample+'/visium/spatial/tissue_positions.csv',index_col=0)
    else :
        visium_allcoords = pd.read_csv("input/"+sample+'/visium/spatial/tissue_positions_list.csv',header=None,index_col=0)
    visium_allcoords.columns = ['in_tissue','array_row','array_col','pxl_row_in_fullres','pxl_col_in_fullres']
    if snakemake.params['only_within_tissue']:
        visium_allcoords = visium_allcoords[visium_allcoords['in_tissue']==1]

    mean_intensity = create_mean_intensity_table(msi_obj,
                                                 visium_allcoords.iloc[:,[4,3]],
                                                 sample,
                                                 verbose = snakemake.params['verbose'])
    create_mock_spaceranger_aggregated_intensity(mean_intensity,
                                        'output/'+sample+'/spaceranger_aggregated/',
                                        'output/'+sample+'/spaceranger/spatial/tissue_hires_image.png',
                                        visium_dir = "input/"+sample+"/visium/",
                                        verbose=snakemake.params['verbose'])
    
if __name__ == "__main__":
    main()