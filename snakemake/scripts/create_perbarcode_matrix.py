import pandas as pd
import scanpy
from sklearn.neighbors import KNeighborsClassifier
from scipy.spatial import distance_matrix
import os
import json
from collections import namedtuple
from scipy.spatial.distance import pdist
from scipy.sparse import csr_matrix
import shutil
from scipy.io import mmwrite
import gzip
import h5py
from pathlib import Path
import numpy as np

def estimate_interspot_distance(data):
#    data = visium_obj.obsm['spatial']
    classes = [1]*data.shape[0]
    knn = KNeighborsClassifier(n_neighbors=1,n_jobs=1,metric='euclidean')
    knn.fit(data, classes)
    dist, ind = knn.kneighbors(data,2)
    distances = [dist[i][1] for i in range(len(dist))]
    return(min(distances)/2)

def group_msi_perspot(msi_obj,visium_allcoords,visium_index):
    min_distance = estimate_interspot_distance(visium_allcoords)
    spot_distances = pd.DataFrame(distance_matrix(visium_allcoords,msi_obj.obsm['spatial']),index=visium_index,columns=msi_obj.obs.index)
    spot_distances['visium_spot']=spot_distances.index
    spot_distances_long = pd.wide_to_long(spot_distances, i="visium_spot",stubnames='MSI_',j='MSI_spot')
    spot_distances_long.columns = ['distance']
    close_points = spot_distances_long[spot_distances_long['distance'] < min_distance].reset_index()
    close_points['MSI_spot'] = ['MSI_'+str(spot_id) for spot_id in close_points['MSI_spot']]
    return(close_points)

def create_mean_intensity_table(msi_obj,visium_allcoords,visium_index):
    close_points = group_msi_perspot(msi_obj,visium_allcoords,visium_index)
    msi_data = msi_obj.X
    msi_data = pd.DataFrame.sparse.from_spmatrix(msi_data)
    msi_data.index = msi_obj.obs.index
    msi_data.columns = msi_obj.var['gene_ids']
    close_points.index = close_points['MSI_spot']
    intensity_matrix = close_points.merge(msi_data,left_index=True,right_index=True)
    def group_mean(group):
        return np.mean(group, axis=0)
    intensity_matrix_mean = intensity_matrix.drop(labels=['distance','MSI_spot'],axis=1).groupby(['visium_spot']).apply(lambda x: group_mean(x.values))
    return(intensity_matrix_mean)

def create_mock_spaceranger_mean_intensity(
    mean_intensity_table,
    output_folder_name="mock_spaceranger",
    msi_image_path="msi_he.png",
    visium_dir=None,
    verbose=True
):
    """
    Creates mock Space Ranger formatted output data for a sample.

    Args:
        input_dir (str): Path to the input directory.
        output_folder_name (str, optional): Name of the output folder. Defaults to "mock_spaceranger".
        msi_image_path (str, optional): Path to the MSI image file. Defaults to "msi_he.png".
        msi_coord_fname (str, optional): Name of the MSI coordinate file. Defaults to "msi_coords.csv".
        msi_coord_original_colnames (tuple, optional): Column names in the MSI coordinate file with original coordinates. Defaults to ("MSI_original_coordinate_x", "MSI_original_coordinate_y").
        msi_coord_new_colnames (tuple, optional): Column names in the MSI coordinate file with new coordinates. Defaults to ("MSI_new_coordinate_y", "MSI_new_coordinate_x").
        msi_spot_prefix (str, optional): Prefix for spot IDs in the MSI data. Defaults to "MSI".
        msi_feat_prefix (str, optional): Prefix for feature IDs in the MSI data. Defaults to "mz".
        visium_sf_json_path (str, optional): Path to the Visium scalefactor JSON file (optional).
        msi_peak_data_path (str, optional): Path to the MSI peak data file (optional).
        verbose (bool, optional): Print informational messages during execution. Defaults to True.

    Returns:
        None
    """

    # Create output folders
    output_path = os.path.join(output_folder_name)
    Path(output_folder_name).mkdir(parents=True, exist_ok=True)
#    os.makedirs(output_path,exist_ok=True)
    spatial_path = os.path.join(output_path, "spatial")
    Path(spatial_path).mkdir(parents=True, exist_ok=True)
#    os.makedirs(spatial_path, exist_ok=True)
    filtered_path = os.path.join(output_path, "filtered_feature_bc_matrix")
    Path(filtered_path).mkdir(parents=True, exist_ok=True)
#    os.makedirs(filtered_path, exist_ok=True)

    with open(visium_dir+"/spatial/scalefactors_json.json", "r") as f:
        st_json = json.load(f)
    scale_factor = st_json['tissue_hires_scalef']
    msi_json = st_json

    # Define namedtuple for tissue positions
    TissuePosition = namedtuple(
        "TissuePosition",
        ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"],
    )

    # Read MSI coordinates
    if verbose:
        print("Reading MSI coordinate file...")
    
    visium_coords = pd.read_csv(visium_dir+"/spatial/tissue_positions_list.csv",header=None,index_col=0)

    msi_tissue_pos = visium_coords.loc[list(mean_intensity_table.index),:]
    # Write tissue_positions.csv
    msi_tissue_pos.to_csv(os.path.join(spatial_path, "tissue_positions_list.csv"),index=True,header=False)
    if verbose:
        print(f"The new MSI coordinate file has been saved, containing {len(msi_tissue_pos)} spots/pixels.")

    # Write new scale factor JSON (if data available)
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
    msi_peaks = mean_intensity_table
   
    # Create feature IDs
    msi_peaks_features = pd.DataFrame(
        {"id1": msi_peaks.columns}
    )
    msi_peaks_features["id2"] = msi_peaks_features["id1"]

    # Create barcode IDs
    msi_peaks_barcodes = msi_peaks.index.astype(str)

    #msi_peaks_barcodes = msi_spot_prefix + "_" + msi_peaks["Unnamed: 0"].astype(str)

    # Prepare MSI peak data matrix
    if verbose:
        print("Preparing MSI peak data...")
    
    msi_peaks_matrix = msi_peaks.transpose().copy()  # Select and transpose data
    msi_peaks_matrix.index = msi_peaks_features["id1"].tolist()  # Set row names
    msi_peaks_matrix.columns = msi_peaks_barcodes.tolist()  # Set column names

    # Convert to sparse matrix (optional)
    # You might need to install `scikit-learn` for sparse matrix functionality
#    from sklearn.preprocessing import csr_matrix
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

    # Gzip compression can be done using external libraries like `shutil`
#    from shutil import make_archive
#    make_archive(os.path.join(filtered_path, "matrix.mtx"),"gz")

    if verbose:
        print("Spaceranger 'filtered_feature_bc_matrix' folder content ready!")

    if verbose:
        print("Writing filtered_feature_bc_matrix.h5 file...")

    from scipy.sparse import coo_array
    msi_peaks_mtx_csr = msi_peaks_mtx.tocoo()
    msi_peaks_mtx_csr = coo_array((msi_peaks_mtx_csr.data, (msi_peaks_mtx_csr.col, msi_peaks_mtx_csr.row)), shape=(msi_peaks_mtx_csr.shape[1],msi_peaks_mtx_csr.shape[0])).tocsr()
    hf = h5py.File(os.path.join(output_path, "filtered_feature_bc_matrix.h5"), 'w')
    g1 = hf.create_group('matrix')
    g1.create_dataset('library_ids', data="utf-8")
    g1.create_dataset('data', data=msi_peaks_mtx_csr.data)
    g1.create_dataset('indices', data=msi_peaks_mtx_csr.indices)
    g1.create_dataset('indptr', data=msi_peaks_mtx_csr.indptr)
    g1.create_dataset('shape', data=msi_peaks_mtx.shape)
    g1.create_dataset('barcodes', data=msi_peaks_barcodes.tolist())
    g2 = g1.create_group('features')
    g2.create_dataset('_all_tags_keys',data='genome')
    g2.create_dataset('id', data=msi_peaks_features["id1"].tolist())
    g2.create_dataset('name', data=msi_peaks_features["id1"].tolist())
    g2.create_dataset('feature_type', data=len(msi_peaks_features["id1"].tolist()) * ['Gene Expression'])
    g2.create_dataset('genome', data=len(msi_peaks_features["id1"].tolist()) * ['unknown'])
    hf.close()

    if verbose:
        print("Spaceranger 'filtered_feature_bc_matrix.h5' file ready!")


    if verbose:
        print("Mock Space Ranger output completed for sample at path "+output_path)

def main():
    sample = snakemake.params['sample']
    msi_obj = scanpy.read_visium('output/'+sample+'/spaceranger/',library_id='myLib')
    msi_obj.obs['filler'] = 1
    visium_obj = scanpy.read_visium("input/"+sample+"/visium/",library_id='myLib')
    visium_allcoords = pd.read_csv("input/"+sample+'/visium/spatial/tissue_positions_list.csv',header=None,index_col=0).iloc[:,[4,3]]
    mean_intensity = create_mean_intensity_table(msi_obj,visium_allcoords,visium_allcoords.index)
    create_mock_spaceranger_mean_intensity(mean_intensity,
                                        'output/'+sample+'/spaceranger_meanIntensity/',
                                        'output/'+sample+'/spaceranger/spatial/tissue_hires_image.png',
                                            visium_dir = "input/"+sample+"/visium/")
    
if __name__ == "__main__":
    main()