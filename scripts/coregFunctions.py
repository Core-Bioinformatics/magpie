import matplotlib.pyplot as plt
import numpy as np
import cv2
import gui
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import scanpy
from pathlib import Path
from argparse import Namespace, ArgumentParser, ArgumentDefaultsHelpFormatter
from pandas import read_csv, DataFrame  # type: ignore
from numpy import ones
from skimage.io import imread, imsave  # type: ignore


def prep_msi_image(output_directory,msi_intensities,msi_coords,setting,peaks=None,size=40):
    fig, ax = plt.subplots(nrows=1, ncols=1 )  # create figure & 1 axis
    ax.margins(x=0,y=0)

    if setting == 'PCA_1':
        from sklearn.decomposition import PCA
        if peaks != None:
            msi_intensities = msi_intensities[peaks]
        pca = PCA(n_components=1)
        reduction = pca.fit_transform(msi_intensities)
        ax.scatter(x=msi_coords['x'], y=msi_coords['y'], c=reduction,marker='.',s=100)

    if setting == 'PCA_3':
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import MinMaxScaler
        if peaks != None:
            msi_intensities = msi_intensities[peaks]
        pca = PCA(n_components=3)
        reduction = pd.DataFrame(pca.fit_transform(msi_intensities))
        scaler = MinMaxScaler()
        reduction_scaled = pd.DataFrame(scaler.fit_transform(reduction), columns=reduction.columns)
        reduction_colours = reduction_scaled.values.tolist()
        ax.scatter(x=msi_coords['x'], y=msi_coords['y'], c=reduction_colours,marker='s',s=200)
    
    elif setting == 'OnePeak':
        ax.scatter(x=msi_coords['x'], y=msi_coords['y'], c=msi_intensities[peaks[0]],marker='.')
        
    ax.set_axis_off()
    print(10*msi_coords['y'].max()/msi_coords['x'].max())
    fig.set_size_inches(10,10*(msi_coords['y'].max()/msi_coords['x'].max()))
    fig.savefig(output_directory + '/MSI_image.png', bbox_inches='tight', pad_inches=0)   # save the figure to file
    plt.close(fig)    # close the figure window


def get_landmarks(visium_img, msi_img, num_landmarks: int):
    """Compare the two images.
    
    Parameters
    ----------
    visium_img
    msi_img
    num_landmarks : int
        The number of landmarks to consider.
    
    Returns
    -------
    coord_df : DataFrame
        The landmark coordinates between both images.
    
    """
    
    visium_coords = np.ones((num_landmarks, 2))
    msi_coords = np.ones((num_landmarks, 2))
    fig, axes = plt.subplots(1, 2)
    axes[0].imshow(visium_img)
    axes[1].imshow(msi_img, cmap="inferno")
    
    for idx in range(num_landmarks):
        for coords, axis in zip([visium_coords, msi_coords], axes.ravel()):
            axis.set_title("Click here")

            point = fig.ginput(1, timeout=0, mouse_add=None, mouse_stop=None)[0]
            coords[idx, :2] = point

            axis.plot(point[0], point[1],color='red')
            axis.annotate(idx, xy = (point[0], point[1]), color='red',textcoords = 'data', fontsize=15)

            plt.draw()

            axis.set_title("")
    
  
    for axis in axes.ravel():
        axis.cla()
        axis.remove()
    plt.close("all")
    
    coord_df = DataFrame()
    coord_df[["visium_col", "visium_row"]] = visium_coords
    coord_df[["msi_col", "msi_row"]] = msi_coords
    return coord_df

import json
from skimage.transform import (  # pylint: disable=no-name-in-module
AffineTransform,
matrix_transform,
ProjectiveTransform,
EuclideanTransform,
SimilarityTransform,
warp)
def do_full_mapping_frominputs(output_directory: Path, 
                    visium_path: Path, 
                    num_landmarks: int,
                    msi_coords: DataFrame, 
                    msi_intensities: DataFrame, 
                    msi_option: str,
                    transform,
                    peak_list = None,
                    msi_he_image: Path = None):

    # Prepare Visium data
    visium_obj = scanpy.read_visium(visium_path)
    visium_he_img = imread(visium_path + '/spatial/tissue_hires_image.png')
    visium_coords = pd.DataFrame(visium_obj.obsm['spatial'])
    f = open(visium_path + 'spatial/scalefactors_json.json')
    scalef = json.load(f)
    f.close()
    scale_factor = scalef['tissue_hires_scalef']

    visium_coords['x']=visium_coords[[0]]*scale_factor
    visium_coords['y']=visium_coords[[1]]*scale_factor
    visium_obj.obsm['spatial_matching']=visium_obj.obsm['spatial']*scale_factor

    msi_coords['x'] = msi_coords['x'].transform(lambda x: (x - x.min()))
    msi_coords['y'] = msi_coords['y'].transform(lambda x: (x - x.min()))
    prep_msi_image(output_directory, msi_intensities, msi_coords, msi_option, peaks = peak_list)
    msi_coords['y'] = msi_coords['y'].transform(lambda x: (x.max() - x))

    msi_image = imread(output_directory + '/MSI_image.png')
    print(round(msi_coords['x'].max()))
    print(round(msi_coords['x'].max()))
    
    ratio = round(msi_coords['y'].max())/round(msi_coords['x'].max())
    print(round(1000*ratio))
    msi_image_scaled = cv2.resize(msi_image,(1000, round(1000*ratio)))
    msi_coords['x'] = msi_coords['x'].transform(lambda x: ((1000/x.max())*x))
    msi_coords['y'] = msi_coords['y'].transform(lambda x: ((1000*ratio/x.max())*x))
    plt.imshow(msi_image_scaled)
    
    # If we have MSI H&E then do that alignment first
    if msi_he_image != None:
            msi_he_img = imread(msi_he_image)
            landmarks_msi_he_msi = get_landmarks(msi_he_img, msi_image_scaled, num_landmarks)
            print(landmarks_msi_he_msi.iloc[:,2:4])
 #           tfm_msi = ProjectiveTransform()
            #tfm_msi = SimilarityTransform()
            tfm_msi = transform()
            tfm_msi.estimate(landmarks_msi_he_msi.iloc[:,2:4],landmarks_msi_he_msi.iloc[:,:2])
        
            new_coords_msi = pd.DataFrame(matrix_transform(
                        msi_coords[['x', 'y']],
                        tfm_msi.params))
            msi_image = msi_he_img
            temp_table = pd.DataFrame({'spot_id':msi_coords['spot_id'],
                                  'MSI_original_coordinate_x':msi_coords['x'],
                                  'MSI_original_coordinate_y':msi_coords['y'],
                                  'MSI_new_coordinate_x':new_coords_msi[0],
                                  'MSI_new_coordinate_y':new_coords_msi[1]})
        
    else:
            msi_image = msi_image_scaled
            temp_table = pd.DataFrame({'spot_id':msi_coords['spot_id'],
                                  'MSI_original_coordinate_x':msi_coords['x'],
                                  'MSI_original_coordinate_y':msi_coords['y'],
                                  'MSI_new_coordinate_x':msi_coords['x'],
                                  'MSI_new_coordinate_y':msi_coords['y']})
            new_coords_msi = pd.DataFrame({0:msi_coords['x'],1:msi_coords['y']})
            landmarks_msi_he_msi = 0
        

    landmarks_msi_he_visium_he = get_landmarks(visium_he_img, msi_image, num_landmarks)
    tfm_visium = transform()
    tfm_visium.estimate(landmarks_msi_he_visium_he.iloc[:,2:4],landmarks_msi_he_visium_he.iloc[:,:2])
    new_new_coords_msi = pd.DataFrame(matrix_transform(
                new_coords_msi,
                tfm_visium.params))
    new_new_coords_msi[0] = new_new_coords_msi[0]/scale_factor
    new_new_coords_msi[1] = new_new_coords_msi[1]/scale_factor
    msi_table = pd.DataFrame({'spot_id':msi_coords['spot_id'],
                          'MSI_original_coordinate_x':msi_coords['x'],
                          'MSI_original_coordinate_y':msi_coords['y'],
                          'MSI_middle_coordinate_x':temp_table['MSI_new_coordinate_x'],
                          'MSI_middle_coordinate_y':temp_table['MSI_new_coordinate_y'],    
                          'MSI_new_coordinate_x':new_new_coords_msi[0],
                          'MSI_new_coordinate_y':new_new_coords_msi[1]})
    rows, cols = visium_he_img.shape[:2]
    tfm_visium = transform()
    tfm_visium.estimate(landmarks_msi_he_visium_he.iloc[:,:2],landmarks_msi_he_visium_he.iloc[:,2:4])
    if msi_he_image != None:
        out_he_image =warp(msi_he_img, tfm_visium,output_shape=(rows,cols))
    else: 
        out_he_image =warp(msi_image, tfm_visium,output_shape=(rows,cols))

    tps = TpsTransform()
    tps.estimate(landmarks_msi_he_visium_he.iloc[:,2:4].to_numpy(),landmarks_msi_he_visium_he.iloc[:,:2].to_numpy())
    trans_coord = pd.DataFrame(tps(new_coords_msi.to_numpy()))
    trans_coord[0] = trans_coord[0]/scale_factor
    trans_coord[1] = trans_coord[1]/scale_factor
    tps_coords = pd.DataFrame({'spot_id':msi_table['spot_id'],
                              'MSI_original_coordinate_x':msi_table['MSI_original_coordinate_x'],
                              'MSI_original_coordinate_y':msi_table['MSI_original_coordinate_x'],
                              'MSI_middle_coordinate_x':msi_table['MSI_middle_coordinate_x'],
                              'MSI_middle_coordinate_y':msi_table['MSI_middle_coordinate_y'],
                              'MSI_new_coordinate_x':trans_coord[0],
                              'MSI_new_coordinate_y':trans_coord[1]})
    tps = TpsTransform()
    tps.estimate(landmarks_msi_he_visium_he.iloc[:,:2].to_numpy(),landmarks_msi_he_visium_he.iloc[:,2:4].to_numpy())
    
    if msi_he_image != None:
        out_he_image_tps = warp(msi_he_img,tps,output_shape=(rows,cols))
    else: 
        out_he_image_tps = warp(msi_image,tps,output_shape=(rows,cols))

    return((msi_table,landmarks_msi_he_visium_he,out_he_image,visium_obj,tps_coords,out_he_image_tps,landmarks_msi_he_msi))

import os
import json
from collections import namedtuple
from scipy.spatial.distance import pdist
from scipy.sparse import csr_matrix
import pandas as pd
import shutil
from scipy.io import mmwrite
import gzip

def create_mock_spaceranger(
    input_dir,
    output_folder_name="mock_spaceranger",
    msi_image_path="msi_he.png",
    msi_coord_fname="msi_coords.csv",
    msi_coord_original_colnames=("MSI_original_coordinate_x", "MSI_original_coordinate_y"),
    msi_coord_new_colnames=("MSI_new_coordinate_y", "MSI_new_coordinate_x"),
    msi_spot_prefix="MSI",
    msi_feat_prefix="mz",
    visium_sf_json_path=None,
    msi_peak_data_path=None,
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

    if verbose:
        print(f"Preparing Space Ranger formatted output data for sample at {input_dir}")

    # Create output folders
    output_path = os.path.join(input_dir, output_folder_name)
    os.makedirs(output_path, exist_ok=True)
    spatial_path = os.path.join(output_path, "spatial")
    os.makedirs(spatial_path, exist_ok=True)
    filtered_path = os.path.join(output_path, "filtered_feature_bc_matrix")
    os.makedirs(filtered_path, exist_ok=True)

    # Define namedtuple for tissue positions
    TissuePosition = namedtuple(
        "TissuePosition",
        ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"],
    )

    # Read MSI coordinates
    if verbose:
        print("Reading MSI coordinate file...")
    msi_coords = pd.read_csv(os.path.join(input_dir, msi_coord_fname), index_col=0)

    # Add "in_tissue" column
    msi_coords["in_tissue"] = 1

    msi_tissue_pos = msi_coords.rename(
    columns={
        "spot_id":"barcode",
        msi_coord_original_colnames[0]: "array_row",
        msi_coord_original_colnames[1]: "array_col",
        msi_coord_new_colnames[0]: "pxl_row_in_fullres",
        msi_coord_new_colnames[1]: "pxl_col_in_fullres",
    }
)

    # Rename "barcodes" so that they are characters
    msi_tissue_pos["barcode"] = msi_spot_prefix + "_" + msi_tissue_pos["barcode"].astype(str)

    msi_tissue_pos = msi_tissue_pos[['barcode','in_tissue','array_row','array_col','pxl_row_in_fullres','pxl_col_in_fullres']]
    # Write tissue_positions.csv
    msi_tissue_pos.to_csv(os.path.join(spatial_path, "tissue_positions_list.csv"),index=False,header=False)
    if verbose:
        print(f"The new MSI coordinate file has been saved, containing {len(msi_tissue_pos)} spots/pixels.")

    # Load matching Visium scalefactor JSON file (if provided)
    if visium_sf_json_path:
        with open(visium_sf_json_path, "r") as f:
            st_json = json.load(f)
        msi_json = st_json

    # Update spot diameter based on MSI center-to-center distance
    if msi_peak_data_path:  # Only calculate diameter if peak data is provided
        # Assuming pxl_row_in_fullres and pxl_col_in_fullres are spatial coordinates
        distances = pdist(msi_tissue_pos[["pxl_row_in_fullres", "pxl_col_in_fullres"]])
        c_c_dist = distances[distances > 0].min()
        msi_json["spot_diameter_fullres"] = c_c_dist
        msi_json["fiducial_diameter_fullres"] = 0

    # Write new scale factor JSON (if data available)
    if (msi_peak_data_path and visium_sf_json_path) or msi_json:
        if verbose:
            print("Creating JSON file...")
        with open(os.path.join(spatial_path, "scalefactors_json.json"), "w") as f:
            json.dump(msi_json, f, indent=4)

    # Copy image into spatial folder
    if verbose:
        print("Copying tissue image...")
    shutil.copyfile(os.path.join(input_dir, msi_image_path), os.path.join(spatial_path, "tissue_hires_image.png"))
    shutil.copyfile(os.path.join(input_dir, msi_image_path), os.path.join(spatial_path, "tissue_lowres_image.png"))

    if verbose:
        print("Spaceranger 'spatial' folder content ready!")

    # Read MSI peak data file
    if verbose:
        print("Reading MSI peak data file...")
    msi_peaks = pd.read_csv(msi_peak_data_path, header=0)
    #print(msi_peaks[:5])

    # Create feature IDs
    msi_peaks_features = pd.DataFrame(
        {"id1": [msi_feat_prefix + "-" + col for col in msi_peaks.columns[1:]]}
    )
    msi_peaks_features["id2"] = msi_peaks_features["id1"]

    # Create barcode IDs
    msi_peaks_barcodes = msi_spot_prefix + "_" + msi_peaks["Unnamed: 0"].astype(str)

    # Prepare MSI peak data matrix
    if verbose:
        print("Preparing MSI peak data...")
    
    msi_peaks_matrix = msi_peaks.iloc[:, 1:].transpose().copy()  # Select and transpose data
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

def estimate_interspot_distance(data):
    from sklearn.neighbors import KNeighborsClassifier
#    data = visium_obj.obsm['spatial']
    classes = [1]*data.shape[0]
    knn = KNeighborsClassifier(n_neighbors=1,n_jobs=1,metric='euclidean')
    knn.fit(data, classes)
    dist, ind = knn.kneighbors(data,2)
    distances = [dist[i][1] for i in range(len(dist))]
    return(min(distances)/2)

def group_msi_perspot(msi_obj,visium_allcoords,visium_index):
    from scipy.spatial import distance_matrix
    min_distance = estimate_interspot_distance(visium_allcoords)
    spot_distances = pd.DataFrame(distance_matrix(visium_allcoords,msi_obj.obsm['spatial']),index=visium_index,columns=msi_obj.obs.index)
    spot_distances['visium_spot']=spot_distances.index
    spot_distances_long = pd.wide_to_long(spot_distances, i="visium_spot",stubnames='MSI_',j='MSI_spot')
    spot_distances_long.columns = ['distance']
    close_points = spot_distances_long[spot_distances_long['distance'] < min_distance].reset_index()
    close_points['MSI_spot'] = ['MSI_'+str(spot_id) for spot_id in close_points['MSI_spot']]
    return(close_points)

def create_mean_intensity_table(msi_obj,visium_obj):
    close_points = group_msi_perspot(msi_obj,visium_obj)
    msi_data = msi_obj.X
    msi_data = pd.DataFrame.sparse.from_spmatrix(msi_data)
    msi_data.index = msi_obj.obs.index
    msi_data.columns = msi_obj.var['gene_ids']
    close_points.index = close_points['MSI_spot']
    intensity_matrix = close_points.merge(msi_data,left_index=True,right_index=True)
    intensity_matrix_mean = intensity_matrix.drop(labels=['distance','MSI_spot'],axis=1).groupby(['visium_spot']).mean()
    return(intensity_matrix_mean)