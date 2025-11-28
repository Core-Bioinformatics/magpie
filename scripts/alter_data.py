"""
Step 1 in MAGPIE pipeline: Mapping coordinates in MSI data to Visium.

This script:
1. Reads user-specified landmark pairs for MSI-> Visium or MSI -> MSI H&E and H&E -> Visium alignment.
2. Applies either affine or thin-plate spline (TPS) transforms at each stage.
3. Produces transformed MSI coordinates and optionally transforms MSI H&E images.
4. Saves diagnostic plots to verify alignment quality.

The pipeline is designed for use within Snakemake and expects the following folder structure:
input/<sample>/msi/
input/<sample>/visium/
output/<sample>/
"""

import pandas as pd
import matplotlib.pyplot as plt
from skimage.io import imread, imsave
from skimage.color import rgb2gray
import os.path
import sys
import glob
from skimage.transform import (  
    AffineTransform,
    ThinPlateSplineTransform,
    matrix_transform,
    warp)


# Identify transform needed to map from MSI dim reduction to Visium H&E and apply it
def map_coords_noHE(sample,
                    transform,
                    verbose=True):
    
    """
    Map MSI coordinates directly to Visium H&E space *without* an MSI H&E image.

    Parameters
    ----------
    sample : str
        Sample ID used to locate input files under input/<sample>/.
    transform : {'affine', 'TPS'}
        The geometric transformation to apply based on landmark pairs.
    verbose : bool
        If True, print status messages during processing.

    Returns
    -------
    dict
        {
            'transformed_coords': pd.DataFrame with transformed MSI coordinates,
            'msi_he_image': None (returned for consistency with other options)
        }

    Notes
    -----
    This branch is used when no MSI H&E image is available. Only coordinate
    mapping is performed; no image warping is done.
    """

    if verbose:
        print("Mapping MSI data to Visium H&E")
    # read landmark and MSI coordinate files
    landmarks = pd.read_csv('input/'+sample+'/landmarks_noHE.csv')
    # if transformed coordinates were saved by shiny app use those instead, otherwise use standard
    if os.path.isfile('input/'+sample+'/msi/MSI_metadata_modified.csv'):
        msi_coords = pd.read_csv('input/'+sample+'/msi/MSI_metadata_modified.csv')
    else:
        msi_coords = pd.read_csv('input/'+sample+'/msi/MSI_metadata.csv')

    # apply affine or TPS transform to coordinates
    if (transform=='affine'):
        if verbose:
            print("Using affine transform to map MSI data to Visium H&E")
        tfm = AffineTransform()
        tfm.estimate(landmarks.iloc[:,:2],landmarks.iloc[:,2:4])
        msi_coords_tfm = pd.DataFrame(matrix_transform(
            msi_coords[['x', 'y']],
            tfm.params))
    elif (transform=='TPS'):
        if verbose:
            print("Using TPS transform to map MSI data to Visium H&E")
        tfm = ThinPlateSplineTransform()
        tfm.estimate((landmarks.iloc[:,:2]).to_numpy(),(landmarks.iloc[:,2:4]).to_numpy())
        msi_coords_tfm = pd.DataFrame(tfm(msi_coords[['x','y']].to_numpy()))

    msi_coords_tfm.index = msi_coords.index
    return {'transformed_coords':msi_coords_tfm,'msi_he_image':None}

# Identify transform needed to map from MSI dim reduction to MSI H&E and apply it
def map_coords_MSI2HE(sample,
                      transform,
                      msi_he_img,
                      verbose=True):
    
    """
    Map MSI pixel coordinates onto their corresponding MSI H&E image.

    Parameters
    ----------
    sample : str
        Sample ID used to locate MSI metadata and landmarks.
    transform : {'affine', 'TPS'}
        Transformation type inferred from MSI→MSI H&E landmark pairs.
    msi_he_img : str or None
        Path to MSI H&E image. A modified version is used if available.
    verbose : bool
        If True, print detailed progress messages.

    Returns
    -------
    pd.DataFrame
        Transformed MSI coordinates in MSI H&E space.

    Notes
    -----
    This step forms the first stage of two-step MSI -> H&E -> Visium coregistration
    when MSI H&E is available. A diagnostic plot with overlaid transformed
    coordinates is saved to allow quality inspection.
    """

    if verbose:
        print("Mapping MSI data to MSI H&E")

    # read landmark and MSI coordinate files
    landmarks = pd.read_csv('input/'+sample+'/landmarks_MSI2HE.csv')
    # if transformed coordinates were saved by shiny app use those instead, otherwise use standard
    if os.path.isfile('input/'+sample+'/msi/MSI_metadata_modified.csv'):
        msi_coords = pd.read_csv('input/'+sample+'/msi/MSI_metadata_modified.csv')
    else:
        msi_coords = pd.read_csv('input/'+sample+'/msi/MSI_metadata.csv')

    # If dimensionality reduction was saved by shiny app this can be used for visualisation
    if os.path.isfile('input/'+sample+'/msi/MSI_dimreduction.csv'):
        msi_dimred = pd.read_csv('input/'+sample+'/msi/MSI_dimreduction.csv')

    # apply affine or TPS transform to coordinates
    if (transform=='affine'):
        if verbose:
            print("Using affine transform to map MSI data to MSI H&E")
        tfm = AffineTransform()
        tfm.estimate(landmarks.iloc[:,:2],landmarks.iloc[:,2:4])
        msi_coords_tfm = pd.DataFrame(matrix_transform(
            msi_coords[['x', 'y']],
            tfm.params))
    elif (transform=='TPS'):
        if verbose:
            print("Using TPS transform to map MSI data to MSI H&E")
        tfm = ThinPlateSplineTransform()
        tfm.estimate((landmarks.iloc[:,:2]).to_numpy(),(landmarks.iloc[:,2:4]).to_numpy())
        msi_coords_tfm = pd.DataFrame(tfm(msi_coords[['x','y']].to_numpy()))

    msi_coords_tfm.index = msi_coords.index

    # plot coordinates over MSI H&E to test success of first stage of coregistration
    if verbose:
        print("Saving original MSI H&E with transformed MSI coordinates overlaid...")
    fig, ax = plt.subplots(nrows=1, ncols=1 )
    
    if os.path.isfile('input/'+sample+'/msi/MSI_HE_modified.jpg'):
        out_image = imread('input/'+sample+'/msi/MSI_HE_modified.jpg')
    else:
        out_image = imread(msi_he_img)

    transformed_coords = msi_coords_tfm.copy()
    transformed_coords.columns = ['x','y']
    transformed_coords['spot_id']=msi_coords['spot_id']
    # plot the dimensionality reduction if available, otherwise use standard colouring
    if os.path.isfile('input/'+sample+'/msi/MSI_dimreduction.csv'):
        transformed_coords = pd.merge(transformed_coords,msi_dimred[['spot_id','color']],on='spot_id')
    else :
        transformed_coords['color']=1
    # plot coordinates on top of H&E image
    plt.imshow(out_image)
    plt.scatter(x=transformed_coords['x'],y=transformed_coords['y'],c=transformed_coords['color'],s=0.1,alpha=0.5)
    fig.savefig('output/' + sample +'/MSI_HE_withMSI2HECoords.png')
    plt.close(fig)

    return msi_coords_tfm

# Identify transform needed to map from MSI H&E to Visium H&E and apply it
def map_coords_HE2HE(sample,
                     msi_coords,
                     transform,
                     msi_he_img,
                     verbose=True):
    
    """
    Map the MSI H&E image and MSI coordinates to the Visium H&E coordinate system.

    Parameters
    ----------
    sample : str
        Sample ID for locating landmark files and Visium images.
    msi_coords : pd.DataFrame
        MSI coordinates from the previous transformation stage.
    transform : {'affine', 'TPS'}
        Transformation type using MSI H&E → Visium H&E landmark pairs.
    msi_he_img : str
        Path to the MSI H&E image to be warped.
    verbose : bool
        If True, print status messages.

    Returns
    -------
    dict
        {
            'transformed_coords': pd.DataFrame (MSI coords in Visium space),
            'msi_he_image': 2D/3D ndarray (warped MSI H&E aligned to Visium)
        }

    Notes
    -----
    Output dimensions match the Visium H&E resolution.
    """

    if verbose:
        print("Mapping MSI H&E to Visium H&E")

    # read landmark and MSI coordinate files as well as Visium H&E (for shape to transform MSI H&E)
    landmarks = pd.read_csv('input/'+sample+'/landmarks_HE2HE.csv')
    visium_he_img = imread('input/'+sample+'/visium/spatial/tissue_hires_image.png')
    # if transformed H&E was saved by shiny app use that instead, otherwise use standard
    if os.path.isfile('input/'+sample+'/msi/MSI_HE_modified.jpg'):
        msi_he_img = imread('input/'+sample+'/msi/MSI_HE_modified.jpg')
    else:
        msi_he_img = imread(msi_he_img)

    rows, cols = visium_he_img.shape[:2]

    # apply affine or TPS transform to coordinates and MSI H&E
    if (transform=='affine'):
        if verbose:
            print("Using affine transform to map MSI H&E to Visium H&E")

        tfm = AffineTransform()
        tfm.estimate(landmarks.iloc[:,:2],landmarks.iloc[:,2:4])
        msi_coords_tfm = pd.DataFrame(matrix_transform(
            msi_coords[[0,1]],
            tfm.params))
        tfm = AffineTransform()
        tfm.estimate(landmarks.iloc[:,2:4],landmarks.iloc[:,:2])
    elif (transform=='TPS'):
        if verbose:
            print("Using TPS transform to map MSI H&E to Visium H&E")

        tfm = ThinPlateSplineTransform()
        tfm.estimate((landmarks.iloc[:,:2]).to_numpy(),(landmarks.iloc[:,2:4]).to_numpy())
        msi_coords_tfm = pd.DataFrame(tfm(msi_coords[[0,1]].to_numpy()))
        tfm = ThinPlateSplineTransform()
        tfm.estimate((landmarks.iloc[:,2:4]).to_numpy(),(landmarks.iloc[:,:2]).to_numpy())

    transformed_image = warp(msi_he_img, tfm,output_shape=(rows,cols))
    msi_coords_tfm.index = msi_coords.index

    # return both the transformed coordinates and transformed H&E image
    return {'transformed_coords':msi_coords_tfm,
            'msi_he_image':transformed_image}

# check whether there is an MSI H&E image and use the transformation pipeline depending on result
def apply_mapping(sample,
                  msi_he_img,
                  verbose=True):
    
    """
    Dispatch function selecting the correct mapping pipeline depending on
    whether an MSI H&E image is present.

    If MSI H&E exists:
        1. MSI → MSI H&E transform
        2. MSI H&E → Visium H&E transform
    Else:
        MSI → Visium H&E transform only

    Returns a dict with transformed coordinates and (if available) warped MSI H&E.
    """

    if not (msi_he_img is None):
        if verbose:
            print("MSI H&E image identified.")
        intermediate_coords = map_coords_MSI2HE(sample,snakemake.params['MSI2HE_transform'],msi_he_img,verbose=verbose)
        return(map_coords_HE2HE(sample,intermediate_coords,snakemake.params['HE2HE_transform'],msi_he_img,verbose=verbose))
    else:
        if verbose:
            print("No MSI H&E image identified.")
        return(map_coords_noHE(sample,snakemake.params['no_HE_transform'],verbose=verbose))

# run full coregistration pipeline on selected sample
def run_coreg(sample,
              verbose=True):

    """
    Run the full MSI -> Visium coregistration pipeline on a given sample.

    Steps
    -----
    1. Detect MSI H&E image (jpg/png/tiff) if available.
    2. Apply appropriate mapping strategy (with or without MSI H&E).
    3. Save:
        - transformed coordinates (output/<sample>/transformed.csv)
        - transformed image (either MSI H&E or Visium H&E)
        - diagnostic images with overlaid coordinates

    Parameters
    ----------
    sample : str
        Sample name passed via Snakemake.
    verbose : bool
        If True, print progress and file-saving messages.

    Notes
    -----
    The function expects Snakemake to provide parameters for:
        snakemake.params['MSI2HE_transform']
        snakemake.params['HE2HE_transform']
        snakemake.params['no_HE_transform']
    """

    # check if there is any image file MSI_HE.jpg/tiff/png
    if glob.glob('input/'+sample+'/msi/MSI_HE.*') != []:
        msi_he_img = glob.glob('input/'+sample+'/msi/MSI_HE.*')
        msi_he_img = [x for x in msi_he_img if any(ext in x for ext in ['tiff','png','jpg'])]
        if msi_he_img == []:
            msi_he_img = None
        else :
            msi_he_img = msi_he_img[0]
    else :
        msi_he_img = None

    # get new MSI coordinates and image
    transformed_result = apply_mapping(sample,
                                       msi_he_img,
                                       verbose=verbose)
    transformed_coords = transformed_result['transformed_coords']
    msi_he_image = transformed_result['msi_he_image']
    transformed_coords.columns = ['x','y']

    # identify if output image is MSI H&E (if available) or Visium H&E
    
    if not (msi_he_img is None):
        out_image = msi_he_image
        # overlay transformed coordinates over MSI H&E
        if verbose:
            print("Saving coregistered MSI image with overlaid new MSI coordinates")
        fig, ax = plt.subplots(nrows=1, ncols=1 )
        plt.imshow(out_image)
        plt.scatter(x=transformed_coords['x'],y=transformed_coords['y'],s=0.1,c='r',alpha=0.7)
        fig.savefig('output/'+sample+'/transformed_withCoords.png')
        plt.close(fig)
    else:
        visium_he_img = imread('input/'+sample+'/visium/spatial/tissue_hires_image.png')
        out_image = visium_he_img

    # save H&E image (either MSI H&E if available or Visium otherwise)
    plt.imsave(arr=out_image,fname='output/'+sample+'/transformed.png')

    # save Visium H&E image with transformed coordinates overlaid
    if verbose:
        print("Saving Visium H&E with new MSI coordinates overlaid...")
    visium_he_img = imread('input/'+sample+'/visium/spatial/tissue_hires_image.png')
    fig, ax = plt.subplots(nrows=1, ncols=1)
    plt.imshow(visium_he_img)
    plt.scatter(x=transformed_coords['x'],y=transformed_coords['y'],s=0.1,c='r',alpha=0.7)
    fig.savefig('output/'+sample+'/transformed_withCoords_VisiumHE.png')

    # save transformed coordinates
    if verbose:
        print("Saving new MSI coordinates...")
    if os.path.isfile('input/'+sample+'/msi/MSI_metadata_modified.csv'):
        msi_coords = pd.read_csv('input/'+sample+'/msi/MSI_metadata_modified.csv')
    else:
        msi_coords = pd.read_csv('input/'+sample+'/msi/MSI_metadata.csv')
    transformed_coords['spot_id']=msi_coords['spot_id']
    transformed_coords = transformed_coords[['spot_id','x','y']]
    transformed_coords.to_csv('output/'+sample+'/transformed.csv',index=False)

def main():
    # run on current sample
    run_coreg(snakemake.params['sample'],verbose=snakemake.params['verbose'])

if __name__ == "__main__":
    main()

