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

    return msi_coords_tfm

# Identify transform needed to map from MSI H&E to Visium H&E and apply it
def map_coords_HE2HE(sample,
                     msi_coords,
                     transform,
                     msi_he_img,
                     verbose=True):

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
    ### WRITE FULL DOCUMENTATION HERE!!!
    Also add verbose mode options
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
    else:
        visium_he_img = imread('input/'+sample+'/visium/spatial/tissue_hires_image.png')
        out_image = visium_he_img

    # save H&E image (either MSI H&E if available or Visium otherwise)
    plt.imsave(arr=out_image,fname='output/'+sample+'/transformed.png')

    # save Visium H&E image with transformed coordinates overlaid
    if verbose:
        print("Saving Visium H&E with new MSI coordinates overlaid...")
    visium_he = imread('input/'+sample+'/visium/spatial/tissue_hires_image.png')
    fig, ax = plt.subplots(nrows=1, ncols=1)
    plt.imshow(visium_he)
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

