import pandas as pd
import matplotlib.pyplot as plt
from skimage.io import imread, imsave
from skimage.color import rgb2gray
from TPS import *
import os.path
import sys

# Snakemake automatically passes input and output files via 'snakemake.input' and 'snakemake.output'
#input_file = snakemake.input[0]

def map_coords_noHE(file,transform):
    from skimage.transform import (  # pylint: disable=no-name-in-module
        AffineTransform,
        matrix_transform)
    """Alter the DataFrame as needed"""
    # Example alteration: Add a new column with doubled values of 'Some_Column'
    landmarks = pd.read_csv(file+'/landmarks_noHE.csv',index_col=0)
    msi_coords = pd.read_csv(file+'/msi/MSI_metadata.csv')
    if (transform=='affine'):
        tfm = AffineTransform()
        tfm.estimate(landmarks.iloc[:,:2],landmarks.iloc[:,2:4])
        msi_coords_tfm = pd.DataFrame(matrix_transform(
            msi_coords[['x', 'y']],
            tfm.params))
    elif (transform=='TPS'):
        tfm = TpsTransform()
        tfm.estimate((landmarks.iloc[:,:2]).to_numpy(),(landmarks.iloc[:,2:4]).to_numpy())
        msi_coords_tfm = pd.DataFrame(tfm(msi_coords[['x','y']].to_numpy()))

    return {'transformed_coords':msi_coords_tfm,'msi_he_image':None}

def map_coords_MSI2HE(file,transform):
    from skimage.transform import (  # pylint: disable=no-name-in-module
        AffineTransform,
        matrix_transform,
        warp)
    """Alter the DataFrame as needed"""
    # Example alteration: Add a new column with doubled values of 'Some_Column'

    landmarks = pd.read_csv(file+'/landmarks_MSI2HE.csv',index_col=0)
    msi_coords = pd.read_csv(file+'/msi/MSI_metadata.csv')
    if (transform=='affine'):
        tfm = AffineTransform()
        tfm.estimate(landmarks.iloc[:,:2],landmarks.iloc[:,2:4])
        msi_coords_tfm = pd.DataFrame(matrix_transform(
            msi_coords[['x', 'y']],
            tfm.params))
    elif (transform=='TPS'):
        tfm = TpsTransform()
        tfm.estimate((landmarks.iloc[:,:2]).to_numpy(),(landmarks.iloc[:,2:4]).to_numpy())
        msi_coords_tfm = pd.DataFrame(tfm(msi_coords[['x','y']].to_numpy()))

    return msi_coords_tfm

def map_coords_HE2HE(file,msi_coords,transform):
    from skimage.transform import (  # pylint: disable=no-name-in-module
        AffineTransform,
        matrix_transform,
        warp)
    """Alter the DataFrame as needed"""
    # Example alteration: Add a new column with doubled values of 'Some_Column'

    visium_he_img = imread(file+'/visium/spatial/tissue_hires_image.png')
    msi_he_img = imread(file+'/msi/MSI_HE.jpg')
    landmarks = pd.read_csv(file+'/landmarks_HE2HE.csv',index_col=0)
    rows, cols = visium_he_img.shape[:2]
    if (transform=='affine'):
        tfm = AffineTransform()
        tfm.estimate(landmarks.iloc[:,:2],landmarks.iloc[:,2:4])
        msi_coords_tfm = pd.DataFrame(matrix_transform(
            msi_coords[['x', 'y']],
            tfm.params))
        tfm = AffineTransform()
        tfm.estimate(landmarks.iloc[:,2:4],landmarks.iloc[:,:2])

    elif (transform=='TPS'):
        tfm = TpsTransform()
        tfm.estimate((landmarks.iloc[:,:2]).to_numpy(),(landmarks.iloc[:,2:4]).to_numpy())
        msi_coords_tfm = pd.DataFrame(tfm(msi_coords[[0,1]].to_numpy()))
        tfm = TpsTransform()
        tfm.estimate((landmarks.iloc[:,2:4]).to_numpy(),(landmarks.iloc[:,:2]).to_numpy())

    transformed_image = warp(msi_he_img, tfm,output_shape=(rows,cols))
    return {'transformed_coords':msi_coords_tfm,'msi_he_image':transformed_image}

def apply_mapping(file):
    if os.path.isfile(file+'/msi/MSI_HE.jpg'):
        intermediate_coords = map_coords_MSI2HE(file,'affine')
        return(map_coords_HE2HE(file,intermediate_coords,snakemake.params['transform_type']))
    else:
        return(map_coords_noHE(file,snakemake.params['transform_type']))
    
def run_coreg(sample):
    transformed_result = apply_mapping('input/'+sample)
    transformed_coords = transformed_result['transformed_coords']
    msi_he_image = transformed_result['msi_he_image']
    visium_he_img = imread('input/'+sample+'/visium/spatial/tissue_hires_image.png')
    fig, ax = plt.subplots(nrows=1, ncols=1 )
    if os.path.isfile('input/'+sample+'/msi/MSI_HE.jpg'):
        out_image = msi_he_image
    else:
        out_image = visium_he_img
    plt.imsave(arr=out_image,fname='output/'+sample+'/transformed.png')
    plt.imshow(out_image)
    plt.scatter(x=transformed_coords[0],y=transformed_coords[1],s=0.1,c='r',alpha=1)

    fig.savefig('output/'+sample+'/transformed_withCoords.png')
    transformed_coords.columns = ['x','y']
    # Step 3: Write modified data to new CSV
    transformed_coords.to_csv('output/'+sample+'/transformed.csv')

def main():
    run_coreg(snakemake.params['sample'])

if __name__ == "__main__":
    main()

