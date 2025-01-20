import pandas as pd
import matplotlib.pyplot as plt
from skimage.io import imread, imsave
from skimage.color import rgb2gray
import os.path
import sys
from skimage.transform import (  
    AffineTransform,
    ThinPlateSplineTransform,
    matrix_transform,
    warp)

# Identify transform needed to map from MSI dim reduction to Visium H&E and apply it
def map_coords_noHE(file,transform):

    # read landmark and MSI coordinate files
    landmarks = pd.read_csv(file+'/landmarks_noHE.csv')
    if os.path.isfile(file+'/msi/MSI_metadata_modified.csv'):
        msi_coords = pd.read_csv(file+'/msi/MSI_metadata_modified.csv')
    else:
        msi_coords = pd.read_csv(file+'/msi/MSI_metadata.csv')

    # apply affine or TPS transform to coordinates
    if (transform=='affine'):
        tfm = AffineTransform()
        tfm.estimate(landmarks.iloc[:,:2],landmarks.iloc[:,2:4])
        msi_coords_tfm = pd.DataFrame(matrix_transform(
            msi_coords[['x', 'y']],
            tfm.params))
    elif (transform=='TPS'):
        tfm = ThinPlateSplineTransform()
        tfm.estimate((landmarks.iloc[:,:2]).to_numpy(),(landmarks.iloc[:,2:4]).to_numpy())
        msi_coords_tfm = pd.DataFrame(tfm(msi_coords[['x','y']].to_numpy()))
    msi_coords_tfm.index = msi_coords.index
    return {'transformed_coords':msi_coords_tfm,'msi_he_image':None}

# Identify transform needed to map from MSI dim reduction to MSI H&E and apply it
def map_coords_MSI2HE(file,transform):

    # read landmark and MSI coordinate files
    landmarks = pd.read_csv(file+'/landmarks_MSI2HE.csv')

    if os.path.isfile(file+'/msi/MSI_metadata_modified.csv'):
        msi_coords = pd.read_csv(file+'/msi/MSI_metadata_modified.csv')
    else:
        msi_coords = pd.read_csv(file+'/msi/MSI_metadata.csv')


    # apply affine or TPS transform to coordinates
    if (transform=='affine'):
        tfm = AffineTransform()
        tfm.estimate(landmarks.iloc[:,:2],landmarks.iloc[:,2:4])
        msi_coords_tfm = pd.DataFrame(matrix_transform(
            msi_coords[['x', 'y']],
            tfm.params))
    elif (transform=='TPS'):
        tfm = ThinPlateSplineTransform()
        tfm.estimate((landmarks.iloc[:,:2]).to_numpy(),(landmarks.iloc[:,2:4]).to_numpy())
        msi_coords_tfm = pd.DataFrame(tfm(msi_coords[['x','y']].to_numpy()))
    msi_coords_tfm.index = msi_coords.index
    return msi_coords_tfm

def map_coords_HE2HE(file,msi_coords,transform):
    # read landmark and MSI coordinate files as well as Visium H&E (for shape to transform MSI H&E)
    visium_he_img = imread(file+'/visium/spatial/tissue_hires_image.png')

    if os.path.isfile(file+'/msi/MSI_HE_modified.jpg'):
        msi_he_img = imread(file+'/msi/MSI_HE_modified.jpg')
    else:
        msi_he_img = imread(file+'/msi/MSI_HE.jpg')
    landmarks = pd.read_csv(file+'/landmarks_HE2HE.csv')
    rows, cols = visium_he_img.shape[:2]

    # apply affine or TPS transform to coordinates and MSI H&E
    if (transform=='affine'):
        tfm = AffineTransform()
        tfm.estimate(landmarks.iloc[:,:2],landmarks.iloc[:,2:4])
        msi_coords_tfm = pd.DataFrame(matrix_transform(
            msi_coords[['x', 'y']],
            tfm.params))
        tfm = AffineTransform()
        tfm.estimate(landmarks.iloc[:,2:4],landmarks.iloc[:,:2])

    elif (transform=='TPS'):
        tfm = ThinPlateSplineTransform()
        tfm.estimate((landmarks.iloc[:,:2]).to_numpy(),(landmarks.iloc[:,2:4]).to_numpy())
        msi_coords_tfm = pd.DataFrame(tfm(msi_coords[[0,1]].to_numpy()))
        tfm = ThinPlateSplineTransform()
        tfm.estimate((landmarks.iloc[:,2:4]).to_numpy(),(landmarks.iloc[:,:2]).to_numpy())

    transformed_image = warp(msi_he_img, tfm,output_shape=(rows,cols))
    msi_coords_tfm.index = msi_coords.index
    # return both the transformed coordinates and transformed H&E image
    return {'transformed_coords':msi_coords_tfm,'msi_he_image':transformed_image}

# check whether there is an MSI H&E image and use the transformation pipeline depending on result
def apply_mapping(file):

    if os.path.isfile(file+'/msi/MSI_HE.jpg'):
        intermediate_coords = map_coords_MSI2HE(file,'affine')
        return(map_coords_HE2HE(file,intermediate_coords,snakemake.params['transform_type']))
    else:
        return(map_coords_noHE(file,snakemake.params['transform_type']))

# run full coregistration pipeline on selected sample
def run_coreg(sample):

    # get new MSI coordinates and image
    transformed_result = apply_mapping('input/'+sample)
    transformed_coords = transformed_result['transformed_coords']
    msi_he_image = transformed_result['msi_he_image']

    # save H&E image (either MSI H&E if available or Visium otherwise)
    fig, ax = plt.subplots(nrows=1, ncols=1 )
    if os.path.isfile('input/'+sample+'/msi/MSI_HE.jpg'):
        out_image = msi_he_image
    else:
        visium_he_img = imread('input/'+sample+'/visium/spatial/tissue_hires_image.png')
        out_image = visium_he_img
    plt.imsave(arr=out_image,fname='output/'+sample+'/transformed.png')

    transformed_coords.columns = ['x','y']
    # plot coordinates on top of H&E image
    plt.imshow(out_image)
    plt.scatter(x=transformed_coords['x'],y=transformed_coords['y'],s=0.1,c='r',alpha=1)
    fig.savefig('output/'+sample+'/transformed_withCoords.png')

    # save transformed coordinates
    if os.path.isfile('input/'+sample+'/msi/MSI_metadata_modified.csv'):
        msi_coords = pd.read_csv('input/'+sample+'/msi/MSI_metadata_modified.csv')
    else:
        msi_coords = pd.read_csv('input/'+sample+'/msi/MSI_metadata.csv')
    transformed_coords['spot_id']=msi_coords['spot_id']
    transformed_coords = transformed_coords[['spot_id','x','y']]
    transformed_coords.to_csv('output/'+sample+'/transformed.csv',index=False)

def main():
    # run on current sample
    run_coreg(snakemake.params['sample'])

if __name__ == "__main__":
    main()

