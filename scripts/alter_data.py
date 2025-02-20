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
import glob

# Identify transform needed to map from MSI dim reduction to Visium H&E and apply it
def map_coords_noHE(sample,transform):

    # read landmark and MSI coordinate files
    landmarks = pd.read_csv('input/'+sample+'/landmarks_noHE.csv')
    if os.path.isfile('input/'+sample+'/msi/MSI_metadata_modified.csv'):
        msi_coords = pd.read_csv('input/'+sample+'/msi/MSI_metadata_modified.csv')
    else:
        msi_coords = pd.read_csv('input/'+sample+'/msi/MSI_metadata.csv')

    # apply affine or TPS transform to coordinates
    if (transform=='affine'):
        print('running affine')
        tfm = AffineTransform()
        tfm.estimate(landmarks.iloc[:,:2],landmarks.iloc[:,2:4])
        msi_coords_tfm = pd.DataFrame(matrix_transform(
            msi_coords[['x', 'y']],
            tfm.params))
    elif (transform=='TPS'):
        print('running TPS')
        tfm = ThinPlateSplineTransform()
        tfm.estimate((landmarks.iloc[:,:2]).to_numpy(),(landmarks.iloc[:,2:4]).to_numpy())
        msi_coords_tfm = pd.DataFrame(tfm(msi_coords[['x','y']].to_numpy()))
    msi_coords_tfm.index = msi_coords.index
    return {'transformed_coords':msi_coords_tfm,'msi_he_image':None}

# Identify transform needed to map from MSI dim reduction to MSI H&E and apply it
def map_coords_MSI2HE(sample,transform,msi_he_img):

    # read landmark and MSI coordinate files
    landmarks = pd.read_csv('input/'+sample+'/landmarks_MSI2HE.csv')

    if os.path.isfile('input/'+sample+'/msi/MSI_metadata_modified.csv'):
        msi_coords = pd.read_csv('input/'+sample+'/msi/MSI_metadata_modified.csv')
    else:
        msi_coords = pd.read_csv('input/'+sample+'/msi/MSI_metadata.csv')

    if os.path.isfile('input/'+sample+'/msi/MSI_dimreduction.csv'):
        msi_dimred = pd.read_csv('input/'+sample+'/msi/MSI_dimreduction.csv')

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

    fig, ax = plt.subplots(nrows=1, ncols=1 )
    out_image = imread(msi_he_img)
    transformed_coords = msi_coords_tfm.copy()
    transformed_coords.columns = ['x','y']
    transformed_coords['spot_id']=msi_coords['spot_id']
    print(transformed_coords[:5])
    if os.path.isfile('input/'+sample+'/msi/MSI_dimreduction.csv'):
        transformed_coords = pd.merge(transformed_coords,msi_dimred[['spot_id','color']],on='spot_id')
    else :
        transformed_coords['color']=1
    print(transformed_coords)
    # plot coordinates on top of H&E image
    plt.imshow(out_image)
    plt.scatter(x=transformed_coords['x'],y=transformed_coords['y'],c=transformed_coords['color'],s=0.1,alpha=0.5)
    fig.savefig('output/' + sample +'/MSI_HE_withCoords.png')

    return msi_coords_tfm

def map_coords_HE2HE(sample,msi_coords,transform,msi_he_img):
    # read landmark and MSI coordinate files as well as Visium H&E (for shape to transform MSI H&E)
    visium_he_img = imread('input/'+sample+'/visium/spatial/tissue_hires_image.png')

    if os.path.isfile('input/'+sample+'/msi/MSI_HE_modified.jpg'):
        msi_he_img = imread('input/'+sample+'/msi/MSI_HE_modified.jpg')
    else:
        msi_he_img = imread(msi_he_img)
    landmarks = pd.read_csv('input/'+sample+'/landmarks_HE2HE.csv')
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
def apply_mapping(sample,msi_he_img):

    if not (msi_he_img is None):
        intermediate_coords = map_coords_MSI2HE(sample,snakemake.params['MSI2HE_transform'],msi_he_img)
        return(map_coords_HE2HE(sample,intermediate_coords,snakemake.params['HE2HE_transform'],msi_he_img))
    else:
        return(map_coords_noHE(sample,snakemake.params['no_HE_transform']))

# run full coregistration pipeline on selected sample
def run_coreg(sample):

    if glob.glob('input/'+sample+'/msi/MSI_HE.*') != []:
        msi_he_img = glob.glob('input/'+sample+'/msi/MSI_HE.*')
        msi_he_img = [x for x in msi_he_img if any(ext in x for ext in ['tiff','png','jpg'])]
        if msi_he_img == []:
            msi_he_img = None
        else :
            msi_he_img = msi_he_img[0]
        print(msi_he_img)
    else :
        msi_he_img = None

    # get new MSI coordinates and image
    transformed_result = apply_mapping(sample,msi_he_img)
    transformed_coords = transformed_result['transformed_coords']
    msi_he_image = transformed_result['msi_he_image']
    transformed_coords.columns = ['x','y']

    # save H&E image (either MSI H&E if available or Visium otherwise)
    if not (msi_he_img is None):
        out_image = msi_he_image
        fig, ax = plt.subplots(nrows=1, ncols=1 )
        plt.imshow(out_image)
        plt.scatter(x=transformed_coords['x'],y=transformed_coords['y'],s=0.1,c='r',alpha=0.7)
        fig.savefig('output/'+sample+'/transformed_withCoords.png')
    else:
        visium_he_img = imread('input/'+sample+'/visium/spatial/tissue_hires_image.png')
        out_image = visium_he_img

    plt.imsave(arr=out_image,fname='output/'+sample+'/transformed.png')
    visium_he = imread('input/'+sample+'/visium/spatial/tissue_hires_image.png')
    fig, ax = plt.subplots(nrows=1, ncols=1)
    plt.imshow(visium_he)
    plt.scatter(x=transformed_coords['x'],y=transformed_coords['y'],s=0.1,c='r',alpha=0.7)
    fig.savefig('output/'+sample+'/transformed_withCoords_VisiumHE.png')

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

