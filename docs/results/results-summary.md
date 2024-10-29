# Pipeline Outputs

After running the pipeline, the following outputs are created:

- **Transformed MSI coordinates**: Comma-separated table with old and new MSI coordinates after mapping onto Visium coordinate space.
- **Spaceranger-style MSI object**: MSI object with newly transformed coordinates which can be read like a Visium object by multiple ecosystems
- **Visium spot-matched MSI object in spaceranger-style format**: MSI object after observations have been collapsed onto Visium spots to ensure multi-omics tools can be readily be applied (optional step)

## Reading spaceranger-style outputs

The resulting spaceranger-style objects can be readily imported by several spatial ecosystems, as shown below:

In the scanpy ecosystem for Python:
```python
import scanpy
msi_obj = scanpy.read_visium('output/'+sample+'/spaceranger/',library_id='myLib')
msi_obj = scanpy.read_visium('output/'+sample+'/spaceranger_meanIntensity/',library_id='myLib')
```

In the Seurat ecosystem for R:
```R
library(Seurat)
msi_obj = Load10X_Spatial(paste0('output/',sample,'/spaceranger'))
msi_obj = Load10X_Spatial(paste0('output/',sample,'/spaceranger_meanIntensity'))
```

Using the semla toolkit, you can even create a multi-modal object directly. For this, either the Visium spot-matched MSI data can be used or the standard spaceranger-style MSI object can be collapsed into Visium spots by *semla* itself.

In the semla ecosystem for R:
```R
library(semla)

# Create MSI object
data_root_directory = paste0('output/',sample,'/spaceranger/')
samples <- Sys.glob(paths = file.path(data_root_directory, 
                                      "filtered_feature_bc_matrix.h5"))
imgs <- Sys.glob(paths = file.path(data_root_directory, 
                                   "spatial", "tissue_hires_image.png"))
spotfiles <- Sys.glob(paths = file.path(data_root_directory, 
                                        "spatial", "tissue_positions_list.csv"))
json <- Sys.glob(paths = file.path(data_root_directory, 
                                   "spatial", "scalefactors_json.json"))
infoTable_msi <- tibble(samples, imgs, spotfiles, json, # Add required columns
                    sample_id = sample) # Add additional column

se_msi <- ReadVisiumData(infoTable = infoTable_msi, 
                         assay = "MSI", 
                         remove_spots_outside_HE = T, 
                         remove_spots_outside_tissue = T)

# Create Visium object
data_root_directory = paste0('input/',sample,'/visium/')
samples <- Sys.glob(paths = file.path(data_root_directory, 
                                      "filtered_feature_bc_matrix.h5"))
imgs <- Sys.glob(paths = file.path(data_root_directory, 
                                   "spatial", "tissue_hires_image.png"))
spotfiles <- Sys.glob(paths = file.path(data_root_directory, 
                                        "spatial", "tissue_positions_list.csv"))
json <- Sys.glob(paths = file.path(data_root_directory, 
                                   "spatial", "scalefactors_json.json"))
infoTable_visium <- tibble(samples, imgs, spotfiles, json, # Add required columns
                    sample_id = sample) # Add additional column

se_visium <- ReadVisiumData(infoTable = infoTable_visium, 
                            assay = "Visium", 
                            remove_spots_outside_HE = T, 
                            remove_spots_outside_tissue = T)

# Create combined multimodal object
se_mmo <- CreateMultiModalObject(object_ref = se_visium, 
                                 object_map = se_msi,
                                 agg_func = "mean",  # multiple MSI pixels per spot are combined by taking mean intensity
                                 new_assay_name = "MSI")
```

