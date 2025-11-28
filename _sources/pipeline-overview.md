# Pipeline Overview

The pipeline consists of the following steps:

1. **Landmark selection**: Manual landmarks are selected using the MAGPIE shiny app
2. **Coregistration**: Coordinates from modalities are set into the same coordinate system
3. **Creation of spaceranger-style object**: A new spaceranger-style object is created for MSI data using coregistered coordinates
4. **Aligning observations across modalities**: Observations are matched between modalities so that multi-omics methods requiring common samples can be applied.

Here is a summary of the MAGPIE pipeline:

![pipeline](figures/pipeline_diagram.jpg)
