# MAGPIE Pipeline Documentation

Welcome to the documentation for our **Snakemake pipeline**. This pipeline is designed to process spatial multi-omics MSI and Visium data, starting from unaligned coordinates for both modalities and outputting aligned coordinates and aligned observations. The outputs are stored in a spaceranger-style format which can be read by a range of spatial ecosystems.

More details about the pipeline can be found in the manuscript:

```{admonition} Manuscript
:class: tip

**Spatially-resolved integrative analysis of transcriptomic and metabolomic changes in tissue injury**

*Eleanor C Williams, Lovisa Franzén, Martina Olsson Lindvall, Gregory Hamm, Steven Oag, 
Muntasir Mamun Majumder, James Denholm, Azam Hamidinekoo, Javier Escudero Morlanes, Marco Vicari, Joakim Lundeberg, Laura Setyo, Aleksandr Zakirov, Jorrit J Hornberg, Marianna Stamou, Patrik L Ståhl, Anna Ollerstam, Jennifer Tan, Irina Mohorianu*

>Recent developments in spatially resolved -omics have enabled studies linking gene expression and metabolite levels to tissue morphology, offering new insights into biological pathways. By capturing multiple modalities on matched tissue sections, one can better probe how different biological entities interact in a spatially coordinated manner. However, such cross-modality integration presents experimental and computational challenges. 
>To align multimodal datasets into a shared coordinate system and facilitate enhanced integration and analysis, we propose MAGPIE (Multi-modal Alignment of Genes and Peaks for Integrated Exploration), a framework for co-registering spatially resolved transcriptomics, metabolomics, and tissue morphology from the same or consecutive sections. 
>We illustrate the generalisability and scalability of MAGPIE on spatial multi-omics data from multiple tissues, combining Visium with both MALDI and DESI mass spectrometry imaging. MAGPIE was also applied to newly generated multimodal datasets created using specialised experimental sampling strategy to characterise the metabolic and transcriptomic landscape in an in vivo model of drug-induced pulmonary fibrosis, to showcase the linking of small-molecule co-detection with endogenous responses in lung tissue.
>MAGPIE highlights the refined resolution and increased interpretability of spatial multimodal analyses in studying tissue injury, particularly in pharmacological contexts, and offers a modular, accessible computational workflow for data integration.

```
This documentation is divided into sections explaining each stage of the pipeline, the rules, the configurations, and the results.

