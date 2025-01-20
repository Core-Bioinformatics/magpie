# MAGPIE Pipeline Documentation

Welcome to the documentation for our **Snakemake pipeline**. This pipeline is designed to process spatial multi-omics MSI and Visium data, starting from unaligned coordinates for both modalities and outputting aligned coordinates and aligned observations. The outputs are stored in a spaceranger-style format which can be read by a range of spatial ecosystems.

More details about the pipeline can be found in the manuscript:

```{admonition} Manuscript
:class: tip

**Spatially-resolved integrative analysis of transcriptomic and metabolomic changes in tissue injury**

*Eleanor C Williams, Lovisa Franzén, Martina Olsson Lindvall, Gregory Hamm, Steven Oag, Jim Denholm, Azam Hamidinekoo, Muntasir Mamun Majumder, Javier Escudero Morlanes, Marco  Vicari, Joakim Lundeberg, Laura Setyo, Aleksandr Zakirov, Jorrit J Hornberg, Anna Ollerstam, Patrik L Ståhl, Marianna Stamou, Jennifer Tan, Irina Mohorianu*

>Recent developments in spatially-resolved -omics have enabled studies linking gene expression and metabolite levels to tissue histology, offering new insights into modulations of expression propagating across modalities; however, cross-modality integration also presents experimental and computational challenges. 
>To align multimodal datasets into a shared coordinate system and facilitate enhanced integration and analysis, we propose MAGPIE (Multi-modal Alignment of Genes and Peaks for Integrated Exploration), a framework for co-registering spatially-resolved transcriptomics and metabolomics measurements on the same or consecutive tissue sections. 
>We illustrate MAGPIE on spatial multi-omics data from lung tissue, an inherently heterogeneous tissue with integrity challenges, for which we developed an experimental sampling strategy to allow multimodal data generation. Moreover, we showcase the linking of pharmaceutical co-detection with endogenous responses in lung tissue and characterise the metabolic and transcriptomic landscape in a model of drug-induced pulmonary fibrosis.
>The generalisability and scalability of MAGPIE were further benchmarked on public datasets across multiple species and tissue types, illustrating its applicability to both DESI and MALDI mass spectrometry imaging coupled with Visium-enabled transcriptomic assessment. MAGPIE highlights the refined resolution and increased interpretability of spatial multimodal analyses in studying tissue injury, particularly in a pharmacological context, and offers a modular, accessible computational workflow for data integration.
```
This documentation is divided into sections explaining each stage of the pipeline, the rules, the configurations, and the results.

