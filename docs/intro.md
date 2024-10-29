# MAGPIE Pipeline Documentation

Welcome to the documentation for our **Snakemake pipeline**. This pipeline is designed to process spatial multi-omics MSI and Visium data, starting from unaligned coordinates for both modalities and outputting aligned coordinates and aligned observations. The outputs are stored in a spaceranger-style format which can be read by a range of spatial ecosystems.

More details about the pipeline can be found in the manuscript:

```{admonition} Manuscript
:class: tip

**Integrative analysis of spatial transcriptomics, metabolomics, and histologic changes illustrated in tissue injury studies**

*Eleanor C Williams, Lovisa Franzén, Martina Olsson Lindvall, Gregory Hamm, Steven Oag, Jim Denholm, Azam Hamidinekoo, Muntasir Mamun Majumder, Javier Escudero Morlanes, Marco  Vicari, Laura Setyo, Aleksandr Zakirov, Jorrit J Hornberg, Anna Ollerstam, Patrik L Ståhl, Marianna Stamou, Jennifer Tan, Irina Mohorianu*

>Recent developments in spatially resolved omics have expanded studies linking gene expression, epigenetic alterations, protein levels, and metabolite intensity to tissue histology. The integration of multiple spatial measurements can offer new insights into alterations modulations of expression propagating across modalities;, however, it also presents experimental and computational challenges.  
>To set the multimodal datasets into a shared coordinate system for enhanced integration and analysis, we propose MAGPIE (Multi-modal Alignment of Genes and Peaks for Integrated Exploration), a framework for co-registering spatially- resolved transcriptomics and spatial metabolomics measurements on the same or consecutive tissue sections, present within their existing histological context. Further, we showcase the utility of the MAGPIE framework on spatial multi-omics data from lung tissue, an inherently heterogeneous tissue type with integrity challenges and for which we developed an experimental sampling strategy to allow multimodal data  generation. In these case studies, we were able to link pharmaceutical co-detection with endogenous responses in rat lung tissue following inhalation of a small molecule, which had previously been stopped during preclinical development with findings of lung irritation, and to characterise the metabolic and transcriptomic landscape in a mouse model of drug-induced pulmonary fibrosis in conjunction with histopathology annotations. 
>The generalisability and scalability of the MAGPIE framework were further benchmarked on public datasets from multiple species and tissue types, demonstrating applicability to both DESI and MALDI mass spectrometry imaging together with Visium-enabled transcriptomic assessment. MAGPIE highlights the refined resolution and increased interpretability of spatial multimodal analyses in studying tissue injury, particularly in a pharmacological context, and offers a modular, accessible computational workflow for data integration
```
This documentation is divided into sections explaining each stage of the pipeline, the rules, the configurations, and the results.

