# MAGPIE Pipeline Documentation

Welcome to the documentation for our **Snakemake pipeline**. This pipeline is designed to process spatial multi-omics MSI and Visium data, starting from unaligned coordinates for both modalities and outputting aligned coordinates and aligned observations. The outputs are stored in a spaceranger-style format which can be read by a range of spatial ecosystems.

More details about the pipeline can be found in the manuscript which has been published in Nature Communications:

```{admonition} Manuscript
:class: tip

**Spatially-resolved integrative analysis of transcriptomic and metabolomic changes in tissue injury**

Eleanor C. Williams, Lovisa Franzén, Martina Olsson Lindvall, Gregory Hamm, Steven Oag,Muntasir Mamun Majumder, James Denholm, Azam Hamidinekoo, Javier Escudero Morlanes, Marco Vicari, Joakim Lundeberg, Laura Setyo, Trevor M. Godfrey, Livia S. Eberlin, Aleksandr Zakirov, Jorrit J. Hornberg, Marianna Stamou, Patrik L. Ståhl, Anna Ollerstam, Jennifer Y. Tan, Irina Mohorianu.
*Nature Communications*, 17, Article 205. (2026)  
https://doi.org/10.1038/s41467-025-68003-w

Recent developments in spatially resolved -omics have enabled the joint study of gene expression, metabolite levels and tissue morphology, offering greater insights into biological pathways. 
Integrating these modalities from matched tissue sections to probe spatially-coordinated processes, however, remains challenging. 
Here we introduce _MAGPIE_, a framework for co-registering spatially resolved transcriptomics, metabolomics, and tissue morphology from the same or consecutive sections. 
We show _MAGPIE_’s generalisability and scalability on spatial multi-omics data from multiple tissues, combining Visium with MALDI and DESI mass spectrometry imaging. 
_MAGPIE_ was also applied to new multi-modal datasets generated with a specialised sampling strategy to characterise the metabolic and transcriptomic landscape in an in vivo model of drug-induced pulmonary fibrosis and to link small-molecule co-detection with endogenous lung responses. 
_MAGPIE_ demonstrates the refined resolution and enhanced interpretability that spatial multi-modal analyses provide for studying tissue injury especially in pharmacological contexts, and delivers a modular, accessible workflow for data integration.


```
This documentation is divided into sections explaining each stage of the pipeline, the rules, the configurations, and the results.

