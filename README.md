# Chromatin_modules  across cell types
This repository contains the code for the paper "Non-coding variants impact cis-regulatory coordination in a cell type-specific manner".

_**Chromatin modules**_ (CMs) represent sub-TAD hubs encompassing interactions
between cis-regulatory elements, such as promoters and enhancers.
CMs are charachterized by
- interaction of active promoters with one or more active enhancers simultaneously
(median of two-four enhancers),
- preferential scaling at a sub-TAD level and therefore often being <100 kb
in size (with median sizes ranging from 32 to 70 kb),
- physical proximity of interacting elements within modules that show
enrichment in 3D contacts,
- deposition of histone modifications and binding of TFs in such regions
in a coordinated fashion with transcriptional
activity being a strong predictor of module formation.

For more details, please refer to [van Mierlo, Pushkarev, Trends in Genetics, 2023](https://www.cell.com/trends/genetics/fulltext/S0168-9525(22)00290-6).

This repository contains:
1. Input data used for CM mapping
2. Scripts for CM mapping with [VCMtools](https://doi.org/10.1016/j.cell.2015.08.001), [Clomics](https://www.science.org/doi/10.1126/science.aat8266) and [PHM](https://www.nature.com/articles/s41588-018-0278-6)
3. Methods for CM reproducibility/similarity scoring
4. CM simulation strategy
 
