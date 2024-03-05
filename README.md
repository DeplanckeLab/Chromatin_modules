This repository contains the code for the paper "Non-coding variants impact cis-regulatory coordination in a cell type-specific manner".
> [!Chromatin modules (CMs)]
> _**Chromatin modules**_ (CMs) represent sub-TAD hubs encompassing interactions
between cis-regulatory elements, such as promoters and enhancers.
> CMs are charachterized by
> - interaction of active promoters with one or more active enhancers simultaneously
(median of two-four enhancers),
> - preferential scaling at a sub-TAD level and therefore often being <100 kb
in size (with median sizes ranging from 32 to 70 kb),
> - physical proximity of interacting elements within modules that show
enrichment in 3D contacts,
> - deposition of histone modifications and binding of TFs in such regions
in a coordinated fashion with transcriptional
activity being a strong predictor of module formation.

> For more details, please refer to [this review](https://www.cell.com/trends/genetics/fulltext/S0168-9525(22)00290-6).

Explore CMs [online](https://chromo.epfl.ch/)!

This repository contains:
1. Input data used for CM mapping
2. Scripts for CM mapping with [VCMtools](https://doi.org/10.1016/j.cell.2015.08.001), [Clomics](https://www.science.org/doi/10.1126/science.aat8266) and [PHM](https://www.nature.com/articles/s41588-018-0278-6)

# Data preparation and CM mapping

**I. Preparation steps:**

1. Clone the Chromatin modules repository
2. Data preparation
    i. Download [test data](https://github.com/DeplanckeLab/Chromatin_modules/tree/main/test_data) for [LCLs](), chr22, or use your own data
    ii. Make sure to have ChIP-seq data (*used in all methods*) in the following format:

| #Chr | start    | end      | pid                     | did                     | strand | sample_id_1 | sample_id_2 | ... |
| ---- | -------- | -------- | ----------------------- | ----------------------- | ------ | ----------- | ----------- | --- |
| 22   | 16192326 | 16192745 | chr22:16192326:16192745 | chr22:16192326:16192745 | +      | -0.212      | -0.175      | ... |
| 22   | 16204838 | 16205246 | chr22:16204838:16205246 | chr22:16204838:16205246 | +      | -0.221      | 0.339       | ... |
| 22   | 16293852 | 16294075 | chr22:16293852:16294075 | chr22:16293852:16294075 | +      | -0.038      | 0.0989      | ... |


3. Genotype data in VCF files (*used in PHM*), see section 3 of [this]() notebook on PHM-specific data preparation

**II. CM mapping:**
1. Clomics
    i. Install [Clomics-2.0](https://github.com/OlgaPushkarev/clomics-2.0)
    ii. Adjust paths and parameters in the clomics_example.sh file
2. VCMtools
3. PHM
    i. Install [PHM](https://github.com/natsuhiko/PHM)
    ii. Prepare data in a PHM-specific format, following the steps explained in [this notebook]().
**III. CM exploration:**
1. Explore CMs from *[...]CM.tracks.bed* in [IGV](https://igv.org). Information on CM peaks can be found in *[...]CM.content.txt*


_**Upcoming updates**_
1. Methods for CM reproducibility/similarity scoring
2. CM simulation strategy
 
