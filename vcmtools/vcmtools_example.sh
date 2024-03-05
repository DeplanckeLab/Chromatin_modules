# +
## Specify output folder
output_path="./vcmtools_test";

## Specify paths to VCMtools scripts
path_to_main="./vcmtools/main";

## define VCMtools-specific parameters
max_peak_dist=0.5;  ## Maximum distance between two peak centers
pv_thresholds="0.0005,0.001"  ## p-values to use for CM mapping

## Specify the number of processs to use
## if parallelizing calculation of correlations (highly recommendded)
n_cores=1;  ## default n_cores=1

## Specify chromosomes to use for CM mapping
chromosomes="22";  ## dash defines a range, default is "1-22" (all chromosomes from chr1 to chr22);
## Alternatively, "chromosomes" parameter can be defined as
## "1,5,10" (without spaces) when using, in this example, three chromosomes;
## or "22", when using a single chromosome

## Specify cell type/dataset name
dataset_name="test_data_chr22";

## Specify paths to count matrices for ChIP-seq/ATAC-seq data
input_path="./test_data/gzipped_bed";
k4me1_matrix="H3K4me1:${input_path}/H3K4me1_chr22.bed.gz";
k27ac_matrix="H3K27ac:${input_path}/H3K27ac_chr22.bed.gz";
# -

sh ${path_to_main}/pyVCMtools.sh \
    -s ${path_to_main}/src \
    -d ${dataset} \
    -m ${max_peak_dist} \
    -p ${pv_thresholds} \
    -f "${k4me1_matrix},${k27ac_matrix}" \
    -c ${chromosomes} \
    -n ${n_cores} \
    -o output_path

