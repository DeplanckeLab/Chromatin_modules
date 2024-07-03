## # Specify path to the installed PHM model (https://github.com/natsuhiko/PHM)
export PHMDIR=/data/software/PHM-0.2;

## Assuming the user runs code from the same folder, where 1.phm_example.sh is
## # Path to PHM helper scripts
path_to_main="./main";

## # Specify dataset name
dataset="test_data";

## # Specify chromosomes and the respective string
chromosomes=( "chr22" );
chromosome_str="22"

## # For a range of chromosomes, we suggest using a list of thresholds as follows:
## chromosomes=( "chr20" "chr21" "chr22" );
## chromosome_str="20-22"
## # Or equivalently:
## chromosome_str="20,21,22"

window=0.5
window_int=$(echo "${window} * 2000000" | bc);
## # threshold for posterior probability
## # of causal interaction for DAG construction:
threshold="0.8";
## # For a list of thresholds, add external loop 
## thresholds=( "0.7" "0.75" "0.8" );

## # Specify path to input data
data_path="../1.mapped_CMs/phm";
## # Specify output path
output_path="../1.mapped_CMs/phm";

## # data_path should point to the folder with files generated
## # with the 0.phm_data_preparation.ipynb notebook.
## # The output folder might be the same or different one,
## # depends on the user's choice.

for chromosome in ${chromosomes[*]};
do
    mkdir -p ${output_path}/${chromosome}/phm_output;
    mkdir -p ${output_path}/${chromosome}/hm_output;

    ## # INPUT DATA
    ## # Log normalised read counts (binary double array)
    FPKM=${data_path}/${chromosome}/normalized_counts_${chromosome}.bin

    ## # VCF file (tabix indexed)
    VCF=${data_path}/${chromosome}/${dataset}_${chromosome}_no_chr.vcf.gz

    ## # Peak annotation (tabix indexed)
    PEAK=${data_path}/${chromosome}/${chromosome}_peak_coordinates.bed.gz
    
    ## # OUTPUT DATA
    ## # Bayes factors calculated by script/bayeslm1.sh
    IN1=${output_path}/${chromosome}/bayeslm_1.gz

    ## # Output directory of hm
    OUT1=${output_path}/${chromosome}/hm_output/

    ## # Bayes factors calculated by script/bayeslm2.sh
    IN2=${output_path}/${chromosome}/bayeslm_2.gz

    ## # Output directory of phm
    OUT2=${output_path}/${chromosome}/phm_output/

    ## # Hyper-parameter estimate of variant level prior
    VL=${output_path}/${chromosome}/hm_output/variant_level.bin

    ## # Prior probability that each peak is a QTL
    PI1=${output_path}/${chromosome}/hm_output/Pi1.bin
    ## # Run main PHM script
    sh ${path_to_main}/PHM.sh -d $PHMDIR -s ${path_to_main}/src -w ${window_int} -p $PEAK -f $FPKM -v $VCF -i $IN1 -n $IN2 -o $OUT1 -u $OUT2 -l $VL -r $PI1
done

## # threshold defines causal interaction threshold for DAG construction
## # Reconstruct DAGs from the pairwise peak associations
python3 ${path_to_main}/src/1.build_DAGs_per_chr.py -d ${dataset} -i ${data_path} -t ${threshold} -c ${chromosome_str} -o ${output_path}

## # Convert DAG files into CM .tracks and .content files
python3 ${path_to_main}/src/2.get_CMs_from_DAGs.py -d ${dataset} -i ${data_path} -t ${threshold} -c ${chromosome_str} -o ${output_path}
