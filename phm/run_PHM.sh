#!/bin/bash
# +
Help()
{
   # Display Help
    echo "Defining PHM input and execuxsting the method."
    echo "Usage: run_PHM.sh [-p|s|d|c|t|o]"
    echo "  -d          Path to input cell type/dataset name"
    echo "  -c          Input chromosome/chromosomes (string)"
    echo "  -p          Path to PHM methods"
    echo "  -s          Path to PHM output processing scripts"
    echo "  -t          Posterior probability thresholds for causal peak interactions"
    echo "-----------------------------------------------------"
    echo "  -o          Path to output directory"
}

# Get the options
while getopts ":p:s:d:c:o:h" option;
do
    case ${option} in
        h) # display Help
            Help
            exit;;
        c) CHROMOSOMES=${OPTARG};;
        d) DATASET=${OPTARG};;
        p) PHMDIR=${OPTARG};;
        s) PATH_TO_HELPER=${OPTARG};;
        t) THRESHOLDS=${OPTARG};;
        o) OUTPUT=${OPTARG};;

        \?) # Invalid option
            echo "Error: Invalid option"
            exit;;
    esac
done


IFS=','

# Create an array to store the split values
read -ra CHROMOSOMES_LIST <<< ${CHROMOSOMES}

IFS=' '


for chromosome_el in ${CHROMOSOMES_LIST[@]};
do

    if [[ ${chromosome_el} != *"chr"* ]];
    then
        chromosome="chr"${chromosome_el}
    else
        chromosome=${chromosome_el}
    fi
    
    cd ${OUTPUT}/${DATASET}/phm_results
    echo ${OUTPUT}/${DATASET}/phm_results

    # Bayes factors calculated by script/bayeslm1.sh
    IN1=${chromosome}/bayeslm_1.gz
    #
    # Output directory of hm
    OUT1=${chromosome}/hm_output/
    #
    # Bayes factors calculated by script/bayeslm2.sh
    IN2=${chromosome}/bayeslm_2.gz
    #
    # Output directory of phm
    OUT2=${chromosome}/phm_output/
    #
    # Log normalised read counts (binary double array)
    FPKM=${chromosome}/normalized_counts_${chromosome}.bin
    #
    # VCF file (tabix indexed)
    VCF=${chromosome}/${DATASET}_${chromosome}_no_chr.vcf.gz
    #
    # Peak annotation (tabix indexed)
    PEAK=${chromosome}/${chromosome}_peak_coordinates.bed.gz
    #
    # Hyper-parameter estimate of variant level prior
    VL=${chromosome}/hm_output/variant_level.bin
    #
    # Prior probability that each peak is a QTL
    PI1=${chromosome}/hm_output/Pi1.bin


    # First stage
    echo "First stage"
    NF=`cat $PEAK | gunzip | wc -l | awk '{print $1}'`
    sh $PHMDIR/script/bayeslm1.sh 1 $NF $FPKM $VCF $PEAK $IN1
    R1=`cat $IN1 | gunzip | wc -l | awk '{print $1}'`
    hm -i $IN1 -c I,S,C1,S,S,S,S,S,C2,C3,S,S,B -r $R1 -f $NF -p -o $OUT1 -v

    # Second stage
    echo "Second stage"
    sh $PHMDIR/script/bayeslm2.sh 1 $NF $FPKM $VCF $PEAK $VL $PI1 $IN2
    R2=`cat $IN2 | gunzip | wc -l | awk '{print $1}'`
    phm -i $IN2 -c J,K,N4,B10 -r $R2 -f $NF -o $OUT2
done
# -

## Convert DAGs into CM format
for threshold in ${THRESHOLDS[*]};
do
    python3.9 ${PATH_TO_HELPER}/1.build_DAGs_per_chr.py \
        -d ${dataset} \
        -i ${OUTPUT}/${dataset}/phm_results \
        -t ${threshold} \
        -c ${CHROMOSOMES}
    python3.9 ${PATH_TO_HELPER}/2.get_CMs_from_DAGs.py \
        -d ${dataset} \
        -i ${OUTPUT}/${dataset}/phm_results \
        -t ${threshold} \
        -c ${CHROMOSOMES}
done
