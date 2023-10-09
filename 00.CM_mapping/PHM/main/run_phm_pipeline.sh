#!/bin/bash
# +
Help()
{
   # Display Help
    echo "Defining PHM input and execuxsting the method."
    echo "Usage: run_PHM_pipeline.sh [-p|s|d|c|t|o]"
    echo "  -p          Path to PHM methods"
    echo "  -s          Path to PHM output processing scripts"
    echo "  -d          Path to input cell type/dataset name"
    echo "  -c          Input chromosome/chromosomes (string)"
    echo "  -t          Time the script (bool 1 (YES) / 0 (NO))"
    echo "-----------------------------------------------------"
    echo "  -o          Path to output directory"
}

TIMING=1;

# Get the options
while getopts ":p:s:d:c:t:o:h" option;
do
    case ${option} in
        h) # display Help
            Help
            exit;;
        # ... other cases ...
        p) PHMDIR=${OPTARG};;
        s) PATH_TO_SCRIPTS=${OPTARG};;
        d) DATASET=${OPTARG};;
        c) CHROMOSOMES=${OPTARG};;
        t) TIMING=${OPTARG};;
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

    ###
    if [[ ${TIMING} == 1 ]];
    then
        echo "Timing the PHM script..."
        /usr/bin/time \
            -f 'Elapsed time: %es\nMemory usage: %M KB\nCPU usage: %P' \
            -o ${OUTPUT}/${chromosome}_memory_and_time_usage.txt \
            sh ${PATH_TO_SCRIPTS}/PHM.sh \
                -d $PHMDIR \
                -p $PEAK \
                -f $FPKM \
                -v $VCF \
                -i $IN1 \
                -n $IN2 \
                -o $OUT1 \
                -u $OUT2 \
                -l $VL \
                -r $PI1
    else
        sh ${PATH_TO_SCRIPTS}/PHM.sh \
            -d $PHMDIR \
            -p $PEAK \
            -f $FPKM \
            -v $VCF \
            -i $IN1 \
            -n $IN2 \
            -o $OUT1 \
            -u $OUT2 \
            -l $VL \
            -r $PI1
    fi
done
