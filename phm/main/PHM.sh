Help()
{
   ## Display Help
   echo "Run PHM pipeline."
   echo "Syntax: PHM.sh [-d|s|w|p|f|v|i|n|o|u|l|r]"
}

## Get the options
while getopts ":d:s:w:p:f:v:i:n:o:u:l:r:" option; do
    case $option in
        h) ## display Help
            Help
            exit;;
        d) PHMDIR=${OPTARG};;
        s) PATH_TO_SRC=${OPTARG};;
        w) WINDOW=${OPTARG};;
        p) PEAK=${OPTARG};;
        f) FPKM=${OPTARG};;
        v) VCF=${OPTARG};;
        i) IN1=${OPTARG};;
        n) IN2=${OPTARG};;
        o) OUT1=${OPTARG};;
        u) OUT2=${OPTARG};;
        l) VL=${OPTARG};;
        r) PI1=${OPTARG};;
        \?) ## Invalid option
            echo "Error: Invalid option"
            exit;;
    esac
done

## First stage
echo "First stage"
NF=`cat $PEAK | gunzip | wc -l | awk '{print $1}'`
sh ${PATH_TO_SRC}/bayeslm1.sh \
    1 \
    $NF \
    $FPKM \
    $VCF \
    $PEAK \
    $IN1 \
    $WINDOW

R1=`cat $IN1 | gunzip | wc -l | awk '{print $1}'`
hm -i $IN1 -c I,S,C1,S,S,S,S,S,C2,C3,S,S,B -r $R1 -f $NF -p -o $OUT1 -v

## Second stage
echo "Second stage"
sh ${PATH_TO_SRC}/bayeslm2.sh \
    1 \
    $NF \
    $FPKM \
    $VCF \
    $PEAK \
    $VL \
    $PI1 \
    $IN2 \
    $WINDOW
    
R2=`cat $IN2 | gunzip | wc -l | awk '{print $1}'`
phm -i $IN2 -c J,K,N4,B10 -r $R2 -f $NF -o $OUT2
