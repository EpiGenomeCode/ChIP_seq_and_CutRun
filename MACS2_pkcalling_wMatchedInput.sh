#!/bin/bash

# Yingqian Zhan
# Jun 19, 2019
# CER, MSKCC

usage()
{
    echo "Usage: MACS2_pkcalling_wMatchedInput.sh --Treatment 'BAM1,BAM2,...' --Control 'BAMc1,BAMc2,...' --gsize --name --broadpk"
    echo ""
    echo "Description: This script will call peaks with matched control(s)"
    echo ""
}

PARAMS=""
while (( "$#" )); do
    case "$1" in
	-h | --help)
	    usage
	    exit 1
	    ;;
    --Treatment)
        Treatment=$2
        shift 2
        ;;
	--Control)
	    control=$2
	    shift 2
	    ;;
	--gsize)
        gsize=$2
        shift 2
        ;;
	--name)
	    name=$2
	    shift 2
	    ;;
	--broadpk)
	    broadpk=true
	    shift 
	    ;;
	--)
	    shift
	    break
	    ;;
	*)
	PARAMS="$PARAMS $1"
	shift
	;;
 esac
done
eval set -- "$PARAMS"

####
# Exit if file doens't exists
#### 
if [[ ${Treatment} == "" ]]; then
    echo "ERROR Treatment bam file(s) doens't exists"
    echo "Check ${Treatment}"
    usage
    exit 1
else
	echo "Treatment: ${Treatment}"
fi

if [[ ${control} == "" ]]; then
	echo "ERROR control bam file(s) doens't exists"
    echo "Check ${control}"
    usage
    exit 1
else
	echo "Control: $control"
fi

echo "gsize: ${gsize}"
echo "name: ${name}"
if [[ "${broadpk}" == true ]]; then
    echo "call broad peaks"
	echo "MACS2_pkcalling_wMatchedInput.sh --Treatment ${Treatment} --Control ${control} -f BAM --broad --gsize ${gsize} --name ${name} --broadpk"
else
    echo "call regular peaks"
	echo "MACS2_pkcalling_wMatchedInput.sh --Treatment ${Treatment} --Control ${control} --gsize ${gsize} --name ${name}"
fi

source /data/cer/shared/venvs/cer_base/bin/activate

module add R/R-3.3.3
samtools=/programs/x86_64-linux/samtools/1.4.1/bin/samtools
picard=/home/kocher/Programs/picard/build/libs/picard.jar
MACS2=/data/cer/shared/venvs/cer_base/bin/macs2

Treatments=`echo ${Treatment}|sed 's/\,/ /g'`
controls=`echo ${control}|sed 's/\,/ /g'`

echo "=========================================================="
echo "Starting on : $(date)"
echo "Current directory : $(pwd)"
echo "=========================================================="

if [[ "${broadpk}" == true ]]; then
	#echo $MACS2 callpeak -t ${Treatments} -c ${controls} -f BAM --broad -g ${gsize} --nomodel --pvalue 0.001 --broad-cutoff 0.005 -n ${name}
	$MACS2 callpeak -t ${Treatments} -c ${controls} -f BAM --broad -g ${gsize} --nomodel --pvalue 0.001 --broad-cutoff 0.001 -n ${name}
	# the --pvalue is for narrow peaks, --broad-cutoff is then a pvalue cut off for broad peaks. 
else
	#echo $MACS2 callpeak -t ${Treatments} -c ${controls} -f BAM -g ${gsize} --nomodel --pvalue 0.001 -n ${name}
	$MACS2 callpeak -t ${Treatments} -c ${controls} -f BAM -g ${gsize} --nomodel --pvalue 0.001 -n ${name}
fi

echo "=========================================================="
echo "Finished on : $(date)"
echo "=========================================================="