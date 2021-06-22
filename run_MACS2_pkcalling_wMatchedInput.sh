#!/bin/bash

# Modified from MSKCC CER
# This script will run peak calling with Matched input/control(s).
# Usage: run_MACS2_pkcalling_wMatchedInput.sh --File <File> --genome <genome> --prjDir <output dir>
#
#    --File: A text file containing three columns 'treatment	control	broadpk'
#		The first column is the treatment file(s), ',' separated if multiple;
#		The second column is the control file(s), ',' separated if multiple;
#		The third column tells whether enable broad peak calling, 'Yes | No';
#		The header should be commented out '#'
#
#   --genome: genome name ('hs' or 'mm')
#
#   --prjDir: Output directory 

usage()
{
    echo "Usage: run_MACS2_pkcalling_wMatchedInput.sh --File <File> --genome <genome> --prjDir <output dir>"
    echo ""
    echo "Description: This script will call peaks from ChIP/ATAC with matched control(s)"
	echo "--File:" 
	echo " A text file containing three columns <treatment control broadpk>"
	echo "  1st column: the treatment file(s), ',' separated if multiple"
	echo "  2nd column: the control file(s), ',' separated if multiple"
	echo "  3rd column: 'Yes or No', tells whether enable broad peak calling"
	echo " --genome: 'hs' or 'mm' (as gsize in MACS2 peak calling)"
	echo "Option:"
	echo " --prjDir: Output directory"
	echo ""
    echo ""
}

PARAMS=""
while (( "$#" )); do
    case "$1" in
	-h | --help)
	    usage
	    exit 1
	    ;;
	--File)
        File=$2
        shift 2
        ;;
    --prjDir)
        prjDir=$2
        shift 2
        ;;
	--genome)
        genome=$2
        shift 2
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
#File=${PARAMS[0]}

####
# Exit if file doens't exists
#### 
if [[ ! -f ${File} ]]; then
    echo "Require a information file"
    usage
    exit 1
else
	echo "Information file: ${File}"
fi

if [[ ${genome} == "" ]]; then
	echo "Genome is not specified"
    usage
    exit 1
else
	echo "Genome: $genome"
fi

if [[ ${prjDir} == "" ]]; then
	${prjDir} = `pwd`
	echo "Output: ${prjDir}"
else
	echo "Output: ${prjDir}"
fi

####
# Main
#### 
source /data/cer/shared/venvs/cer_base/bin/activate
pkcall=/data/cer/shared/code/CER_scripts/pipeline/MACS2_pkcalling_wMatchedInput.sh
outdir=${prjDir}/MACS2_MatchedInput
mkdir -p ${outdir}
cd ${outdir}

echo
echo "========================================================================="
echo "submitting MACS2 peak calling with control peaks: $(date)"
echo "========================================================================="
while IFS= read -r line
do
	treat=`echo $line|awk '{print $1}'`
	bambase=$(basename "$treat" .bam)
	macsout=MACS2_vsMatchedInput_${bambase}
	ctrl=`echo $line|awk '{print $2}'`
	broad_tag=`echo $line|awk '{print $3}'`

	orig_dir=`dirname $treat`
	orig_dirname=`basename ${orig_dir}`
	mkdir -p ${orig_dirname}
	cd ${orig_dirname}

	if [[ $broad_tag == "Yes" ]];then
		echo 'submitting ${macsout} ...'
		bsub -o ${macsout}.log -n 2 -W 2:00 -R "rusage[mem=10]" -J ${macsout} \
        	$pkcall --Treatment ${treat} \
                --Control ${ctrl} \
                --gsize $genome \
                --name ${macsout} \
				--broadpk
	else
		echo 'submitting ${macsout} ...'
		bsub -o ${macsout}.log -n 2 -W 2:00 -R "rusage[mem=10]" -J ${macsout} \
        	$pkcall --Treatment ${treat} \
                --Control ${ctrl} \
                --gsize $genome \
                --name ${macsout} 
	fi
	cd ${outdir}
done < <(grep -v '#' ${File})

echo
echo "========================================================================="
echo "Complete: $(date)"
echo "========================================================================="
