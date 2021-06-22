#!/bin/bash
# Modified from MSKCC CER
# File: run_norm.sh
#
# Usage: run_norm.sh <readsfile> <assay> <genome> <scale-factor>
#
# Arguments:
#    1. readsfile -- Path to alignment file.
#    2. assay -- Name of assay (ChIP or ATAC)
#    3. genome -- Name of genome assembly to use.
#    4. scale10mil -- Factor by which to multiple fragment count to get to 10 million fragments.
#
# Outputs:
#    1. bigWig file corresponding to given alignment file.
#    2. Outputs from run_tss.sh
#        i. TSS featureCounts file
#        ii. Normalized TSS featureCounts file


bamsortrmdup=$1
assay=$2
genome=$3
scale10mil=$4


# Determine read extension.
if [ "$assay" == "ChIP" ]; then
    read_ext=200
elif [ "$assay" == "ATAC" ]; then
    read_ext=0
elif [ "$assay" == "RNA" ]; then
    read_ext=0
else
    echo -e "Unrecognized assay name (choose 'ChIP'| 'ATAC' | 'RNA): $assay"
    exit
fi


### Build bigWig.
echo -e "\n============================================="
echo -e "RUNNING: bigWig construction"
echo "Current time : $(date)"

chromsizes=$ANNOTATIONS/$genome/$genome.chromInfo.txt
if [ ! -e $chromsizes ]; then
    echo -e "Make sure RESOURCES variable points to annotation files and folders."
    exit
fi

echo -e "Scaling factor (to 10M reads): ${scale10mil}\n"
bwout=${bamsortrmdup/%.bam/.10mNorm.bw}
${CER_CODE}/CER_scripts/pipeline/make_bigwig.sh $bamsortrmdup $chromsizes $read_ext $scale10mil $bwout

# Count fragments by TSS.
${CER_CODE}/CER_scripts/pipeline/run_tss.sh ${bamsortrmdup} ${genome} ${scale10mil}

echo -e "Finished normalization : $(date)\n"
