#!/bin/bash

# Modified from MSKCC CER
# File: merge_and_filter_peaks.sh
#
# Usage:
#    merge_and_filter_peaks.sh -B <blacklist.bed> -D <merge-dist> -O </path/to/outfolder> [peakfile1 peakfile2 peakfile3 ...]


# Get CER shell tools.
if [ ! -e ${NGSTK}/shell_tools.sh ]; then
    echo "Set NGSTK environment variable to path to CER-NGStk code repository."
    exit
fi
. ${NGSTK}/shell_tools.sh


# Define and parse command-line interface.
PARAMS=""
while (( "$#" )); do
    case "$1" in
        -B|--blacklist)
            poor_map_regions=$2
            shift 2
            ;;
        -D|--distance)
            distance=$2
            shift 2
            ;;
        -O|--outdir)
            outdir=$2
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


# Validate command-line inputs.

if [ ! "$distance" ]; then
    echo -e "Must specify distance for merge with -D or --distance"
    exit
fi

if [ ! "$poor_map_regions" ]; then
    echo -e "Must specify blacklist/mask regions with -B or --blacklist"
    exit
elif [ ! -f "$poor_map_regions" ]; then
    echo -e "Alleged blacklist/mask regions isn't a file: ${poor_map_regions}"
    exit
fi

if [ ! "$outdir" ]; then
    echo -e "Must specify output folder with -O or --outdir"
    exit
elif [ ! -d "$outdir" ]; then
    echo -e "Output folder isn't a directory: ${outdir}"
    exit
fi


# Build the filepath and then create the result (BED).
merge_base=$(merged_peaks_path $outdir $distance)
merged_peaks_file=${merge_base}.bed
cat ${PARAMS[@]} | grep -v "rand\|chrUn" | cut -f1-3 | sort -k1,1 -k2,2n | bedtools intersect -v -a stdin -b ${poor_map_regions} | bedtools merge -d ${distance} -i stdin > ${merged_peaks_file}

# Annotations file in SAF format (featureCounts and HOMER, e.g. findMotifsGenome.pl)
annsFile=${merge_base}_SAF.txt
cat ${merged_peaks_file} | awk -v OFS='\t' ' BEGIN{print "GeneID\tChr\tStart\tEnd\tStrand"}; { posID = $1"_"$2"_"$3; print posID, $0, "+" }' > ${annsFile}
