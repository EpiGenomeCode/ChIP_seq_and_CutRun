#!/bin/bash

# 1. Sorted deduplicated BAM
# 2. chromsizes
# 3. Read extension (200 for ChIP, 0 for ATAC)
# 4. scale10mil (scaling factor to 10 million reads)
# 5. output file

bamsortrmdup=$1
chromsizes=$2
read_ext=$3
scale10mil=$4
bwout=$5

bamToBed -split -i $bamsortrmdup | slopBed -i stdin -g $chromsizes -s -l 0 -r $read_ext | grep -v "rand\|Un" | genomeCoverageBed -split -g $chromsizes -i stdin -bg -scale $scale10mil | wigToBigWig stdin $chromsizes $bwout
