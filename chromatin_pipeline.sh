#!/bin/bash

# Pipeline for a single ChIP- or ATAC-seq sample, single- or paired-end.
#
# Usage: /path/to/this/script <assay> <readtype> <genome> <read1> <read2>
#
# Arguments
#     1. assay -- the name of the experiment type ('ATAC' or 'ChIP')
#     2. readtype -- indication of sequencing protocol type ('PE' for paired-end, 'SE' for single-end)
#     3. genome -- name of the genomic assembly ('mm10' or 'hg19')
#     4. read1 -- path to read1 file
#     5. read2 -- path to read2 file (can use arbitrary dummy value for single-end as this is ignored then)
#     6. Output directory


# Check if key environment variables are set. If not, set to default.

if [ -z ${NGSTK} ]; then
    NGSTK=/data/cer/shared/code/CER-NGStk
    echo "NGSTK environment variable empty, setting to ${NGSTK}"
fi

if [ -z ${CER} ]; then
    CER=/data/cer/shared
    echo "CER environment variable empty, setting to ${CER}"
fi

if [ -z ${BT2_IDX} ]; then
    BT2_IDX=/data/cer/shared/resources/Genomes_BT2
    echo "BT2_IDX environment variable empty, setting to ${BT2_IDX}"
fi

if [ -z ${ANNOTATIONS} ]; then
    ANNOTATIONS=/data/cer/shared/resources/annotations
    echo "ANNOTATIONS environment variable empty, setting to ${ANNOTATIONS}"
fi

if [ -z ${PICARD} ]; then
    PICARD=/home/kocher/Programs/picard/build/libs/picard.jar
    echo "PICARD environment variable empty, setting to ${PICARD}"
fi

if [ -z ${CER_CODE} ]; then
    CER_CODE=/data/cer/shared/code
    echo "CER_CODE environment variable empty, setting to ${CER_CODE}"
fi



# Get CER shell tools.
if [ ! -e ${NGSTK}/shell_tools.sh ]; then
    echo "Set NGSTK environment variable to path to CER-NGStk code repository."
    exit
fi
. ${NGSTK}/shell_tools.sh

## Usually run from python/perl wrapper with system call, a la:
## 
## wrapper run example:
## bsub -n4 -W 08:00 -R "rusage[mem=10]" -o %J.stdout_RNA_STAR_hg19PE_fullPL_outfile.txt ./run_AlignRNA_STAR_hg19_PE_HTseq_bwEtc_lilac_2017-12-10.sh test_hs_R1_b2min.fastq.gz  test_hs_R2_b2min.fastq.gz

numthreads=6

# module dependency below: R
module add R/R-3.5.1 

# use CER python virtual environment
source /data/cer/shared/venvs/cer_base/bin/activate

# some path discrepancies on computenodes, be sure to include:
#source ~/.bash_profile

#personal + HPC installs
#trimG=/home/kocher/Programs/TrimGalore-0.4.5/trim_galore
#cutadapt=/home/kocher/.local/bin/cutadapt
#cutadapt=/data/cer/shared/venvs/cer_base/bin/cutadapt
#bowtie2=/home/kocher/Programs/bowtie2-2.3.4-linux-x86_64/bowtie2
#samtools=/programs/x86_64-linux/samtools/1.4.1/bin/samtools
#picard=/home/kocher/Programs/picard/build/libs/picard.jar


echo "=========================================================="
echo "Starting on : $(date)"
echo "=========================================================="

### saba pipeline ported to lilac
### PE RNA-seq for mm9
### two fastq.gz inputs specified

### usually run with 6 threads/10 Gb each

PARAMS=""
while (( "$#" )); do
    case "$1" in
        --assay)
            assay=$2
            shift
            ;;
        --readtype)
            readtype=$2
            shift
            ;;
        --genome)
            genome=$2
            shift
            ;;
        --rone0)
            rone0=$2
            shift
            ;;
        --rtwo0)
            rtwo0=$2
            shift 2
            ;;
        --outDir)
            outDir=$2
            shift 2
            ;;
        --spikein)
            spikein=$2
            shift 2
            ;;
        --stranded)
            stranded=$2
            shift 2
            ;;
        --numthreads)
            numthreads=$2
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



echo -e "assay = $assay\nreadtype = $readtype\ngenome = $genome\nread1 = $rone0\nread2 = $rtwo0\noutDir = $outDir\nspikein = $spikein"
#exit
# Parameters
#assay=$1
#readtype=$2
#genome=$3
#
## INPUT = *.fastq.gz
#rone0=$4
#rtwo0=$5
#
##Output
#outDir=$6

outPath=$(readlink -f ${outDir})
fastqPath=$(readlink -f ./)
# Make sure that current path isn't the same as outDir
if [ "$fastqPath" == "$outPath" ]; then
        echo -e "fastq folder is the same as output. Specify other output folder"
        exit
fi

baseR1=$(basename ${rone0})
ln -s ${rone0} ${outDir}/${baseR1} 
rone0="${outDir}/${baseR1}"


# Ensure that we have a valid read type scenario.
if [ "$readtype" == "PE" ]; then
	baseR2=$(basename ${rtwo0})
	ln -s ${rtwo0} ${outDir}/${baseR2}
    if [ ! -e "${rtwo0}" ]; then
        echo -e "Selected paired-end but R2 doesn't exist: ${rtwo0}"
        exit
    fi
	rtwo0="${outDir}/${baseR2}"
elif [ "$readtype" != "SE" ]; then
    echo -e "Invalid read type (choose 'SE' or 'PE'): ${readtype}"
    exit
fi


# Jump to outDir for processing
cd ${outDir}



# Start with a clean sample stats file, so that we can append throught the duration of new processing without getting duplication.
if [ -e stats.tsv ]; then
    echo -e "Removing existing samples stats file."
    rm stats.tsv
fi


# Validate choice of assay / experiment type.
if [ "$assay" == "ChIP" ]; then
    read_ext=200
elif [ "$assay" == "ATAC" ]; then
    read_ext=0
else
    echo -e "Unrecognized assay name (choose 'ChIP' or 'ATAC'): $assay"
    exit
fi


#Name sample
sampleName=$(basename ${outDir}| sed 's/Sample_//g')


# Parse the genome and determine its variable setting implications.
human_bg_data="${CER}/resources/BackgroundData_ChIPinput/humanCD34_input.$genome.sorted.RmDup.bam" # TODO Change to ${BACKGROUND} when ready
mouse_bg_data="${CER}/resources/BackgroundData_ChIPinput/mouseESC_input.$genome.sorted.RmDup.bam"
if [[ "$spikein" == 1 ]]; then
    bt2_idx="/data/cer/shared/resources/Genome_BT2_Drosophila/$genome"
else
    bt2_idx="${BT2_IDX}/$genome/$genome"
fi

if [[ "$genome" == "hg19" || "$genome" == "hg38" ]]; then
    organism=hs
    bg_data=$human_bg_data
elif [[ "$genome" == "mm9" || "$genome" == "mm10" ]]; then
    organism=mm
    bg_data=$mouse_bg_data
else
    echo -e "Invalid genome choice: $genome; only mm9, mm10, hg19, and hg38 supported"
    exit
fi

echo $bt2_idx


echo -e "\n============================================="
echo -e "RUNNING: trimming and alignment"
echo "Current time : $(date)"


# clean up path on input files and name based on trim_galore output:
rone0nodir=${rone0##*/}
rtwo0nodir=${rtwo0##*/}



#/home/kocher/Programs/TrimGalore-0.4.5/trim_galore --quality 15 --fastqc --gzip --path_to_cutadapt /home/kocher/.local/bin/cutadapt --paired $2 $3
if [ ! -e "${NGSTK}/run_trim_align.py" ]; then
    echo -e "Cannot find script for trimming and alignment; NGSTK=${NGSTK}"
    exit
fi

if [ "$readtype" == "PE" ]; then
    rone1=${rone0nodir/%.fastq.gz/_val_1.fq.gz}
    rtwo1=${rtwo0nodir/%.fastq.gz/_val_2.fq.gz}
    alignout=${sampleName}.$genome.bam
    ${NGSTK}/run_trim_align.py --reads-files $rone0 $rtwo0 --threads $numthreads --trim-qual 15 --bt2-idx $bt2_idx --output $outDir/$alignout
elif [ "$readtype" == "SE" ]; then
    rone1=${rone0nodir/%.fastq.gz/_trimmed.fq.gz}
    alignout=${sampleName}.$genome.bam
    ${NGSTK}/run_trim_align.py --readsfile $rone0 --threads $numthreads --trim-qual 15 --bt2-idx $bt2_idx --output $outDir/$alignout
else
    echo -e "Invalid read type: ${readtype}"
    exit
fi

if [ ! -e "$alignout" ]; then
    echo -e "No alignment file: ${alignout}"
    exit
fi
alignment_filesize=$(stat -c%s ${alignout})
if [ "0" == "$alignment_filesize" ]; then
    echo -e "Empty alignment file: ${alignment_filesize}"
    exit
fi

R1reads=$(grep "Total reads processed" *_R1_*.gz_trimming_report.txt|  awk '{print $NF}'|sed 's/,//g')


echo -e "\n============================================="
echo -e "CLEARING: Intermediate files ${rone1}, ${rtwo1}"
echo "Current time : $(date)"
rm $outDir/${rone1}
if [ "$readtype" == "PE" ]; then
    rm $outDir/${rtwo1}
fi


#trim_galore --quality 15 --fastqc --gzip --path_to_cutadapt $cutadapt --paired $rone0 $rtwo0
#echo -e "\n============================================="
#echo -e "RUNNING: bowtie2 alignment to ${genome}"
#echo "Current time : $(date)"
#bowtie2 --local -p $numthreads -x $bt2_idx -1 $rone1 -2 $rtwo1 | samtools view -@ $numthreads -bS - > $alignout
if [[ "$spikein" == 1 ]]; then
    echo -e "\n============================================="
    echo -e "RUNNING: Drosophila spike-in filter and counting"
    echo "Current time : $(date)"
    samtools sort -@ $numthreads -o ${genome}.tmp.sorted.bam $alignout
    samtools index ${genome}.tmp.sorted.bam
    
    genomeChrom=$(cat $ANNOTATIONS/${genome}/${genome}.chromInfo.txt |awk '{print $1}'|grep -v _random| tr '\n' ' ')
    Dm6Chrom=$(cat /data/cer/shared/resources/Genome_BT2_Drosophila/dm6.chrom.sizes |awk '{print $1}'| tr '\n' ' ')
    
    
    samtools view -H ${genome}.tmp.sorted.bam > ${genome}.header.sam
    
    if [[ "$readtype" == "PE" ]]; then
        samtools view ${genome}.tmp.sorted.bam ${Dm6Chrom}  | awk '$7 !~/^chr/' > dm6.chrom.sam
        cat ${genome}.header.sam dm6.chrom.sam |samtools view -b - -o dm6.tmp.bam
        rm dm6.chrom.sam
    else
        samtools view -b ${genome}.tmp.sorted.bam ${Dm6Chrom}  > dm6.tmp.bam
    fi
    
    
    
    samtools sort  -@ $numthreads -o dm6.tmp.sorted.bam dm6.tmp.bam

    samtools index dm6.tmp.sorted.bam
    java -Xms500m -Xmx6g -jar ${PICARD} MarkDuplicates REMOVE_DUPLICATES=true METRICS_FILE=Dm6dup.txt AS=true INPUT=dm6.tmp.sorted.bam OUTPUT=dm6.rmDup.bam

    dm6reads=$( samtools view  -q 10 -F 0x904 dm6.rmDup.bam| wc -l)
    echo "Drosophila reads (Dm6)= ${dm6reads}"
    echo -e "${dm6reads}\t${sampleName}" > dm6reads.tsv 
    
    
    if [[ "$readtype" == "PE" ]]; then
        samtools view ${genome}.tmp.sorted.bam ${genomeChrom} | awk '$7 !~/^Dm/' > ${genome}.tmpNoSpikeIn.sam
        cat ${genome}.header.sam ${genome}.tmpNoSpikeIn.sam | samtools view -b - -o tmp.bam
        rm ${genome}.tmpNoSpikeIn.sam
    else
        samtools view -b ${genome}.tmp.sorted.bam ${genomeChrom} -o tmp.bam
    fi
    
    
    samtools view -H $alignout |grep -vf   <(awk '{print $1}' /data/cer/shared/resources/Genome_BT2_Drosophila/dm6.chrom.sizes) - > header.sam
    
    samtools reheader -P -i header.sam tmp.bam > $alignout.tmp
    
    mv $alignout.tmp  $alignout
    
    

    rm tmp.bam
    rm dm6.tmp.bam
    rm dm6.tmp.bam.bai
    rm dm6.tmp.sorted.bam
    rm dm6.tmp.sorted.bam.bai
    rm tmp.bam
    rm header.sam
    rm ${genome}.tmp.sorted.bam
    rm dm6.rmDup.bam
fi

echo -e "\n============================================="
echo -e "RUNNING: position-sort bam, index, and get pre-markdup read summary"
echo "Current time : $(date)"

bamsort=${alignout/%.bam/.sorted.bam}

samtools sort -@ $numthreads -o $bamsort $alignout

samtools index $bamsort

samtools idxstats $bamsort

echo -e "\n============================================="
echo -e "CLEARING: Intermediate files $bamsort"
echo "Current time : $(date)"
if [ ! -e "$bamsort" ]; then
    echo -e "No sorted alignment file: $bamsort"
    exit
fi
rm $outDir/${alignout}


mtReads=$(samtools view -c ${bamsort} chrM)


echo -e "\n============================================="
echo -e "RUNNING: Picard mark duplicates -- remove"
echo "Current time : $(date)"

#bamsortrmdup = output from picard rmdup: 
bamsortrmdup=${bamsort/%.bam/.RmDup.bam}

java -Xms500m -Xmx6g -jar ${PICARD} MarkDuplicates REMOVE_DUPLICATES=true METRICS_FILE=dup.txt AS=true INPUT=$bamsort OUTPUT=$bamsortrmdup

echo -e "\n============================================="
echo -e "RUNNING: Alignment file indexing and statistics"
echo "Current time : $(date)"

samtools index $bamsortrmdup
echo -e "\nsamtools idxstats:\n"
samtools idxstats $bamsortrmdup
echo -e "\nsamtools flagstat:\n"
samtools flagstat $bamsortrmdup


echo -e "\n============================================="
echo -e "CLEARING: Intermediate files ${bamsort}"
echo "Current time : $(date)"
if [ ! -e "$bamsortrmdup" ]; then
    echo -e "No sorted alignment file: $bamsortrmdup"
    exit
fi
rm $outDir/${bamsort}
rm $outDir/${bamsort}.bai

echo -e "\n============================================="
echo -e "RUNNING: Picard insert size metrics\n"
echo "Current time : $(date)"

### Get insert size estimates from Picard 
picout=PicardHistForInsertMetrics_${bamsortrmdup}.pdf
picouttxt=PicardHistForInsertMetrics_${bamsortrmdup}.txt
java -Xms500m -Xmx2g -jar ${PICARD} CollectInsertSizeMetrics HISTOGRAM_FILE=$picout INPUT=$bamsortrmdup AS=true OUTPUT=$picouttxt MINIMUM_PCT=0.2



### Count fragments for normalization multiplier.
echo -e "\n============================================="
echo -e "RUNNING: fragment counting\n"
echo "Current time : $(date)"

# Use 10 as the MAPQ threshold for counting fragments.
# Definition of count_frags function comes from a shell utilities file in the
# CER's NGStk (CER-NGStk repository, shell_tools.sh file).
n_frags=$(count_frags $bamsortrmdup 10)
echo -e "Counting with records that are mapped, primary, non-supplementary alignments"
echo -e "Fragments: $n_frags\n"
echo -e "Counting with samtools idxstats"
echo -e "Fragments (assuming single-end):\n"
samtools idxstats $bamsortrmdup | grep -v "chrM" | awk '{s+=$3} END {print s}'



##### Peaks (calling with MACS2)
echo -e "\n============================================="
echo -e "RUNNING: peak calling\n"
echo "Current time : $(date)"

# Clean up path and suffix to leave basename for output.
bambase=$(basename "$bamsortrmdup" .bam)
macsout=MACS2_vsInput_${bambase}

if [ ! -e $bg_data ]; then
    echo -e "Cannot locate input/background data; is BACKGROUND variable accurate?\n"
    echo -e "BACKGROUND: ${BACKGROUND}\n"
    exit
fi

# Call raw peaks, using standard input as the "control/treatment" for comparison.
macs2 callpeak -t $bamsortrmdup -c $bg_data -f BAM -g $organism --nomodel --pvalue 0.001 -n $macsout
n_peaks=$(wc -l ${macsout}_peaks.narrowPeak | cut -f1 -d" ")


# Do the normalization-dependent steps (i.e., bigWig construction and counting of fragments by TSS).
# The script specific for the normalization steps doesn't make an assumption about the working
# directory from which it's run, so pass the current (sample folder) directory along with the 
# name of the sorted, deduplicated aligned reads file when it's called here.
N_norm=10000000
scale10mil=$(awk "BEGIN {print ${N_norm}/${n_frags}}")

${CER_CODE}/CER_scripts/pipeline/run_norm.sh ${bamsortrmdup} ${assay} ${genome} ${scale10mil}

mtReadsProp=$(echo "scale=3; $mtReads/${R1reads}" | bc)
dupReads=$(cat dup.txt | grep -w "Unknown Library" | awk -F '\t' '{print $9}' |  cut -c -4)

samtools view -b -q 10 -F 0x904 -o "${sampleName}.${genome}.RmDup.filtered.bam" $bamsortrmdup 

readsInPeaks=$(bedtools intersect -sorted -c -a ${macsout}_peaks.narrowPeak -b "${sampleName}.${genome}.RmDup.filtered.bam" | awk '{sum+=$NF} END {print sum}')
rm "${sampleName}.${genome}.RmDup.filtered.bam"


#SPOT (Signal Portion of Tag) #Rough estimate by collapsing number of paired ended reads by 2. 
if [ "$readtype" == "PE" ]; then
    readsInPeaks=$(echo "${readsInPeaks}/2"| bc )
fi

SPOT=$(echo "scale=2; $readsInPeaks/${n_frags}" | bc)

propFrags=$(echo "scale=2; ${n_frags}/${R1reads}" | bc)
peaksNearTSS=$(bedops -u --range 5000 ${ANNOTATIONS}/${genome}/gencode/gencode.*.TSS.bed| bedmap --echo --echo-map --skip-unmapped ${macsout}_peaks.narrowPeak - | wc -l)
peaksNearTSSprop=$(echo "scale=2; $peaksNearTSS/${n_peaks}" | bc)

Rscript /data/cer/shared/software/phantompeakqualtools/run_spp_nodups.R -p=${numthreads} -c=${bamsortrmdup} -savp -out=CrossCorr.txt

NSC=$(awk '{print $9}' CrossCorr.txt)
RSC=$(awk '{print $10}' CrossCorr.txt)
QualityTags=$(awk '{print $11}' CrossCorr.txt)
rm $outDir/CrossCorr.txt


echo -e "\n============================================="
echo -e "RUNNING: STATs\n"
echo "Current time : $(date)"
echo -e "reads\t${R1reads}\t${sampleName}" | tee -a  stats.tsv
echo -e "dupReads_Prop\t${dupReads}\t${sampleName}" | tee -a  stats.tsv
echo -e "mtReads\t${mtReadsProp}\t${sampleName}" | tee -a  stats.tsv
echo -e "fragments\t${n_frags}\t${sampleName}" | tee -a  stats.tsv
echo -e "fragments_Prop\t${propFrags}\t${sampleName}" | tee -a  stats.tsv
echo -e "peaks\t${n_peaks}\t${sampleName}" | tee -a  stats.tsv
echo -e "SPOT\t${SPOT}\t${sampleName}" | tee -a  stats.tsv
echo -e "peaksNearTSS\t${peaksNearTSSprop}\t${sampleName}" | tee -a  stats.tsv
echo -e "NSC\t${NSC}\t${sampleName}" | tee -a  stats.tsv
echo -e "RSC\t${RSC}\t${sampleName}" | tee -a  stats.tsv
echo -e "QualityTags\t${QualityTags}\t${sampleName}" | tee -a  stats.tsv
echo -e "scale\t${scale10mil}\t${sampleName}" | tee -a  stats.tsv
if [[ "$spikein" == 1 ]]; then
    echo -e "DmSpikeIn\t${dm6reads}\t${sampleName}"| tee -a  stats.tsv
fi
echo -e "\n============================================="
echo -e "RUNNING: Quality control through ngsplot and cross correlation\n"
echo "Current time : $(date)"
mkdir -p NGSplot
module unload R/R-3.3.3
#TSS
ngs.plot.r -G ${genome} -R tss -C ${bamsortrmdup} -RR 90 -P $numthreads -O NGSplot/NGSplot.TSS -T "${sampleName}_TSS" -L 3000 
#Gene body
ngs.plot.r -G ${genome} -R genebody -C  ${bamsortrmdup} -RR 90 -P $numthreads -O NGSplot/NGSplot.genebody -T "${sampleName}_Genebody" -L 3000 
#MACS2 peaks
#cp ${macsout}_peaks.narrowPeak ${macsout}_peaks.narrowPeak.bed
ngs.plot.r -G ${genome} -R bed -C  ${bamsortrmdup} -RR 90 -P $numthreads -O NGSplot/NGSplot.Peaks -T "${sampleName}_Peaks" -L 3000 -E ${macsout}_peaks.narrowPeak 
#rm ${macsout}_peaks.narrowPeak.bed
#Fantom5 enhancers
ngs.plot.r -G ${genome} -R bed -C  ${bamsortrmdup} -RR 90 -P $numthreads -O NGSplot/NGSplot.fantom5 -T "${sampleName}_Fantom5" -L 3000 -E "${ANNOTATIONS}/${genome}/fantom5/${genome}.permissive_enhancers_phase_1_and_2.bed"
#Cross corr

### Clean up all intermediate files
#echo -e "\n============================================="
#echo -e "RUNNING: Final cleanup\n"
#echo "Current time : $(date)"
#





echo "=========================================================="
echo "Finished on : $(date)"
echo "=========================================================="
