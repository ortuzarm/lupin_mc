#!/bin/bash
# scripts for 16S rRNA profiling data 
# demultiplexing reads by the alignment of primer sequences
# Author garridoo@mpipz.mpg.de (https://github.com/garridoo/ljsphere), modified from maiteortuzar@usal.es

source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate qiime1

# Demultiplexing
DATA_DIR=/home/ortuzar
R1=$DATA_DIR/SC_Undetermined_S0_L001_R1_001.fastq.gz
R2=$DATA_DIR/SC_Undetermined_S0_L001_R2_001.fastq.gz
I1=$DATA_DIR/SC_Undetermined_S0_L001_I1_001.fastq.gz
I2=$DATA_DIR/SC_Undetermined_S0_L001_I2_001.fastq.gz
MAPPING_FILE=$DATA_DIR/SynCom_mapping_barcodes_2025.tsv

DEMUX_OUT=$DATA_DIR/SC_demultiplexed_samples_$(date +%Y%m%d_%H%M%S)
mkdir -p $DEMUX_OUT

BARCODES_FASTQ=$DEMUX_OUT/barcodes.fastq.gz
echo
paste <(zcat $I1) <(zcat $I2) | \
awk 'NR%4==1 {header=$1} NR%4==2 {seq=$1$2; print header"\n"seq"\n+\n"gensub(/./,"I","g",seq)}' | \
gzip > $BARCODES_FASTQ

echo 
BARCODE_LEN=24

split_libraries_fastq.py \
    -i $R1 \
    -b $BARCODES_FASTQ \
    -m $MAPPING_FILE \
    --rev_comp_mapping_barcodes \
    --barcode_type $BARCODE_LEN \
    --store_demultiplexed_fastq \
    -o $DEMUX_OUT/forward_temp

split_libraries_fastq.py \
    -i $R2 \
    -b $BARCODES_FASTQ \
    -m $MAPPING_FILE \
    --rev_comp_mapping_barcodes \
    --barcode_type $BARCODE_LEN \
    --store_demultiplexed_fastq \
    -o $DEMUX_OUT/reverse_temp

echo 
FINAL_OUT=$DEMUX_OUT/final_samples
mkdir -p $FINAL_OUT

if [[ -f $DEMUX_OUT/forward_temp/seqs.fastq ]]; then
    gzip -c $DEMUX_OUT/forward_temp/seqs.fastq > $FINAL_OUT/forward.fastq.gz
else
    echo 
    exit 1
fi

if [[ -f $DEMUX_OUT/reverse_temp/seqs.fastq ]]; then
    gzip -c $DEMUX_OUT/reverse_temp/seqs.fastq > $FINAL_OUT/reverse.fastq.gz
else
    echo 
    exit 1
fi

echo 


# Organize by SampleID
echo 
OUT_BASE=$DATA_DIR/SC_demultiplexed_by_sample_$(date +%Y%m%d_%H%M%S)
mkdir -p $OUT_BASE

FORWARD_FASTQ=$OUT_BASE/forward.fastq
REVERSE_FASTQ=$OUT_BASE/reverse.fastq

gunzip -c $FINAL_OUT/forward.fastq.gz > $FORWARD_FASTQ
gunzip -c $FINAL_OUT/reverse.fastq.gz > $REVERSE_FASTQ

FORWARD_SPLIT=$OUT_BASE/forward_split
mkdir -p $FORWARD_SPLIT
split_sequence_file_on_sample_ids.py --file_type fastq -i $FORWARD_FASTQ -o $FORWARD_SPLIT

REVERSE_SPLIT=$OUT_BASE/reverse_split
mkdir -p $REVERSE_SPLIT
split_sequence_file_on_sample_ids.py --file_type fastq -i $REVERSE_FASTQ -o $REVERSE_SPLIT

FINAL_BY_SAMPLE=$OUT_BASE/samples_per_id
mkdir -p $FINAL_BY_SAMPLE

for fwd_file in $FORWARD_SPLIT/*.fastq; do
    SAMPLE=$(basename $fwd_file .fastq)
    SAMPLE_DIR=$FINAL_BY_SAMPLE/$SAMPLE
    mkdir -p $SAMPLE_DIR
    cp $fwd_file $SAMPLE_DIR/${SAMPLE}_R1.fastq
done

for rev_file in $REVERSE_SPLIT/*.fastq; do
    SAMPLE=$(basename $rev_file .fastq)
    SAMPLE_DIR=$FINAL_BY_SAMPLE/$SAMPLE
    mkdir -p $SAMPLE_DIR
    cp $rev_file $SAMPLE_DIR/${SAMPLE}_R2.fastq
done

for SAMPLE_DIR in $FINAL_BY_SAMPLE/*; do
    gzip -f $SAMPLE_DIR/*.fastq
done


# Export
EXPORT_DIR=$DATA_DIR/SC_export_$(date +%Y%m%d_%H%M%S)
mkdir -p $EXPORT_DIR

echo 
cp -v $FINAL_BY_SAMPLE/*/*.fastq.gz $EXPORT_DIR/
OUT_FILE=$EXPORT_DIR/samples_list.tsv
echo -e "SampleID\tForward\tReverse" > $OUT_FILE

for fwd in $EXPORT_DIR/*_R1.fastq.gz; do
    sample=$(basename "$fwd" | sed 's/_R1.fastq.gz//')
    rev="$EXPORT_DIR/${sample}_R2.fastq.gz"
    if [[ -f "$rev" ]]; then
        echo -e "${sample}\t${fwd}\t${rev}" >> $OUT_FILE
    else
        echo 
    fi
done

echo 
