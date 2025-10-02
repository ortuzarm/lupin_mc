#!/bin/bash
# maiteortuzar@usal.es
set -e

BASE_DIR="/home/ortuzar/02.amplicons"
WORK_DIR="$BASE_DIR/1.analysis"
RAW_DIR="$BASE_DIR/rawdata"
METADATA="$BASE_DIR/sample-metadata.tsv"
CLASSIFIER="$BASE_DIR/classifier-ezbio-v5v7.qza"
FORWARD_PRIMER="AACMGGATTAGATACCCKG"
REVERSE_PRIMER="ACGTCATCCCCACCTTCC"

eval "$(conda shell.bash hook)"
conda activate qiime2-2023.2

qiime taxa filter-table \
    --i-table "$BASE_DIR/table-dada2.qza" \
    --i-taxonomy "$BASE_DIR/taxonomy_silva1381.qza" \
    --o-filtered-table "$BASE_DIR/rep-seq-taxonomy-filter.qza" \
    --p-exclude Mitochondria,Chloroplast,Unassigned,Eukaryota

qiime taxa filter-seqs \
    --i-sequences "$BASE_DIR/rep-seqs-dada2.qza" \
    --i-taxonomy "$BASE_DIR/taxonomy_silva1381.qza" \
    --o-filtered-sequences "$BASE_DIR/rep-seqs-dada2-filter.qza" \
    --p-exclude Mitochondria,Chloroplast,Unassigned,Eukaryota

mkdir -p "$BASE_DIR/2.results" "$BASE_DIR/3.export"

qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences "$BASE_DIR/rep-seqs-dada2-filter.qza" \
    --o-alignment "$BASE_DIR/2.results/rep-seqs-dada2-filter-align.qza" \
    --o-masked-alignment "$BASE_DIR/2.results/rep-seqs-dada2-masked.qza" \
    --o-tree "$BASE_DIR/2.results/unrooted-tree.qza" \
    --o-rooted-tree "$BASE_DIR/2.results/rooted-tree.qza"

qiime tools export \
    --input-path "$BASE_DIR/rep-seq-taxonomy-filter.qza" \
    --output-path "$BASE_DIR/temp_table_export"

min_reads=$(biom summarize-table -i "$BASE_DIR/temp_table_export/feature-table.biom" | awk '/Min:/ {print int($2)}')
sampling_depth=$((min_reads * 90 / 100))
[[ $sampling_depth -le 0 ]] && sampling_depth=$min_reads

qiime feature-table rarefy \
    --i-table "$BASE_DIR/rep-seq-taxonomy-filter.qza" \
    --p-sampling-depth $sampling_depth \
    --o-rarefied-table "$BASE_DIR/2.results/rep-seq-taxonomy-filter_rarefy.qza"

rm -r "$BASE_DIR/temp_table_export"

qiime diversity core-metrics-phylogenetic \
    --i-table "$BASE_DIR/2.results/rep-seq-taxonomy-filter_rarefy.qza" \
    --i-phylogeny "$BASE_DIR/2.results/rooted-tree.qza" \
    --m-metadata-file "$METADATA" \
    --p-sampling-depth $sampling_depth \
    --output-dir "$BASE_DIR/2.results/core-metrics-results"

qiime tools export \
    --input-path "$BASE_DIR/table-dada2.qza" \
    --output-path "$BASE_DIR/3.export/table"

biom convert \
    -i "$BASE_DIR/3.export/table/feature-table.biom" \
    -o "$BASE_DIR/3.export/feature-table.tsv" \
    --to-tsv

qiime tools export \
    --input-path "$BASE_DIR/taxonomy_silva1381.qza" \
    --output-path "$BASE_DIR/3.export/"

mv "$BASE_DIR/3.export/taxonomy.tsv" "$BASE_DIR/3.export/taxonomy.tsv" 2>/dev/null || true

qiime tools export \
    --input-path "$BASE_DIR/2.results/unrooted-tree.qza" \
    --output-path "$BASE_DIR/3.export/tree"

mv "$BASE_DIR/3.export/tree/tree.nwk" "$BASE_DIR/3.export/unrooted-tree.nwk"

cp "$METADATA" "$BASE_DIR/3.export/sample-metadata.tsv"

for level in {1..7}; do
    qiime taxa collapse \
        --i-table "$BASE_DIR/table-dada2.qza" \
        --i-taxonomy "$BASE_DIR/taxonomy_silva1381.qza" \
        --p-level $level \
        --o-collapsed-table "$BASE_DIR/table-L${level}.qza"

    qiime tools export \
        --input-path "$BASE_DIR/table-L${level}.qza" \
        --output-path "$BASE_DIR/temp_export_L${level}"

    biom convert \
        -i "$BASE_DIR/temp_export_L${level}/feature-table.biom" \
        -o "$BASE_DIR/3.export/feature-table-L${level}.tsv" \
        --to-tsv

    rm -r "$BASE_DIR/temp_export_L${level}" "$BASE_DIR/table-L${level}.qza"
done

echo 
