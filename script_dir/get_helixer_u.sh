#!/usr/bin/env bash

genome=$1
gff=$2

gffcompare \
    -r braker.gff3 \
    -o helixer_vs_braker \
    $gff 

cut -f3 helixer_vs_braker.${gff}.tmap | tail -n +2 | sort | uniq -c | sort -nr > gff_compare_stats.txt

awk '$3=="u"' helixer_vs_braker.${gff}.tmap > helixer_u.tsv
awk '{print $4}' helixer_u.tsv | sed 's/\.t[1-9]//' > ids_helixer_u_gene.txt

# get gff3 helixer.u
agat_sp_filter_feature_from_keep_list.pl \
    --gff $gff \
    --keep_list ids_helixer_u_gene.txt \
    --output helixer.u.gff3

agat_sp_statistics.pl --gff helixer.u.gff3 -o helixer.u.stats.txt

agat_sp_extract_sequences.pl \
    --gff helixer.u.gff3 \
    --fasta $genome \
    --protein \
    -o helixer.u.faa

