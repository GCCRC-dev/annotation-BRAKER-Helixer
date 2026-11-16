#!/usr/bin/env bash

# change type of read!
#

strandedness="RF"

bam=$1
gff=$2
base=$(basename "$bam" .bam)

# get exons from helixer gff
awk 'BEGIN{OFS="\t"} $3=="exon"' $gff > helixer.u.exons.gff3

bedtools coverage \
    -split \
    -a helixer.u.exons.gff3 \
    -b $bam \
    > helixer.u.exons.coverage.tsv

# get introns from helixer gff
agat_sp_add_introns.pl \
    --gff $gff \
    --out helixer.u_wINTRONS.gff3
    
awk '$3=="intron"' helixer.u_wINTRONS.gff3 > helixer.u.introns.gff3

# convert gff3 coords to bed
awk 'BEGIN{OFS="\t"}
$3=="intron" {
    id="."
    if (match($9, /Parent=([^;]+)/, m)) id=m[1]

    print $1, $4-1, $5, id, ".", $7
}' helixer.u.introns.gff3 \
> helixer.u.introns.bed

sort -k1,1V -k2,2n helixer.u.introns.bed > helixer.u.introns_sorted.bed
mv helixer.u.introns_sorted.bed helixer.u.introns.bed


### JUNCTIONS ANALYSIS
# get junctions
regtools junctions extract -m 20 -M 500000 -s $strandedness $bam -o ${base}.junctions.bed

# this file is BED12, so it shows only the skipped introns, let's convert it
# getting the introns from these junctions
awk 'BEGIN{OFS="\t"}
{
    split($11, sizes, ",")
    intron_start = $2 + sizes[1]
    intron_end   = $3 - sizes[2]

    print $1, intron_start, intron_end, $4, $5, $6
}' ${base}.junctions.bed \
> ${base}.introns.bed

sort -k1,1V ${base}.introns.bed > ${base}.introns_sorted.bed
mv ${base}.introns_sorted.bed ${base}.introns.bed

# compare introns from helixer with introns from RNA-seq junctions.
bedtools intersect \
    -a helixer.u.introns.bed \
    -b ${base}.introns.bed \
    -s \
    -f 1.0 \
    -r \
    -wa -wb \
    > helixer.u.introns.exact.tsv


### DONOR/ACCEPTOR analysis
# starts
awk 'BEGIN{OFS="\t"}
{
    intron_id=$1 ":" $2 ":" $3 ":" $6 ":" $4
    print $1,$2,$2+1,intron_id,$5,$6
}' helixer.u.introns.bed > helixer.u.starts.bed

awk 'BEGIN{OFS="\t"} {
    print $1,$2,$2+1,$4,$5,$6
}' ${base}.introns.bed > rnaseq.starts.bed

bedtools intersect \
    -a helixer.u.starts.bed \
    -b rnaseq.starts.bed \
    -s -wa -wb \
> helixer.u.shared_starts.tsv

# ends
awk 'BEGIN{OFS="\t"}
{
    intron_id=$1 ":" $2 ":" $3 ":" $6 ":" $4
    print $1,$3-1,$3,intron_id,$5,$6
}' helixer.u.introns.bed > helixer.u.ends.bed

awk 'BEGIN{OFS="\t"} {
    print $1,$3-1,$3,$4,$5,$6
}' ${base}.introns.bed > rnaseq.ends.bed

# compare starts/ends from helixer GFF to the RNA-seq ones
bedtools intersect \
    -a helixer.u.ends.bed \
    -b rnaseq.ends.bed \
    -s -wa -wb \
> helixer.u.shared_ends.tsv
