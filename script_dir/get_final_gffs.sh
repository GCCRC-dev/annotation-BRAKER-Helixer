#!/usr/bin/env bash

set -euo pipefail

# Usage

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <genome.fasta> <braker.gff3> <helixer.u.gff3>" >&2
    echo >&2
    echo "Example:" >&2
    echo "  $0 ../helixer/Vtub.hap1.v1.softmasked.fasta \\" >&2
    echo "     ../braker/braker.gff3 \\" >&2
    echo "     ../helixer/helixer.u.gff3" >&2
    exit 1
fi

#source "$(conda info --base)/etc/profile.d/conda.sh"

genome="$1"
braker_gff="$2"
helixer_gff="$3"
script_dir="$4"

threads=24

# Extract prefix from genome filename
#
# Vtub.hap1.v1.softmasked.fasta
# -> Vtub.hap1.v1

filename=$(basename "$genome")
prefix=$(echo "$filename" | cut -d'.' -f1-3)

echo "Genome: $genome"
echo "Prefix: $prefix"
echo

# Input checks
for f in "$genome" "$braker_gff" "$helixer_gff"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: file not found: $f" >&2
        exit 1
    fi
done

# Evidence ID lists

braker_hc='../braker/braker.HC.gene_ids.txt'
braker_lc='../braker/braker.LC.gene_ids.txt'

helixer_hc='../helixer/helixer.u.HC.ids.txt'
helixer_lc='../helixer/helixer.u.LC.ids.txt'

for f in "$braker_hc" "$braker_lc" "$helixer_hc" "$helixer_lc"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: ID list not found: $f" >&2
        exit 1
    fi
done


# Output filenames

merged="${prefix}.annotation.merged.gff3"
completeC="${prefix}.annotation.completeCDS.gff3"
renamed="${prefix}.annotation.renamed.gff3"

hc_ids="${prefix}.HC.ids.txt"
lc_ids="${prefix}.LC.ids.txt"

hc_gff="${prefix}.HC.gff3"
lc_gff="${prefix}.LC.gff3"

hc_long="${prefix}.HC.longestIsoform.gff3"
lc_long="${prefix}.LC.longestIsoform.gff3"

hc_faa="${prefix}.HC.longestIsoform.faa"
lc_faa="${prefix}.LC.longestIsoform.faa"

validation_dir="${prefix}.validation"

mkdir -p "$validation_dir"


# 1. Merge BRAKER + Helixer unique annotations
# Standardize transcript -> mRNA

echo "### 1. Merging BRAKER and Helixer"

{
    echo '##gff-version 3'
    for gff in "$braker_gff" "$helixer_gff"; do

        awk '
        BEGIN {
            FS=OFS="\t"
        }

        /^#/ {
            next
        }

        {
            if ($3 == "transcript")
                $3 = "mRNA"

            print
        }
        ' "$gff"

    done

} > "$merged"


# ============================================================
# 2. Check / remove incomplete CDS models
# ============================================================

echo "### 2. Checking CDS completeness"

agat_sp_filter_incomplete_gene_coding_models.pl \
    --gff "$merged" \
    --cpu "$threads" \
    --fasta "$genome" \
    --output "$completeC" \
    2>&1 |
    tee "$validation_dir/incomplete_CDS.log"

grep "Number of gene affected" \
    "$validation_dir/incomplete_CDS.log" \
    > "$validation_dir/number_of_incomplete_CDS_removed.txt" \
    || true


# ============================================================
# 3. AGAT statistics for complete merged annotation
# ============================================================

echo "### 3. AGAT statistics"

agat_sp_statistics.pl \
    --gff "$completeC" \
    -o "$validation_dir/${prefix}.agat_statistics.txt"

#################
# RENAME GENES.
# get old ids

awk -F'\t' '
$3=="gene" {
    if (match($9,/ID=([^;]+)/,a))
        print $1, $4, $5, a[1]
}
' OFS='\t' "$completeC" > "${prefix}.genes_ordered.tsv"

# get mapping file between old and new ids
awk -v prefix="$prefix" '
BEGIN {
    FS=OFS="\t"
    print "chrom","start","end","old_gene_id","new_gene_id"
}
{
    chr=$1
    start=$2
    end=$3
    old=$4

    if (chr != prev_chr) {
        n=10
        prev_chr=chr
    } else {
        n+=10
    }

    new=sprintf("%s.%s.g%06d", prefix, chr, n)

    print chr,start,end,old,new
}
' "${prefix}.genes_ordered.tsv" \
> "${prefix}.gene_ID_mapping.tsv"

# rename them
perl $script_dir/rename_gff_id_mapping.pl $completeC ${prefix}.gene_ID_mapping.tsv > $renamed


# ============================================================
# 4. Check that AGAT can parse/reconstruct the GFF hierarchy
# ============================================================

echo "### 4. Checking GFF hierarchy"

agat_convert_sp_gxf2gxf.pl \
    --gff "$renamed" \
    -o "$validation_dir/${prefix}.agat_parsed.gff3" \
    > "$validation_dir/agat_parser.log" 2>&1


# ============================================================
# 5. Basic feature counts
# ============================================================

echo "### 5. Counting features"

awk -F'\t' '
BEGIN {
    OFS="\t"
}

$0 !~ /^#/ {

    if ($3=="gene")
        gene++

    if ($3=="mRNA")
        mrna++

    if ($3=="exon")
        exon++

    if ($3=="CDS")
        cds++

    if ($3=="intron")
        intron++

    if ($3=="five_prime_UTR")
        utr5++

    if ($3=="three_prime_UTR")
        utr3++
}

END {
    print "genes", gene+0
    print "mRNAs", mrna+0
    print "exons", exon+0
    print "CDS", cds+0
    print "introns", intron+0
    print "five_prime_UTR", utr5+0
    print "three_prime_UTR", utr3+0
}
' "$renamed" \
> "$validation_dir/feature_counts.tsv"


# ============================================================
# 6. Genes per chromosome
# ============================================================

echo "### 6. Counting genes per chromosome"

awk -F'\t' '
$0 !~ /^#/ && $3=="gene" {
    n[$1]++
}

END {
    for (chr in n)
        print chr "\t" n[chr]
}
' "$renamed" |
sort -V \
> "$validation_dir/genes_per_chromosome.tsv"


# ============================================================
# 7. Same-strand overlapping genes
# ============================================================

echo "### 7. Checking overlapping genes"

awk '
BEGIN {
    FS=OFS="\t"
}

$0 !~ /^#/ && $3=="gene" {

    id="."

    if (match($9, /ID=([^;]+)/, a))
        id=a[1]

    print $1, $4-1, $5, id, ".", $7
}
' "$renamed" |
sort -k1,1 -k2,2n \
> "$validation_dir/genes.bed"


bedtools intersect \
    -a "$validation_dir/genes.bed" \
    -b "$validation_dir/genes.bed" \
    -s \
    -wa -wb |
awk -F'\t' '$4 != $10' \
> "$validation_dir/overlapping_genes.tsv"

# ============================================================
# 8. Generate HC annotation
# ============================================================

echo "### 8. Generating HC annotation"

cat "$braker_hc" "$helixer_hc" | sed -E 's/\.t[0-9]+$//' | sort -u > tmp.hc.ids.txt

awk 'NR==FNR && FNR>1 {map[$4]=$5; next}
     {
         if ($1 in map)
             print map[$1]
         else
		 print "HC removed/not present: " $1 > "/dev/stderr"
     }' \
    "${prefix}.gene_ID_mapping.tsv" tmp.hc.ids.txt \
    > "$hc_ids"

agat_sp_filter_feature_from_keep_list.pl \
    --cpu 12 \
    --gff "$renamed" \
    --keep_list "$hc_ids" \
    --output "$hc_gff"


# ============================================================
# 9. Generate LC annotation
# ============================================================

echo "### 9. Generating LC annotation"

cat "$braker_lc" "$helixer_lc" | sed -E 's/\.t[0-9]+$//' | sort -u > tmp.lc.ids.txt

awk 'NR==FNR && FNR>1 {map[$4]=$5; next}
     {
         if ($1 in map)
             print map[$1]
         else
                 print "LC removed/not present: " $1 > "/dev/stderr"
     }' \
    "${prefix}.gene_ID_mapping.tsv" tmp.lc.ids.txt \
    > "$lc_ids"

agat_sp_filter_feature_from_keep_list.pl \
    --cpu 12 \
    --gff "$renamed" \
    --keep_list "$lc_ids" \
    --output "$lc_gff"


# ============================================================
# 10. AGAT stats for HC and LC separately
# ============================================================

echo "### 10. HC / LC statistics"

agat_sp_statistics.pl \
    --gff "$hc_gff" \
    -o "$validation_dir/${prefix}.HC.agat_statistics.txt"

agat_sp_statistics.pl \
    --gff "$lc_gff" \
    -o "$validation_dir/${prefix}.LC.agat_statistics.txt"


# ============================================================
# 11. Longest isoform GFFs
# ============================================================

echo "### 11. Selecting longest isoforms"

agat_sp_keep_longest_isoform.pl \
    --cpu 12 \
    --gff "$hc_gff" \
    -o "$hc_long"

agat_sp_keep_longest_isoform.pl \
    --cpu 12 \
    --gff "$lc_gff" \
    -o "$lc_long"


# ============================================================
# 12. Representative proteins
#
# Temporary/review proteins.
# Final release FASTAs can later be generated with gffread
# after gene renaming.
# ============================================================

echo "### 12. Extracting longest-isoform proteins"

agat_sp_extract_sequences.pl \
    --cpu 12 \
    --gff "$hc_long" \
    --fasta "$genome" \
    --protein \
    -o "$hc_faa"

agat_sp_extract_sequences.pl \
    --cpu 12 \
    --gff "$lc_long" \
    --fasta "$genome" \
    --protein \
    -o "$lc_faa"


# ============================================================
# 13. Protein length summaries
# ============================================================

echo "### 13. Protein length statistics"

seqkit fx2tab -n -i -l "$hc_faa" \
    > "$validation_dir/HC.protein_lengths.tsv"

seqkit fx2tab -n -i -l "$lc_faa" \
    > "$validation_dir/LC.protein_lengths.tsv"


# overwrite, do not append from previous runs
: > "$validation_dir/protein_size_summary.txt"


for class in HC LC; do

    if [[ "$class" == "HC" ]]; then
        sizes="$validation_dir/HC.protein_lengths.tsv"
    else
        sizes="$validation_dir/LC.protein_lengths.tsv"
    fi

    echo "$class" >> "$validation_dir/protein_size_summary.txt"

    awk '
    {
        n++

        sum += $2

        if ($2 < 50)
            a++

        else if ($2 < 100)
            b++

        else
            c++
    }

    END {
        print "proteins:", n+0
        print "<50 aa:", a+0
        print "50-99 aa:", b+0
        print ">=100 aa:", c+0

        if (n > 0)
            print "mean length:", sum/n

        print ""
    }
    ' "$sizes" \
    >> "$validation_dir/protein_size_summary.txt"

done


# ============================================================
# 14. Compact validation summary
# ============================================================

echo "### 14. Writing summary"

{
    echo "========================================"
    echo "ANNOTATION SUMMARY: $prefix"
    echo "========================================"
    echo

    echo "### CDS COMPLETENESS"
    cat "$validation_dir/number_of_incomplete_CDS_removed.txt" 2>/dev/null \
        || echo "No incomplete-CDS count found"
    echo

    echo "### FINAL FEATURE COUNTS"
    cat "$validation_dir/feature_counts.tsv"
    echo

    echo "### HC / LC GENE COUNTS"

    printf "HC genes\t"
    awk -F'\t' '$3=="gene"{n++} END{print n+0}' "$hc_gff"

    printf "LC genes\t"
    awk -F'\t' '$3=="gene"{n++} END{print n+0}' "$lc_gff"

    echo
    echo "### PROTEIN LENGTHS"
    cat "$validation_dir/protein_size_summary.txt"

    echo
    echo "### GENES PER CHROMOSOME"
    cat "$validation_dir/genes_per_chromosome.tsv"

} > "$validation_dir/validation_summary.txt"


# generate final fastas.
conda activate rnaseq
gffread "$hc_gff" -g "$genome" -w ${prefix}.HC.transcripts.fa
gffread "$hc_gff" -g "$genome" -x ${prefix}.HC.cds.fa
gffread "$hc_gff" -g "$genome" -y ${prefix}.HC.proteins.fa

gffread "$lc_gff" -g "$genome" -w ${prefix}.LC.transcripts.fa
gffread "$lc_gff" -g "$genome" -x ${prefix}.LC.cds.fa
gffread "$lc_gff" -g "$genome" -y ${prefix}.LC.proteins.fa

echo
echo "========================================"
echo "DONE"
echo "========================================"
echo
echo "Merged annotation:"
echo "  $completeC"
echo
echo "HC annotation:"
echo "  $hc_gff"
echo
echo "LC annotation:"
echo "  $lc_gff"
echo
echo "Longest isoforms:"
echo "  $hc_long"
echo "  $lc_long"
echo
echo "Transcript, CDS and protein FASTAs ready."
echo
echo "Validation summary:"
echo "  $validation_dir/validation_summary.txt"
