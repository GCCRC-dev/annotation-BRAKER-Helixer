#!/usr/bin/env bash

set -euo pipefail

# change DB paths
#trembl='~/db/uniprot_trembl_plants.dmnd'
#swissprot='~/db/swissprot_plants.dmnd'
#ptrep='~/db/trep-db_proteins_Rel-25.dmnd'


if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <prefix> <query_proteins.fa(.gz)>" >&2
    echo >&2
    echo "Arguments:" >&2
    echo "  prefix              Prefix for all output files (e.g. braker or helixer.u)" >&2
    echo "  query_proteins.fa   Protein FASTA used as DIAMOND query (plain or .gz)" >&2
    echo >&2
    echo "Example:" >&2
    echo "  $0 braker braker.longest.aa.gz" >&2
    echo "  $0 helixer.u helixer.u.faa.gz" >&2
    exit 1
fi

prefix=$1
query=$2
threads=16

fields=(
    qseqid sseqid pident length qlen slen
    qstart qend sstart send
    qcovhsp scovhsp evalue bitscore stitle
)

header='qseqid	sseqid	pident	length	qlen	slen	qstart	qend	sstart	send	qcovhsp	scovhsp	evalue	bitscore	stitle'


run_diamond() {
    local database="$1"
    local output_prefix="$2"

    diamond blastp \
        --query "$query" \
        --db "$database" \
        --out "${output_prefix}.raw.tsv" \
        --outfmt 6 "${fields[@]}" \
        --evalue 1e-10 \
        --max-target-seqs 5 \
        --sensitive \
        --threads "$threads"

    {
        printf '%s\n' "$header"
        cat "${output_prefix}.raw.tsv"
    } > "${output_prefix}.tsv"
}


# Run searches
run_diamond $swissprot ${prefix}.swissprot
run_diamond $trembl ${prefix}.trembl
run_diamond $ptrep ${prefix}.ptrep


# Swiss-Prot:
# First retain hits with >=80% query and subject coverage,
# then choose the highest-bitscore qualifying hit per query.
awk -F'\t' '
BEGIN { OFS = FS }
$11 >= 80 && $12 >= 80 && $13 <= 1e-10
' ${prefix}.swissprot.raw.tsv |
sort -t $'\t' -k1,1 -k14,14gr -k13,13g |
awk -F'\t' '!seen[$1]++' \
> ${prefix}.swissprot.strong.raw.tsv

{
    printf '%s\n' "$header"
    cat ${prefix}.swissprot.strong.raw.tsv
} > ${prefix}.swissprot.strong.tsv


# TrEMBL:
# Same >=80% query and subject coverage criteria.
awk -F'\t' '
BEGIN { OFS = FS }
$11 >= 80 && $12 >= 80 && $13 <= 1e-10
' ${prefix}.trembl.raw.tsv |
sort -t $'\t' -k1,1 -k14,14gr -k13,13g |
awk -F'\t' '!seen[$1]++' \
> ${prefix}.trembl.strong.raw.tsv

{
    printf '%s\n' "$header"
    cat ${prefix}.trembl.strong.raw.tsv
} > ${prefix}.trembl.strong.tsv


# PTREP:
# Require >=70% coverage of the predicted protein.
awk -F'\t' '
BEGIN { OFS = FS }
$11 >= 70 && $13 <= 1e-10
' ${prefix}.ptrep.raw.tsv |
sort -t $'\t' -k1,1 -k14,14gr -k13,13g |
awk -F'\t' '!seen[$1]++' \
> ${prefix}.ptrep.strong.raw.tsv

{
    printf '%s\n' "$header"
    cat ${prefix}.ptrep.strong.raw.tsv
} > ${prefix}.ptrep.strong.tsv
