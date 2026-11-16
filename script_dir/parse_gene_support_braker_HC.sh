#!/usr/bin/env bash

input="gene_support.tsv"
#suffix=$1
suffix="HC_candidates"
output="out_gene_support_${suffix}.out"

if [[ -z "$input" ]]; then
    echo "Usage: $0 <evidence.tsv> [output.tsv]" >&2
    exit 1
fi

awk '
BEGIN {
    FS = OFS = "\t"
}

# Keep comment lines
#/^#/ {
#    print
#    next
#}

# Keep the non-comment header
$1 == "gene_id" {
    print
    next
}

{
    n_introns = $5 + 0

    intron_rna  = ($7 == "NA" ? 0 : $7 + 0)       #$7=pct_introns_sup_rnaseq
    intron_prot = ($9 == "NA" ? 0 : $9 + 0)   #$9=pct_introns_sup_protein
    exon_rna    = ($14 == "NA" ? 0 : $14 + 0)    #$14=pct_exons_sup_rnaseq
    exon_prot   = ($16 == "NA" ? 0 : $16 +0)    #$16=pct_exons_sup_protein

    hc = 0

    # Single-exon transcript:
    # HC if its exon is fully supported by RNA-seq OR protein.
    if (n_introns == 0) {
        if (exon_rna == 100 || exon_prot == 100)
            hc = 1
    }

    # Multi-exon transcript with substantial intron support:
    # retain as HC if either evidence type supports >=50% of introns.

    else {
	    if (intron_rna >= 50 || intron_prot >= 50)
		    hc = 1

    # No intron support from either source:
    # rescue as HC if all exons are supported by either RNA-seq or protein.
    else if (intron_rna == 0 && intron_prot == 0 &&
	    (exon_rna == 100 || exon_prot == 100))
    	hc = 1
    }
	if (hc)
		print
}
	

' "$input" > "$output"
