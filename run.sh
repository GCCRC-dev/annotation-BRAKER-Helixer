# pipeline to generate merged HC and LC from both BRAKER4 and Helixer
# each script calls different software so check the dependencies and include them in a conda environment or define container files

set -euo pipefail

#### 1. Change paths ####
workdir="$HOME/annotation/vtub/hap1"
path_braker_res="$HOME/annotation/vtub/braker/output/vtub/results"
path_helixer_res="$HOME/annotation/vtub/helixer/output"

# path to bam file with all merged rna-seq alignments agaisnt assembly
bam="$HOME/rnaseq/vtub/merged_rnaseq_vtub_hap1.bam"

# path to genome fasta
genome="$HOME/assembly/vtub/Vtub.hap1.v1.softmasked.fasta"

helixer_gff="$HOME/annotation/vtub/helixer/output/Vtub.hap1.gff3.gz"
braker_gff="$HOME/annotation/vtub/braker/output/vtub/results/braker.gff3.gz"
braker_longest="$HOME/annotation/vtub/braker/output/vtub/results/braker.longest.aa.gz"
braker_gene_sup="$HOME/annotation/vtub/braker/output/vtub/results/gene_support.tsv"

# path to where the scripts are
script_dir="$HOME/script_dir"

cd ${workdir}
mkdir -p curating && cd curating
mkdir -p braker helixer final


#### 2. BRAKER curation ####
cd ${workdir}/curating/braker

ln -sfn $braker_longest .
ln -sfn $braker_gene_sup .
gzip -dc "$braker_gff" > braker.gff3

# PROTEIN SUPPORT
${script_dir}/run_diamond_3dbs.sh braker braker.longest.aa.gz

# RNA SUPPORT
${script_dir}/parse_gene_support_braker_HC.sh

# merge everything
python ${script_dir}/classify_braker_confidence.py \
    --gff braker.gff3 \
    --structural-hc out_gene_support_HC_candidates.out \
    --swissprot braker.swissprot.strong.tsv \
    --trembl braker.trembl.strong.tsv \
    --ptrep braker.ptrep.strong.tsv \
    --output braker.confidence.tsv \
    --hc-ids braker.HC.gene_ids.txt \
    --lc-ids braker.LC.gene_ids.txt


#### 3. HELIXER curation ####
cd ${workdir}/curating/helixer
ln -sfn ../braker/braker.gff3 braker.gff3
gzip -dc "$helixer_gff" > "$(basename "${helixer_gff%.gz}")"

# obtain Helixer.u
${script_dir}/get_helixer_u.sh $genome "$(basename "${helixer_gff%.gz}")"

# PROTEIN SUPPORT
${script_dir}/run_diamond_3dbs.sh helixer helixer.u.faa

# RNA SUPPORT
mkdir -p ${workdir}/curating/helixer/rna && cd ${workdir}/curating/helixer/rna

gff_u='../helixer.u.gff3' 
${script_dir}/obtain_rnaseq_support_helixer.sh $bam $gff_u

Rscript ${script_dir}/merge_rnaseq_helixer.Rscript 

cd ${workdir}/curating/helixer

# merge everything
python ${script_dir}/classify_helixer_confidence.py \
    --proteins helixer.u.faa \
    --rnaseq rna/helixer.u_rnaseq_support.tsv \
    --swissprot helixer.swissprot.strong.tsv \
    --trembl helixer.trembl.strong.tsv \
    --ptrep helixer.ptrep.strong.tsv \
    --output helixer.u.confidence.tsv \
    --hc-ids helixer.u.HC.ids.txt \
    --lc-ids helixer.u.LC.ids.txt

#### 4. Merge HC and LC, validate, rename genes, generate final files ####
cd ${workdir}/curating/final

braker_gff="${workdir}/curating/braker/braker.gff3"
helixer_gff="${workdir}/curating/helixer/helixer.u.gff3"
${script_dir}/get_final_gffs.sh $genome $braker_gff $helixer_gff

