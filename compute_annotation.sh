#!/bin/bash
#$ -S /bin/bash
#$ -M amine.ghozlane@pasteur.fr
#$ -m bea
#$ -q fast
#$ -p common,dedicated
#$ -pe thread 1
#$ -l mem_total=2G
#$ -J metahit_

source /local/gensoft2/adm/etc/profile.d/modules.sh
module purge
module add blast+/2.6.0 Python/3.6.0 taxadb/0.6.0
## -evalue 1E-3
blastn -query toannot/metahit_v2_${SLURM_ARRAY_TASK_ID}.fna \
    -db /local/databases/fasta/nt -num_threads 1 \
    -out blast3/metahit_v2_${SLURM_ARRAY_TASK_ID}_nt.tsv \
    -max_target_seqs 1 -max_hsps 1 \
    -outfmt '6 qseqid sseqid qlen length mismatch gapopen qstart qend sstart send pident qcovs evalue bitscore'\
    -perc_identity 50 -qcov_hsp_perc 50
#-word_size 100
$HOME/split_fasta/get_taxonomy3.py -i blast3/metahit_v2_${SLURM_ARRAY_TASK_ID}_nt.tsv \
    -d /pasteur/services/policy01/banques/prod/rel/taxadb/taxadb_2018-05-01/db/taxadb_full.sqlite \
    -o annotation3/metahit_v2_${SLURM_ARRAY_TASK_ID}_taxonomy.tsv
$HOME/split_fasta/ExtractNCBIDB2.py -f blast3/metahit_v2_${SLURM_ARRAY_TASK_ID}_nt.tsv \
    -g annotation3/metahit_v2_${SLURM_ARRAY_TASK_ID}_taxonomy.tsv \
    -nb 1 -o annotation3/metahit_v2_${SLURM_ARRAY_TASK_ID}_annotation.tsv

