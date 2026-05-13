#!/bin/bash

function analyze {
    python generate_mini_dataset.py \
        -p ${PROTEIN} \
        -c ${CHROM} \
        -g ${GENOME} \
        --prefix ${PREFIX}
}

PROTEIN=/lisc/scratch/molevo/dts/ODP_genomes/GenDB_annotated_chr/odp_ncbi_genome_db/output/source_data/annotated_genomes/GCF_022113875.1/GCF_022113875.1.chrFilt.pep.gz
CHROM=/lisc/scratch/molevo/dts/ODP_genomes/GenDB_annotated_chr/odp_ncbi_genome_db/output/source_data/annotated_genomes/GCF_022113875.1/GCF_022113875.1.chrFilt.chrom.gz
GENOME=/lisc/scratch/molevo/dts/ODP_genomes/GenDB_annotated_chr/odp_ncbi_genome_db/output/source_data/annotated_genomes/GCF_022113875.1/GCF_022113875.1.chr.fasta.gz
PREFIX=mini_hydra/mini_hydra
mkdir -p mini_hydra
analyze

PROTEIN=/lisc/scratch/molevo/dts/ODP_genomes/GenDB_annotated_chr/odp_ncbi_genome_db/output/source_data/annotated_genomes/GCF_018143015.1/GCF_018143015.1.chrFilt.pep.gz
CHROM=/lisc/scratch/molevo/dts/ODP_genomes/GenDB_annotated_chr/odp_ncbi_genome_db/output/source_data/annotated_genomes/GCF_018143015.1/GCF_018143015.1.chrFilt.chrom.gz
GENOME=/lisc/scratch/molevo/dts/ODP_genomes/GenDB_annotated_chr/odp_ncbi_genome_db/output/source_data/annotated_genomes/GCF_018143015.1/GCF_018143015.1.chr.fasta.gz
PREFIX=mini_urchin/mini_urchin
mkdir -p mini_urchin
analyze


