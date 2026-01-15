#!/bin/bash

# This helper script dereplicates previously downloaded genomes using skDER.

echo Preparing to dereplicate genomes...

# If the derep_out already exists, delete it:
if [ -d derep_out ]; then
    rm -rf derep_out
fi

pi_cutoff=$1
pc_cutoff=$2
nb_cores=$3
genome_folder=$4
low_mem=$5

echo Dereplicating genomes in "$genome_folder" with identity cutoff of $pi_cutoff % and coverage cutoff of $pc_cutoff %
if [ $low_mem = "low_mem" ]; then
    echo Using low-memory mode
    mode="low_mem_greedy"
else
    mode="greedy"
fi

# Run skDER on the downloaded genomes and enable secondary clustering (-n flag):
echo -e "Starting skDER\n"

# Necessary for the regular expression to work in the following command
shopt -s nullglob

# Pass the genome files to skder
skder -g $genome_folder/*.{fna,fa,fasta,fna.gz,fa.gz,fasta.gz} -o derep_out -i $pi_cutoff -f $pc_cutoff -c $nb_cores -d $mode -n

# skDER stores the dereplicated genomes in its own output folder. Compare the amount of files in skder_out folder with initial folder where
# all genomes reside.
echo -e "\nDereplication done! $(ls "$genome_folder" | grep -E '.fasta|.fna|.fa|.fna.gz|.fasta.gz|.fa.gz' | wc -w) genomes were reduced to $(ls derep_out/Dereplicated_Representative_Genomes | wc -w) genomes"
