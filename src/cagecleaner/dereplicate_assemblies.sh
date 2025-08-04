#!/bin/bash

# This helper script dereplicates previously downloaded genomes using skDER.

echo Preparing to dereplicate genomes...

# Create the folder for skder output:
mkdir -p data

# If the download folder still exists, it will be deleted:
rm -rf data/downloads

pi_cutoff=$1
nb_cores=$2
genome_folder=$3

echo Dereplicating genomes in "$genome_folder"

# Run skDER on the downloaded genomes and enable secondary clustering (-n flag):
echo Dereplicating genomes with percent identity cutoff of $pi_cutoff.
echo -e "Starting skDER\n"

# Necessary for the regular expression to work in the following command
shopt -s nullglob

# Pass the genome files to skder
skder -g $genome_folder/*.{fna,fa,fasta,fna.gz,fa.gz,fasta.gz} -o data/skder_out -i $pi_cutoff -c $nb_cores -n


# skDER stores the dereplicated genomes in its own output folder. Compare the amount of files in skder_out folder with initial folder where
# all genomes reside.
echo -e "\nDereplication done! $(ls "$genome_folder" | grep -E '.fasta|.fna|.fa|.fna.gz|.fasta.gz|.fa.gz' | wc -w) genomes were reduced to $(ls data/skder_out/Dereplicated_Representative_Genomes | wc -w) genomes"
