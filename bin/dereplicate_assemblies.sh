#!/bin/bash

# This helper script dereplicates previously downloaded genomes using skDER.

echo Preparing to dereplicate genomes...

# If the download folder still exists, it will be deleted:
rm -rf data/downloads

pi_cutoff=$1
nb_cores=$2

# Run skDER on the downloaded genomes and enable secondary clustering (-n flag):
echo Dereplicating genomes with percent identity cutoff of $pi_cutoff.
echo -e "Starting skDER\n"
skder -g data/genomes/* -o data/skder_out -i $pi_cutoff -c $nb_cores -n

# skDER stores the dereplicated genomes in its own output folder. Compare the amount of files in skder_out folder with initial folder where
# all genomes reside.
echo -e "\nDereplication done! $(ls data/genomes | wc -w) genomes were reduced to $(ls data/skder_out/Dereplicated_Representative_Genomes | wc -w) genomes"
