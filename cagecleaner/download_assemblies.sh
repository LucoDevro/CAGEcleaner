#!/bin/bash

# This helper script reads a formatted text file with batches of assembly IDs and downloads the corresponding genomes.

echo Preparing to download assemblies...

# If intermediary/output files and folders already exist, they will be deleted:
rm -rf data

# Set some useful variables:
assembly_count=$(cat download_batches.txt | wc -w)  # Amount of words corresponds to amount of assembly IDs.
batch_count=$(cat download_batches.txt | wc -l)  # Each line is one batch. Enforced by the parent python script.

echo Got $assembly_count assembly IDs. Downloading genomes in $batch_count batches...

# Create the directory to store all genomes
mkdir genomes downloads

# Here we loop over each line in the 'download_batches.txt' file:
batch_counter=1
while read line; do
    echo -e "\nBatch $batch_counter/$batch_count"

    # Download the dehydrated genome package. xargs removes the trailing whitespace:
    datasets download genome accession $(echo "$line" | xargs) --dehydrated

    # Unzip the files
    unzip -d downloads ncbi_dataset.zip && rm ncbi_dataset.zip

    # Rehydrate
    datasets rehydrate --directory downloads --gzip

    # Put all genomes from this batch into the collection directory
    mv downloads/ncbi_dataset/data/GC*/* genomes/
    
    # Clean up
    rm -rf downloads
    
    batch_counter=$((batch_counter+1))
done < download_batches.txt

echo -e "\nDownloading finished!"

# Print amount of files in the genomes folder. This should equal the initial amount of assembly IDs
echo $(ls genomes | wc -w) genomes in $(pwd)/genomes.
echo Directory size: $(du --human-readable -s genomes/).  # Size of the genomes folder

# Cleaning up
rm -rf data/downloads
