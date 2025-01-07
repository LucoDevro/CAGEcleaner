#!/bin/bash

# This helper script retrieves the NCBI Assembly IDs linked to a list of Nucleotide IDs using the NCBI E-utilities.
# It discerns between non-WGS entries and WGS entries. While the former contain sequence data and are directly linked to an assembly ID, 
# the latter are only indirectly linked to an assembly ID via its master record.
# Queries are launched in batches of 5000 in a subshell to avoid the 'too many arguments' error.

scaffolds_list=$1
batch_size=5000

## non-WGS entries
echo "--> Checking for non-WGS entries"
echo "Found $(cat $scaffolds_list | grep -E '_|^[A-Z|0-9]{8}\.[1-9]' | wc -l) non-WGS Nucleotide accession codes."
echo "Querying NCBI..."
# non-WGS entries can be directly fetched via a Nucleotide-Assembly link.
cat $scaffolds_list | grep -E '_|^[A-Z|0-9]{8}\.[1-9]' | xargs -n $batch_size bash -c 'elink -db nucleotide -target assembly -id "$@" | efetch -format docsum | xtract -pattern DocumentSummary -element AssemblyAccession >> non_wgs_assembly_accessions'
echo "Got $(cat non_wgs_assembly_accessions | wc -l) non-WGS Assembly accession codes"

## WGS entries
echo "--> Checking for WGS entries"
echo "Found $(cat $scaffolds_list | grep -vE '_|^[A-Z|0-9]{8}\.[1-9]' | wc -l) WGS Nucleotide accession codes."
# WGS entries need to be redirected to the master record holding all the separate scaffold records in Nucleotide
# before linking to Assembly and fetching the assembly
echo "Redirecting to WGS master record accession codes..."
cat $scaffolds_list | grep -vE '_|^[A-Z|0-9]{8}\.[1-9]' | 
sed -E 's/[1-9][0-9]{5}\.[1-9]/000000/g' | 
sed -E 's/[1-9][0-9]{4}\.[1-9]/00000/g' |
sed -E 's/[1-9][0-9]{3}\.[1-9]/0000/g' |
sed -E 's/[1-9][0-9]{2}\.[1-9]/000/g' |
sed -E 's/[1-9][0-9]\.[1-9]/00/g' |
sed -E 's/[1-9]\.[1-9]/0/g' > wgs_masters
echo "Querying NCBI..."
cat wgs_masters | xargs -n $batch_size bash -c 'elink -db nucleotide -target assembly -id "$@" | efetch -format docsum | xtract -pattern DocumentSummary -element AssemblyAccession >> wgs_assembly_accessions'
echo "Got $(cat wgs_assembly_accessions | wc -l) WGS Assembly accession codes"

## concatenating the retrieved accession codes
echo "Merging results..."
cat non_wgs_assembly_accessions wgs_assembly_accessions | sort -u > assembly_accessions
echo "Found $(cat assembly_accessions | wc -l) accession codes after removing duplicates"

## Cleaning up
rm *wgs_*
