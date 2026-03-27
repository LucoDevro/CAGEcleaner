#!/bin/bash

# This helper script retrieves the NCBI Assembly IDs linked to a list of Nucleotide IDs using the NCBI E-utilities.
# It discerns between non-WGS entries and WGS entries. While the former contain sequence data and are directly linked to an assembly ID, 
# the latter are only indirectly linked to an assembly ID via its master record.
# Queries are launched in batches of 5000 in a subshell to avoid the 'too many arguments' error.

scaffolds_list=$1
mode=$2
batch_size=5000

## non-WGS entries
echo "--> Checking for non-WGS entries"
nb_non_wgs=$(cat $scaffolds_list | grep -vE '([A-Z]{4,}[0-9]{8,})\.[0-9]+' | wc -l)
echo "Found $nb_non_wgs non-WGS Nucleotide accession codes."

# Continue if not zero
if [[ $nb_non_wgs -gt 0 ]]
then
    echo "Querying NCBI..."
    # non-WGS entries can be directly fetched via a Nucleotide-Assembly link.
    case $mode in
        Genbank)
            command='elink -db nucleotide -target assembly -id "$@" | efetch -format docsum | xtract -pattern DocumentSummary -element Genbank >> non_wgs_assembly_accessions'
            ;;
        RefSeq)
            command='elink -db nucleotide -target assembly -id "$@" | efetch -format docsum | xtract -pattern DocumentSummary -element RefSeq >> non_wgs_assembly_accessions'
            ;;
        *)
            echo "Unknown fetching mode"
            ;;
    esac
    
    cat $scaffolds_list | grep -vE '([A-Z]{4,}[0-9]{8,})\.[0-9]+' | xargs -n $batch_size bash -c "$command"
    
    echo "Got $(cat non_wgs_assembly_accessions | wc -l) non-WGS Assembly accession codes"
else
    echo "Skipping NCBI querying"
fi

## WGS entries
echo "--> Checking for WGS entries"
nb_wgs=$(cat $scaffolds_list | grep -E '([A-Z]{4,}[0-9]{8,})\.[0-9]+' | wc -l)
echo "Found $nb_wgs WGS Nucleotide accession codes."

# Continue if not zero
if [[ $nb_wgs -gt 0 ]]
then
    # WGS entries need to be redirected to the master record holding all the separate scaffold records in Nucleotide
    # before linking to Assembly and fetching the assembly
    echo "Redirecting to WGS master record accession codes..."
    cat $scaffolds_list | grep -E '([A-Z]{4,}[0-9]{8,})\.[0-9]+' | 
    sed -E 's/[1-9][0-9]{5}\.[1-9]/000000/g' | 
    sed -E 's/[1-9][0-9]{4}\.[1-9]/00000/g' |
    sed -E 's/[1-9][0-9]{3}\.[1-9]/0000/g' |
    sed -E 's/[1-9][0-9]{2}\.[1-9]/000/g' |
    sed -E 's/[1-9][0-9]\.[1-9]/00/g' |
    sed -E 's/[1-9]\.[1-9]/0/g' > wgs_masters
    
    echo "Querying NCBI..."
    case $mode in
        Genbank)
            command='elink -db nucleotide -target assembly -id "$@" | efetch -format docsum | xtract -pattern DocumentSummary -element Genbank >> wgs_assembly_accessions'
            ;;
        RefSeq)
            command='elink -db nucleotide -target assembly -id "$@" | efetch -format docsum | xtract -pattern DocumentSummary -element RefSeq >> wgs_assembly_accessions'
            ;;
        *)
            echo "Unknown fetching mode"
            ;;
    esac
    
    cat wgs_masters | xargs -n $batch_size bash -c "$command"
    
    echo "Got $(cat wgs_assembly_accessions | wc -l) WGS Assembly accession codes"
else
    echo "Skipping NCBI querying"
fi

## concatenating the retrieved accession codes
echo "Merging results..."
cat *wgs_assembly_accessions | sort -u > assembly_accessions.txt
echo "Found $(cat assembly_accessions.txt | wc -l) accession codes after removing duplicates"

## Cleaning up
rm *wgs_*
