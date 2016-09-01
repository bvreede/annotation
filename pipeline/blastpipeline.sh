#!/bin/bash

#This script is a pipeline to blast a set of sequences, parse the output
#files, and collect the appropriate sequences.
#The script works with NCBI blast 2.2.29+ in python 2.7.
#Author: Barbara Vreede
#Contact: b.vreede@gmail.com
#Date: 17 August 2014

if [ -z $5 ];

then
echo ""
echo "Script usage: bash blastpipeline.sh species infolder blasttype extra_nt max_scaf"
echo "--------------------------------------------------------------------------------"
echo "species = 4 letter abbreviation of species (make sure the scripts 'blast.py' and 'getseq.py' contain custom path info!)"
echo "infolder = the folder containing the sequences for blasting"
echo "blasttype = the type of blast required (eg. tblastn)"
echo "extra_nt = the number of additional nucleotides to add to retrieved sequences"
echo "max_scaf = the maximum number of sequences (scaffold) per blast hit to retrieve"
echo "";

else
species=$1
infolder=$2
blasttype=$3
extrant=$4
maxscaf=$5

outfolder=$2"_blast2"$1

time python blast.py $species $infolder $blasttype
echo "Blast algorithm complete! Analysing data..."
time python parse_data.py $outfolder
echo "Data parsed. Retrieving sequences..."
time python parse_seq.py $species $outfolder $extrant $maxscaf;

fi

