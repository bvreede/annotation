#!/bin/bash

#This script is a pipeline to blast a set of sequences, parse the output
#files, and collect the appropriate sequences.
#The script works with NCBI blast 2.2.29+ in python 2.7.
#Author: Barbara Vreede
#Contact: b.vreede@gmail.com
#Date: 17 August 2014

echo "Script usage: sh blastpipeline.sh species infolder blasttype extrant"

species=$1
infolder=$2
blasttype=$3
extrant=$4

outfolder=$2"_blast2"$1

time python blast.py $species $infolder $blasttype
echo "Blast algorithm complete! Analysing data..."
time python parse_data.py $outfolder
echo "Data parsed. Retrieving sequences..."
time python parse_seq.py $species $outfolder $extrant


