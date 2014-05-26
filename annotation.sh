#!/bin/bash

exon=${1%.csv}"-exons.csv"
blast=${1%.csv}"-blast.csv"

python exonfinder_prot.py $1
echo "All exons found! Moving on to blasting stage..."
python exonblast.py $exon $2
echo "All exons blasted! Now to align them..."
python blastalign.py $blast $2
