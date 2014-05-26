#!/bin/bash

exon=${1%.csv}"-exons.csv"
blast=${1%.csv}"-blast.csv"

python exonfinder_prot.py $1
echo "Done finding exons! Moving on to blasting stage..."
python exonblast.py $exon $2
echo "Done blasting exons! Now to align them..."
python blastalign.py $blast $2
