#!/bin/bash

exon=${1%.csv}"-exons.csv"
blast=${1%.csv}"-blast.csv"

python exonfinder_prot.py $1
python exonblast.py $exon $2
python blastalign.py $blast $2
