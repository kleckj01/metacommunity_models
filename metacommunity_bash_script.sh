#!/bin/bash

module load octave
bsub -J "meta[1-1000]" -o %J.out -e %J.err -n 1 -R "rusage[mem=1024]" -W 1:00 octave --eval "metacommunity('metacommunity_\$LSB_JOBINDEX.mat')"

