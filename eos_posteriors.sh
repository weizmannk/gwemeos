#!/bin/bash
source /cvmfs/oasis.opensciencegrid.org/ligo/sw/conda/etc/profile.d/conda.sh
conda activate nmma
python eos_posteriors.py $@
conda deactivate 
