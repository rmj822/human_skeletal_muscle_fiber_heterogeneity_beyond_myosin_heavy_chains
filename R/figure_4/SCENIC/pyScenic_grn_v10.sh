#!/bin/bash
#$ -pe serial 12
#$ -l h_vmem=12G
#$ -N GRN_pySCENIC

singularity exec /group/irc/shared/pyScenic/pySCENIC-0.12.1.sif \
    pyscenic grn \
    --num_workers 12 \
    -o /home/robinb/fiber_types/scenic-output/fiber_types_v10.adjacencies.tsv \
    /home/robinb/fiber_types/data/fiber_types.loom \
    /home/robinb/pyScenic_Robin/databases/allTFs_hg38.txt
