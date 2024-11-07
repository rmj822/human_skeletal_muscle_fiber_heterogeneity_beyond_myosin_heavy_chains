#!/bin/bash 
#$ -pe serial 6 
#$ -l h_vmem=10G
#$ -N AUC_pySCENIC

singularity exec /group/irc/shared/pyScenic/pySCENIC-0.12.1.sif \
    pyscenic aucell \
        /home/robinb/fiber_types/data/fiber_types.loom \
        /home/robinb/fiber_types/scenic-output/regulons_raw_v10.csv \
        -o /home/robinb/fiber_types/scenic-output/fiber_types_pySCENIC.loom \
        --num_workers 6
