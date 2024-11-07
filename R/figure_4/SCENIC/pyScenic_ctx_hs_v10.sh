#!/bin/bash 
#$ -pe serial 12
#$ -l h_vmem=12G
#$ -N CTX_pySCENIC

singularity exec /group/irc/shared/pyScenic/pySCENIC-0.12.1.sif \
    pyscenic ctx \
        /home/robinb/fiber_types/scenic-output/fiber_types_v10.adjacencies.tsv \
        /home/robinb/pyScenic_Robin/databases/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
        /home/robinb/pyScenic_Robin/databases/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
        --annotations_fname /home/robinb/pyScenic_Robin/databases/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
        --expression_mtx_fname /home/robinb/fiber_types/data/fiber_types.loom \
        --mode "dask_multiprocessing" \
        --output /home/robinb/fiber_types/scenic-output/regulons_raw_v10.csv \
        --num_workers 12

