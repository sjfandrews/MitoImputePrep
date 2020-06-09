'''Snakefile for Hyperparameter gridsearch in UMAP'''

shell.prefix('module load R/3.6.2; ')

COHORT = 'mtref'
PCS = [2,3,4,5,6,7,8,9,10]
N_COMPONENTS = [2]
MIN_DIST = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]
N_NEIGHBOR = [15, 25, 50, 75, 100]

# COHORT = 'mtref'
# PCS = [2,10]
# N_COMPONENTS = [2]
# MIN_DIST = [0.5]
# N_NEIGHBOR = [15]

rule all:
    input:
        expand('data/test/hyperparameters/umpap_{cohort}_{pcs}pcs_{n_components}components_{min_dist}distance_{n_neighbor}neigbor.txt', cohort=COHORT, pcs=PCS, n_components=N_COMPONENTS, min_dist=MIN_DIST, n_neighbor=N_NEIGHBOR)

rule umap:
    input:
        infile = 'data/test/mtref_mtped.txt'
    output:
        outfile = 'data/test/hyperparameters/umpap_{cohort}_{pcs}pcs_{n_components}components_{min_dist}distance_{n_neighbor}neigbor.txt'
    params:
        pcs = '{pcs}',
        n_components = '{n_components}',
        min_dist = '{min_dist}',
        n_neighbor = '{n_neighbor}'
    script:
        'scripts/R/umap_tune.R'
