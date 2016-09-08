#! /bin/bash

# Etienne et al. 2014:
#
#     we simulated 1000 phylo- genetic trees under the protracted speciation
#     model for various sets of parameters (b=0.5, λ=0.1,0.3,1, μ1 =μ2 =μ= 0,
#     0.1, 0.2), and a fixed crown age of 5, 10, or 15 My, condi- tional on the
#     realized tree retaining the initial root (i.e., survival of both original
#     crown lineages).
#
#   b1 = b2 = {0.5}
#   c1 = 0.1, 1.0, 10, 100, 100000
#   e1 = e2 = {0, 0.1, 0.2}
#
# Etienne et al. 2014, supp. mat.:
#
#   S3:
#   speciation completion rate = {1e-03, 1e-01, 1e+01, 1e+03}

SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
TREE_SCRIPT="${SCRIPT_DIR}/e01_generate_containing_trees.py"
BPP_SCRIPT="${SCRIPT_DIR}/e02_generate_bpp_analyses.py"
CUR_DIR=${PWD}

mkdir -p 01-trees && cd 01-trees || exit
for birth_rate in 0.5 ; do
    for extinction_rate in 0.0 0.2 ; do
        for conversion_rate in 0000.001 0000.100 0001.000 0010.000 1000.000; do
            for sim_time in 10 ; do
                title="b${birth_rate}_e${extinction_rate}_c${conversion_rate}_t${sim_time}"
                echo $title
                python ${TREE_SCRIPT} --title ${title} --b1 ${birth_rate} --b2 ${birth_rate} --c1 ${conversion_rate} --e1 ${extinction_rate} --e2 ${extinction_rate} --max-time ${sim_time} --nreps 20 || exit
            done
        done
    done
done
cd ..
mkdir -p 02-bpp && cd 02-bpp || exit
# cd 02-bpp || exit
# for mut_rate in 1e-06 1e-08; do
for mut_rate in 1e-08; do
    subtitle="m${mut_rate}"
    python ${BPP_SCRIPT} --title ${subtitle} --population-size 10000 --num-individuals-per-population 4 --num-loci-per-individual 10 --num-characters-per-locus 1000 --mutation-rate-per-site ${mut_rate} ../01-trees/*.lineages.tre
done

