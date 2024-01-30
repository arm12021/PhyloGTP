#!/bin/bash

n_genes=100
seq_length=400

tree_path = f'sage_data/{inp}/gene_trees/gene.{i}.{rate}.pruned.tree'
out_path = iq_data/{inp}/genes/{i}

out_path=out_inferred
mkdir $out_path
for ((i=1;i<=$n_genes;i++)); do
    tree_path=out_true/gene_trees/gene.${i}.pruned.tree
    iqtree2 --alisim ${out_path}/${i}.${seq_length} -m GTR -t $tree_path --length $seq_length
    iqtree2 -s ${out_path}/${i}.${seq_length}.phy -m JC --prefix ${out_path}/${i}.${seq_length}.phy
done
