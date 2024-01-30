#!/bin/bash

sage_path=sagephy/sagephy-1.0.0.jar
out_path=out_true
n_leaf=50
n_genes=100
birth=5.0
death=2.5

d=0.3
l=0.6
t=0.6

mkdir $out_path
java -jar $sage_path HostTreeGen -nox -min $n_leaf -max $n_leaf -a 100000 1.0 $birth $death ${out_path}/spec

mkdir ${out_path}/gene_trees
for ((i=1;i<=$n_genes;i++)); do
    java -jar $sage_path GuestTreeGen -min 10 -nox  ${out_path}/spec.pruned.tree $d $l $t ${out_path}/gene_trees/gene.${i}
done