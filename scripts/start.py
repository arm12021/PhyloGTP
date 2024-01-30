from ete3 import Tree
import os
import numpy as np
from itertools import repeat
from utilities import *

def get_taxa_freqs(gene_tree_str):
    gene_trees = [Tree(s, format=1) for s in gene_tree_str]

    taxa = {}
    for gt in gene_trees:
        added = []
        for leaf in gt:
            s = leaf.name
            if s not in taxa:
                taxa[s] = 0
            if s not in added:
                taxa[s] += 1
                added.append(s)

    taxa = list(taxa.items())
    taxa.sort(key = lambda x: x[1], reverse=True)

    return taxa

def init_start_trees(taxa):
    staxa = taxa[:3]
    a,b,c = staxa[0][0], staxa[1][0], staxa[2][0]
    t1 = f'(({a},{b}),{c});'
    t2 = f'(({a},{c}),{b});'
    t3 = f'(({b},{c}),{a});'
    return [t1, t2, t3]

def get_pruned_genetrees(keep, gene_tree_str, weights):
    pruned_gts = []
    cur_weights = []
    for i,gts in enumerate(gene_tree_str):
        t = Tree(gts, format=1)
        common = [l for l in t if l.name in keep]
        if len(common) >= 3:
            t.prune(common)
            pruned_gts.append(t.write(format=9))
            if weights:
                cur_weights.append(weights[i])
    return pruned_gts, cur_weights

def write_gene_trees(pruned_gts, temp_dir, unrooted):
    temp_gene_trees_path = os.path.join(temp_dir, f'GENETREES-temp')
    with open(temp_gene_trees_path, 'w+') as f:
        for gt in pruned_gts:
            if unrooted:
                f.write(f'[&U]{gt}\n')
            else:
                f.write(f'{gt}\n')
    return temp_gene_trees_path


def additive_start_tree(gene_trees_path, ranger_args, weights, unrooted, temp_dir, num_cores):
    with open(gene_trees_path) as f:
        if unrooted:
            gene_tree_str = [line.strip()[line.index(']')+1:] for line in f.readlines()]
        else:
            gene_tree_str = [line.strip() for line in f.readlines()]
    
    taxa = get_taxa_freqs(gene_tree_str)
    spec_trees = init_start_trees(taxa)
    keep = [x[0] for x in taxa[:3]]
    taxa = taxa[3:]
    pruned_gts, cur_weights = get_pruned_genetrees(keep, gene_tree_str, weights)

    if len(pruned_gts) > 0:
        temp_gene_trees_path = write_gene_trees(pruned_gts, temp_dir, unrooted)
        scores = [getTotalRecCost(ranger_args, spec_tree, temp_gene_trees_path, cur_weights, temp_dir) for spec_tree in spec_trees]
        start_tree = spec_trees[np.argmin(scores)]
    else:
        while len(pruned_gts) == 0 and taxa:
            s = taxa.pop(0)
            keep.append(s[0])
        start_tree = Tree()
        start_tree.populate(len(names), names_library=names)
        start_tree = start_tree.write(format=9)
    
    mt = MasterTree(tree=Tree(start_tree))
    mt.set_node_names()

    while taxa:
        s = taxa.pop(0)[0]
        keep.append(s)
        pruned_gts, cur_weights = get_pruned_genetrees(keep, gene_tree_str, weights)
        temp_gene_trees_path = write_gene_trees(pruned_gts, temp_dir, unrooted)

        neighborhood = []
        order = [node for node in mt.tree.traverse()]
        p_node = Tree(name=s)
        w_node_name = f'n{len(mt.tree)}'
        for r_node in order:
            mt.regraft(p_node, r_node, w_node_name)
            neighborhood.append(mt.write())
            mt.prune(p_node)
        
        if num_cores > 1:
            p = multiprocessing.Pool(processes=num_cores)
            scores = p.starmap(getTotalRecCost, zip(repeat(ranger_args), neighborhood, repeat(temp_gene_trees_path), repeat(cur_weights), repeat(temp_dir)))
        else:
            scores = [getTotalRecCost(ranger_args, spec_tree, temp_gene_trees_path, cur_weights, temp_dir) for spec_tree in neighborhood]
    
        min_score = min(scores)
        optimal_trees = [neighborhood[i] for i in range(len(scores)) if scores[i] == min_score]
        index = random.randint(0, len(optimal_trees)-1)
        mt = MasterTree(tree=Tree(optimal_trees[index], format=8))
    
    return mt, min_score
