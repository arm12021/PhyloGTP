# region Imports

from ete3 import Tree
import os
import numpy as np
from itertools import repeat
import utilities
import time

# endregion Imports

# region Helper Functions


def get_taxa_freqs(gene_tree_str):
    gene_trees = [Tree(s, format=1) for s in gene_tree_str]

    taxa = {}
    for gt in gene_trees:
        # Remove the leaf names suffix from each leaf
        for leaf in gt:
            if "_" in leaf.name:
                leaf.name = leaf.name[: leaf.name.index("_")]

        added = []
        for leaf in gt:
            s = leaf.name
            if s not in taxa:
                taxa[s] = 0
            if s not in added:
                taxa[s] += 1
                added.append(s)

    taxa = list(taxa.items())
    taxa.sort(key=lambda x: x[1], reverse=True)

    return taxa


def get_pruned_genetrees(keep, gene_tree_str, weights):
    pruned_gts = []
    cur_weights = []
    for i, gts in enumerate(gene_tree_str):
        t = Tree(gts, format=1)
        # Find leaves whose base names (before underscore) match the keep list
        common = []
        for leaf in t:
            leaf_base_name = (
                leaf.name.split("_")[0] if "_" in leaf.name else leaf.name
            )
            if leaf_base_name in keep:
                common.append(leaf)
        # common = [l for l in t if l.name in keep]
        if len(common) >= 3:
            t.prune(common)
            pruned_gts.append(t.write(format=9))
            if weights:
                cur_weights.append(weights[i])
        else:
            pruned_gts.append(";")
    return pruned_gts, cur_weights


def write_gene_trees(pruned_gts, temp_dir, unrooted, use_ecceTERA=False):
    temp_gene_trees_path = os.path.join(temp_dir, f"GENETREES-temp")
    with open(temp_gene_trees_path, "w+") as f:
        for gt in pruned_gts:
            if unrooted and not use_ecceTERA:
                f.write(f"[&U]{gt}\n")
            else:
                f.write(f"{gt}\n")
    return temp_gene_trees_path


# endregion Helper Functions

# region Starting Tree Generation


# Build a starting tree with a random topology
def buildRandomStartingTree(names):
    t = Tree()
    t.populate(len(names), names_library=names)
    return t


# This function generates a starting tree for the given taxa.
def init_start_trees(taxa):
    staxa = taxa[:3]
    a, b, c = staxa[0][0], staxa[1][0], staxa[2][0]
    t1 = f"(({a},{b}),{c});"
    t2 = f"(({a},{c}),{b});"
    t3 = f"(({b},{c}),{a});"
    return [t1, t2, t3]


def additive_start_tree(
    gene_trees_path,
    gene_family_trees_path,
    dtl_args,
    weights,
    unrooted,
    temp_dir,
    num_cores,
    keep_temp=False,
):
    # region Init Start Tree

    getRecCost_time = 0
    localWrite_time = 0

    # If unrooted, get gene trees without [&U] prefix
    # [&U](((a,b),c),((d,a),((e,f),e))) -> (((a,b),c),((d,a),((e,f),e)))
    with open(gene_trees_path) as f:
        if unrooted and not dtl_args["use_ecceTERA"]:
            gene_tree_str = [
                line.strip()[line.index("]") + 1 :] for line in f.readlines()
            ]
        else:
            gene_tree_str = [line.strip() for line in f.readlines()]

    # Get a list of all taxa and their frequencies
    # taxa = [("a", 3), ("b", 3), ("c", 3), ("d", 3), ("e", 3)]
    taxa = get_taxa_freqs(gene_tree_str)

    # Initialize the starting trees from the first three (3) taxa in the taxa list
    # spec_trees = [ ((a,b),c), ((a,c),b), ((b,c),a) ]
    spec_trees = init_start_trees(taxa)

    # Keep the first three (3) taxa for the starting tree
    # keep = [ "a", "b", "c" ]
    keep = [x[0] for x in taxa[:3]]

    # Remove the first three (3) taxa from the taxa list
    # taxa = [("d", 3), ("e", 3), ("f", 3)]
    taxa = taxa[3:]

    # Prune the gene trees to only the leaves from the initial three (3) taxa (the keep list)
    # keep = [ "a", "b", "c" ]
    # (((a,b),c),((d,a),((e,f),e))) -> (((a,b),c),a)
    # pruned_gts = [ (((a,b),c),a), (((a,b),c),a), (((a,b),c),a) ]
    # Also trim the weights list (if necessary) to only the weights for the pruned gene trees
    pruned_gts, cur_weights = get_pruned_genetrees(keep, gene_tree_str, weights)

    # endregion Init Start Tree

    # region First Start Tree

    # If there are any pruned gene trees
    if len(pruned_gts) > 0:
        # Write the pruned gene trees to a temporary file (GENETREES-temp)
        # pruned_gts = [ (((a,b),c),a), (((a,b),c),a), (((a,b),c),a) ]
        temp_gene_trees_path = write_gene_trees(
            pruned_gts,
            temp_dir,
            unrooted,
            use_ecceTERA=dtl_args["use_ecceTERA"],
        )

        # For EACH of the starting/species trees (one at a time)
        # Use RANGER-DTL to calculate the reconciliation cost (score) to the pruned gene trees (all together)
        # spec_tree = "((a,b),c)" or "((a,c),b)" or "((b,c),a)"
        # pruned_gts = [ (((a,b),c),a), (((a,b),c),a), (((a,b),c),a) ]
        scores = [
            utilities.getTotalRecCost(
                dtl_args,
                spec_tree,
                temp_gene_trees_path,
                gene_family_trees_path,
                cur_weights,
                temp_dir,
                keep_temp,
            )
            for spec_tree in spec_trees
        ]

        if sum(scores) < 0:
            raise ValueError(
                "The scores are negative. Please check the input gene trees."
            )

        # Find the minimum score (cost) and the corresponding starting/species tree
        start_tree = spec_trees[np.argmin(scores)]
    else:
        while len(pruned_gts) == 0 and taxa:
            s = taxa.pop(0)
            keep.append(s[0])
        start_tree = Tree()
        start_tree.populate(len(names), names_library=names)
        start_tree = start_tree.write(format=9)

    # Create a MasterTree object from the starting tree and set the node names
    # mt = ((a,b)n1,c)
    mt = utilities.MasterTree(tree=Tree(start_tree))
    mt.set_node_names()

    # endregion First Start Tree

    print(f"Starting additive tree search with {len(taxa)} taxa")

    # region Start Tree Loop

    while taxa:
        s = taxa.pop(0)[0]
        keep.append(s)
        pruned_gts, cur_weights = get_pruned_genetrees(
            keep, gene_tree_str, weights
        )

        start_time = time.time()
        temp_gene_trees_path = write_gene_trees(
            pruned_gts,
            temp_dir,
            unrooted,
            use_ecceTERA=dtl_args["use_ecceTERA"],
        )
        localWrite_time += time.time() - start_time

        neighborhood = []
        order = [node for node in mt.tree.traverse()]
        p_node = Tree(name=s)
        w_node_name = f"n{len(mt.tree)}"
        for r_node in order:
            mt.regraft(p_node, r_node, w_node_name)
            neighborhood.append(mt.write())
            mt.prune(p_node)

        print(
            f"Adding {len(taxa)} taxa to tree with {len(neighborhood)} neighbors"
        )

        start_time = time.time()
        if num_cores > 1:
            p = utilities.multiprocessing.Pool(processes=num_cores)
            scores = p.starmap(
                utilities.getTotalRecCost,
                zip(
                    repeat(dtl_args),
                    neighborhood,
                    repeat(temp_gene_trees_path),
                    repeat(gene_family_trees_path),
                    repeat(cur_weights),
                    repeat(temp_dir),
                    repeat(keep_temp),
                ),
            )
        else:
            scores = [
                utilities.getTotalRecCost(
                    dtl_args,
                    spec_tree,
                    temp_gene_trees_path,
                    gene_family_trees_path,
                    cur_weights,
                    temp_dir,
                    keep_temp,
                )
                for spec_tree in neighborhood
            ]
        getRecCost_time += time.time() - start_time

        # print(f"localWrite Time: {localWrite_time}")
        print(f"getRecCost Time: {getRecCost_time}")
        if dtl_args["use_ecceTERA"]:
            print(f"ecceTERA Time: {utilities.eccetera_time}")
            print(f"Read Time Fam: {utilities.read_time_fam}")
            print(f"Read Time Genes: {utilities.read_time_gene}")
            print(f"Write Time Species: {utilities.write_time_spec}")
            print(f"Write Time Genes: {utilities.write_time_gene}")
            print(f"Total ecceTERA calls: {utilities.total_eccetera_calls}")
        else:
            print(f"RANGER-DTL Time: {utilities.ranger_time}")
            print(f"Write Time Species: {utilities.write_time_spec}")
            print(f"Total RANGER-DTL calls: {utilities.total_ranger_calls}")
        print(f"")

        min_score = min(scores)
        optimal_trees = [
            neighborhood[i]
            for i in range(len(scores))
            if scores[i] == min_score
        ]
        index = np.random.randint(len(optimal_trees))
        mt = utilities.MasterTree(tree=Tree(optimal_trees[index], format=8))

    # endregion Start Tree Loop

    return mt, min_score


# endregion Starting Tree Generation
