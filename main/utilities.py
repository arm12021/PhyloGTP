# region Imports

from ete3 import Tree
import subprocess
import shutil
import time
import multiprocessing
import os
import numpy as np

# endregion Imports

# region MasterTree Class


class MasterTree:
    def __init__(self, tree=Tree()):
        self.tree = tree

    # Print the tree structure in ASCII format
    def print_structure(self):
        print(self.tree.get_ascii(show_internal=True))

    # Print the tree in newick format
    # The default format is 9 (newick with internal node names)
    # The format 8 is newick with internal node names and the root name
    def print_newick(self, form=9):
        if form == 8:
            root_name = self.tree.name
            tree_str = self.tree.write(format=8)
            tree_str = tree_str[:-1] + root_name + ";"
            print(tree_str)
        else:
            print(self.tree.write(format=form))

    # Write the tree to a string in the specified format
    # The default format is 8 (newick with internal node names)
    def write(self, form=8):
        if form == 8:
            root_name = self.tree.name
            tree_str = self.tree.write(format=8)
            tree_str = tree_str[:-1] + root_name + ";"
            return tree_str
        else:
            return self.tree.write(format=form)

    def save(self, path, format=9):
        self.tree.write(outfile=path, format=format)

    def n(self, name):
        return self.tree.search_nodes(name=name)[0]

    def copy(self):
        return MasterTree(self.tree.copy("newick"))

    def get_sibling(self, node):
        if node.is_root():
            return None
        else:
            for c in node.up.children:
                if c != node:
                    return c

    def add_node(self, locationNode, newNodeName, internalName):
        node = self.n(locationNode)

        if node.is_root():
            new_internal = Tree(name=internalName)
            new_internal.add_child(self.tree)
            new_internal.add_child(name=newNodeName)
            self.tree = new_internal
        else:
            parent = node.up
            node.detach()
            new_internal = parent.add_child(name=internalName)
            new_internal.add_child(node)
            new_internal.add_child(name=newNodeName)

    def prune(self, p_node):
        # def prune(self, p_node_node):
        # p_node = self.tree&p_node_node
        assert not p_node.is_root(), "pruned node cannot be the root"
        w_node = p_node.up
        w_node_name = w_node.name

        p_node.detach()
        sibling = w_node.children[0]

        if w_node.is_root():
            # sibling = w_node.get_children()[0].detach()
            sibling.detach()
            self.tree = sibling
        # else:
        # w_node.delete()
        w_node.delete()

        return w_node_name, sibling
        # return p_node, w_node_name

    def regraft(self, p_node, r_node, w_node_name):
        if r_node.up:
            newNode = r_node.up.add_child(name=w_node_name)
            r_node.detach()
            newNode.add_child(r_node)
            newNode.add_child(p_node)
        else:
            newNode = Tree(name=w_node_name)
            newNode.add_child(self.tree)
            newNode.add_child(p_node)
            self.tree = newNode

    def set_node_names(self):
        c = 0
        for node in self.tree.traverse():
            if not node.is_leaf():
                node.name = "n" + str(c)
                c += 1

    def unique_leafs(self):
        # Returns a list of unique leaf names in the tree
        # This is useful for checking if the tree has duplicate leaves
        return list(set([leaf.name for leaf in self.tree.get_leaves()]))

    def duplicate_leaves(self):
        # Returns a list of leaf names that appear more than once in the tree
        leaf_names = [leaf.name for leaf in self.tree.get_leaves()]
        name_counts = {}

        # Count occurrences of each name
        for name in leaf_names:
            if name in name_counts:
                name_counts[name] += 1
            else:
                name_counts[name] = 1

        # Return only names that appear multiple times
        return [name for name, count in name_counts.items() if count > 1]

    def check_duplicate_leaves(self):
        # Check if the tree has duplicate leaves
        # Returns True if there are duplicate leaves, False otherwise
        return len(self.tree.get_leaves()) != len(self.unique_leafs())


# endregion MasterTree Class

# region ecceTERA Calls


def callEcceTERA(ecce_args, temp_path_species, temp_path_genes):
    global total_eccetera_calls
    # Make system call to ecceTERA
    d_cost = ecce_args["D"]
    t_cost = ecce_args["T"]
    l_cost = ecce_args["L"]
    params = [
        ecce_args["path"],
        f"species.file={temp_path_species}",
        f"gene.file={temp_path_genes}",
        f"dupli.cost={d_cost}",
        f"HGT.cost={t_cost}",
        f"loss.cost={l_cost}",
        f"dated=0",
        # f"resolve.trees={ecce_args['unrooted']}",
        f"amalgamate={ecce_args['amalgamate']}",
    ]
    # print(f"Calling ecceTERA with params: {params[1:]}")

    try:
        raw_out = subprocess.run(
            params,
            check=True,
            capture_output=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as e:
        print(f"Error calling ecceTERA!")
        print(f"{e}")
        print(f"{e.stderr}")
        return -1
    except Exception as e:
        print(f"Error calling ecceTERA!")
        print(f"{e}")
        return -1

    # print(raw_out)
    # Expected output: "Cost of a most parsimonious reconciliation: 1"
    ran_out = raw_out.stdout.strip()
    total_eccetera_calls += 1

    ran_out = raw_out.stdout.strip().split("\n")
    costs = []
    for line in ran_out:
        if "Cost of a most parsimonious reconciliation" in line:
            c = int(line[line.index(":") + 2 :])
            costs.append(c)

    # print(f"Reconciliation cost: {sum(costs)}")
    return sum(costs)


def callEcceTERAWeighted(
    ecce_args, temp_path_species, temp_path_genes, weights
):
    global total_eccetera_calls
    # Make system call to ecceTERA
    d_cost = ecce_args["D"]
    t_cost = ecce_args["T"]
    l_cost = ecce_args["L"]
    params = [
        ecce_args["path"],
        f"species.file={temp_path_species}",
        f"gene.file={temp_path_genes}",
        f"dupli.cost={d_cost}",
        f"HGT.cost={t_cost}",
        f"loss.cost={l_cost}",
        f"dated=0",
        # f"resolve.trees={ecce_args['unrooted']}",
        f"amalgamate={ecce_args['amalgamate']}",
    ]
    # print(f"Calling ecceTERA with params: {params[1:]}")

    try:
        raw_out = subprocess.run(
            params,
            check=True,
            capture_output=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as e:
        print(f"Error calling ecceTERA!")
        print(f"{e}")
        print(f"{e.stderr}")
        return -1
    except Exception as e:
        print(f"Error calling ecceTERA!")
        print(f"{e}")
        return -1

    # print(raw_out)
    ran_out = raw_out.stdout.strip().split("\n")
    total_eccetera_calls += 1

    costs = []
    for line in ran_out:
        if "Cost of a most parsimonious reconciliation" in line:
            c = int(line[line.index(":") + 2 :])
            costs.append(c)
    weighted_costs = [costs[i] * weights[i] for i in range(len(costs))]
    return sum(weighted_costs)


# endregion ecceTERA Calls

# region Ranger Calls


# Makes system call to ranger and returns DTL score
def callRanger(ranger_args, temp_path):
    global total_ranger_calls
    # Make system call to ranger
    params = [
        ranger_args["path"],
        "-i",
        temp_path,
        "-D",
        ranger_args["D"],
        "-T",
        ranger_args["T"],
        "-L",
        ranger_args["L"],
        "-q",
        "-s",
    ]
    # print(f"Calling Ranger with params: {params[1:]}")
    raw_out = subprocess.run(
        params,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        universal_newlines=True,
    )
    # print(raw_out)
    # raw_out = subprocess.run(params, check=True, stdout=subprocess.PIPE, universal_newlines=True)
    ran_out = raw_out.stdout.strip()
    total_ranger_calls += 1
    # print(ran_out)
    # Expected output: "Total reconciliation cost: 12 (Duplications: 0, Transfers: 3, Losses: 0)"
    rec_cost = int(ran_out[ran_out.index(":") + 2 : ran_out.index("(") - 1])
    return rec_cost


def callRangerWeighted(ranger_args, temp_path, weights):
    global total_ranger_calls
    # Make system call to ranger
    params = [
        ranger_args["path"],
        "-i",
        temp_path,
        "-D",
        ranger_args["D"],
        "-T",
        ranger_args["T"],
        "-L",
        ranger_args["L"],
        "-q",
    ]
    raw_out = subprocess.run(
        params,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        universal_newlines=True,
    )
    # print(raw_out)
    ran_out = raw_out.stdout.strip().split("\n")
    total_ranger_calls += 1
    costs = []
    for line in ran_out:
        if "The minimum reconciliation cost" in line:
            c = int(line[line.index(":") + 2 : line.index("(") - 1])
            costs.append(c)
    weighted_costs = [costs[i] * weights[i] for i in range(len(costs))]
    return sum(weighted_costs)


# endregion Ranger Calls

# region getTotalRecCost

reconciliation_cache = {}
eccetera_time = 0.0
ranger_time = 0.0
read_time_gene = 0.0
read_time_fam = 0.0
write_time_spec = 0.0
write_time_gene = 0.0
total_eccetera_calls = 0
total_ranger_calls = 0


# Takes ranger_args, species tree, and gene trees path as input and returns the total reconciliation cost
def getTotalRecCost(
    dtl_args,
    spec_tree,
    gene_trees_path,
    gene_family_trees_path,
    weights,
    temp_dir,
    keep_temp=False,
):
    global eccetera_time
    global ranger_time
    global read_time_gene
    global read_time_fam
    global write_time_spec
    global write_time_gene
	
    # Create a cache key from inputs
    cache_key = (
        spec_tree,
        hash(gene_trees_path),
        tuple(weights) if weights else None,
    )

    # Check cache first
    if cache_key in reconciliation_cache:
        print(f"Cache hit for {cache_key}")
        return reconciliation_cache[cache_key]

    temp_path = os.path.join(temp_dir, f"TEMPFILE-{os.getpid()}")
    temp_fam_path = os.path.join(temp_dir, f"TEMPFAMILYFILE-{os.getpid()}")

    gene_family_indices = []
    start_time = time.time()
    if dtl_args["use_ecceTERA"]:
        with open(gene_family_trees_path, "r") as f:
            for line in f:
                gene_family_indices.append(int(line.strip()))
    read_time_fam += time.time() - start_time

    # Write species tree to temp file
    start_time = time.time()
    with open(temp_path, "w+") as f1, open(gene_trees_path, "r") as f2:
        f1.write(f"{spec_tree}\n")
        # Then copy contents of gene trees file ONLY if ecceTERA is NOT used
        # This is because ranger needs the gene trees in the same file as the species tree
        # and ecceTERA needs them in separate files
        if not dtl_args["use_ecceTERA"]:
            shutil.copyfileobj(f2, f1)
    write_time_spec += time.time() - start_time

    if len(gene_family_indices) and dtl_args["use_ecceTERA"]:
        gene_trees_all = []
        score = 0
        start_time = time.time()
        with open(gene_trees_path, "r") as f:
            for line in f:
                line = line.strip()
                if len(line) > 0:
                    if line[0] == "[":
                        end = line.index("]")
                        w = float(line[1:end])
                        weights.append(w)
                        gene_trees_all.append(line[end + 1 :])
                    else:
                        gene_trees_all.append(line)
        read_time_gene += time.time() - start_time

        tree_i = 0
        for family_i in gene_family_indices:
            family_trees = gene_trees_all[tree_i : tree_i + family_i]

            # Write family trees to temp file
            trees = False
            start_time = time.time()
            with open(temp_fam_path, "w+") as f:
                for tree in family_trees:
                    if tree != ";":
                        f.write(f"{tree}\n")
                        trees = True
            write_time_gene += time.time() - start_time

            # Call ecceTERA with family trees
            if trees:
                start_time = time.time()
                score += callEcceTERA(dtl_args, temp_path, temp_fam_path)
                eccetera_time += time.time() - start_time
            else:
                print(
                    f"Skipping family {tree_i}:{tree_i + family_i} with no trees. Score: {score}"
                )
                score += 0
            tree_i += family_i

        # print(f"Score: {score}")
        reconciliation_cache[cache_key] = score
        return score
    else:
        if dtl_args["use_ecceTERA"]:
            if len(weights) == 0:
                score = callEcceTERA(dtl_args, temp_path, gene_trees_path)
            else:
                score = callEcceTERAWeighted(
                    dtl_args, temp_path, gene_trees_path, weights
                )
        else:
            start_time = time.time()
            if len(weights) == 0:
                score = callRanger(dtl_args, temp_path)
            else:
                score = callRangerWeighted(dtl_args, temp_path, weights)
            ranger_time += time.time() - start_time

    if not keep_temp:
        if os.path.exists(temp_fam_path):
            os.remove(temp_fam_path)
        if os.path.exists(temp_path):
            os.remove(temp_path)
    # print(f"Worker: {os.getpid()}, {score}")
    reconciliation_cache[cache_key] = score
    return score


# endregion getTotalRecCost

# region Neighborhood Generation


# Takes starting tree as input and returns full SPR neighborhood as a list of newick strings
def genFullNeighborhood(template):
    neighborhood = [template.write()]

    order1 = [node for node in template.tree.iter_descendants()]
    for p_node in order1:
        w_name, sibling = template.prune(p_node)
        order2 = [node for node in template.tree.traverse()]
        for r_node in order2:
            if r_node != sibling:
                template.regraft(p_node, r_node, w_name)
                neighborhood.append(template.write())
                # revert back to outer loop template
                template.prune(p_node)
        # revert back to template tree
        template.regraft(p_node, sibling, w_name)

    return neighborhood


# Generates the neighborhood of possible regraft positions of p_node
def genRegraftNeighborhood(template, p_node_name):
    neighborhood = []
    p_node = template.n(p_node_name)
    w_name, sibling = template.prune(p_node)
    order2 = [node for node in template.tree.traverse()]
    for r_node in order2:
        if r_node != sibling:
            template.regraft(p_node, r_node, w_name)
            neighborhood.append(template.write())
            # revert back to outer loop template
            template.prune(p_node)
    # revert back to template tree
    template.regraft(p_node, sibling, w_name)
    return neighborhood


# endregion Neighborhood Generation

# region SPR Search


# Takes a template tree and ranger_args as input and returns the optimal trees and their score
# Single iteration of SPR search
def SPRsearch(
    template,
    start_score,
    dtl_args,
    gene_trees_path,
    gene_family_trees_path,
    num_cores,
    weights,
    temp_dir,
    keep_temp=False,
):
    neighborhood = genFullNeighborhood(template)
    # print(neighborhood)

    global getTotalRecCost_helper

    def getTotalRecCost_helper(spec_tree):
        score = getTotalRecCost(
            dtl_args,
            spec_tree,
            gene_trees_path,
            gene_family_trees_path,
            weights,
            temp_dir,
            keep_temp,
        )
        return score

    if num_cores > 1:
        p = multiprocessing.Pool(processes=num_cores)
        scores = p.map(getTotalRecCost_helper, neighborhood)
    else:
        scores = [
            getTotalRecCost(
                dtl_args,
                spec_tree,
                gene_trees_path,
                gene_family_trees_path,
                weights,
                temp_dir,
                keep_temp,
            )
            for spec_tree in neighborhood
        ]

    min_score = min(scores)
    optimal_trees = [
        neighborhood[i] for i in range(len(scores)) if scores[i] == min_score
    ]

    return optimal_trees, min_score


# Choosest highest score among all regraft positions of a single subtree if improvement, otherwise iterates to next subtree
def greedySPRsearch(
    template,
    start_score,
    dtl_args,
    gene_trees_path,
    gene_family_trees_path,
    num_cores,
    weights,
    temp_dir,
    keep_temp=False,
):
    global getTotalRecCost_helper

    def getTotalRecCost_helper(spec_tree):
        score = getTotalRecCost(
            dtl_args,
            spec_tree,
            gene_trees_path,
            gene_family_trees_path,
            weights,
            temp_dir,
            keep_temp,
        )
        return score

    n = len(template.tree)
    equiv_class = [template.write()]

    subtrees = [
        node.name
        for node in template.tree.iter_descendants()
        if len(node) < n - 1
    ]
    np.random.shuffle(subtrees)
    it = 1
    for p_node_name in subtrees:
        neighborhood = genRegraftNeighborhood(template, p_node_name)
        print(
            f"Regrafting node: {it}/{len(subtrees)} with {len(neighborhood)} neighbors"
        )
        it += 1

        if num_cores > 1:
            p = multiprocessing.Pool(processes=num_cores)
            scores = p.map(getTotalRecCost_helper, neighborhood)
        else:
            scores = [
                getTotalRecCost(
                    dtl_args,
                    spec_tree,
                    gene_trees_path,
                    gene_family_trees_path,
                    weights,
                    temp_dir,
                    keep_temp,
                )
                for spec_tree in neighborhood
            ]

        min_score = min(scores)
        optimal_trees = [
            neighborhood[i]
            for i in range(len(scores))
            if scores[i] == min_score
        ]
        if min_score < start_score:
            return optimal_trees, min_score
        elif min_score == start_score:
            equiv_class.extend(optimal_trees)

    return equiv_class, start_score


# endregion SPR Search
