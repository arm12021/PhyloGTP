# region Imports

import argparse
import os
import time
import shutil
import numpy as np

from pathlib import Path
from sys import platform
from ete3 import Tree

import utilities
import start

# endregion Imports

# region Arguments

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i",
    "--input-file",
    type=str,
    required=True,
    help="Path to input file containing input data of (optional) weights and gene trees. ",
)
parser.add_argument(
    "-s",
    "--start-tree-file",
    type=str,
    default=None,
    help="Path to species tree used as the starting point in the search. If none is provided, will use an additive taxon procedure.",
)
parser.add_argument(
    "-x",
    "--random-start",
    action="store_true",
    help="Toggle to start tree search from a random topology instead of the additive taxon procedure.",
)
parser.add_argument(
    "-o", "--output-dir", type=str, default="out", help="Path to output dir."
)
parser.add_argument(
    "-r",
    "--ranger-path",
    type=str,
    default=None,
    help="Path to Ranger-DTL-Fast executable. Defaults to the platform-compatible executable located in the binaries directory.",
)
parser.add_argument(
    "-D",
    "--dup-cost",
    type=str,
    default="2",
    help="Duplication cost used in Ranger.",
)
parser.add_argument(
    "-T",
    "--transfer-cost",
    type=str,
    default="4",
    help="Transfer cost used in Ranger.",
)
parser.add_argument(
    "-L", "--loss-cost", type=str, default="1", help="Loss cost used in Ranger."
)
parser.add_argument(
    "-t",
    "--threshold",
    type=float,
    default=0,
    help="Terminates loop after difference in scores falls below the given threshold.",
)
parser.add_argument(
    "-m",
    "--max-iter",
    type=int,
    default=-1,
    help="Max number of iterations before termination. Runs until convergence by default.",
)
parser.add_argument(
    "-f", "--full-SPR", action="store_true", help="Use the full SPR search."
)
parser.add_argument(
    "-c",
    "--num-cores",
    type=int,
    default=1,
    help="Number of CPUs to use in parallel.",
)
parser.add_argument(
    "-u",
    "--unrooted",
    action="store_true",
    help="Input gene trees are unrooted.",
)
parser.add_argument(
    "-k",
    "--keep-temp",
    action="store_true",
    default=False,
    help="Do not remove temporary files when finished.",
)


parser.add_argument(
    "-ue",
    "--use-ecceTERA",
    action="store_true",
    help="Use ecceTERA for reconciliation instead of RANGER-DTL.",
)
parser.add_argument(
    "-e",
    "--eccetera-path",
    type=str,
    default=None,
    help="Path to ecceTERA executable. Defaults to the platform-compatible executable located in the binaries directory.",
)
parser.add_argument(
    "-a",
    "--amalgamate",
    action="store_true",
    help="Amalgamate the given set of gene trees.",
)

args = parser.parse_args()

if not args.keep_temp:
    if os.path.isdir(args.output_dir):
        shutil.rmtree(args.output_dir)
os.makedirs(args.output_dir, exist_ok=True)
log_file_path = os.path.join(args.output_dir, "log.txt")

temp_dir = os.path.join(args.output_dir, "__TEMP__")
if not os.path.isdir(temp_dir):
    os.mkdir(temp_dir)
part_tree_dir = os.path.join(args.output_dir, "partial_trees")
if not os.path.isdir(part_tree_dir):
    os.mkdir(part_tree_dir)
gene_trees_path = os.path.join(temp_dir, f"GENETREES-{os.getpid()}")

if args.ranger_path is None:
    script_dir_path = Path(os.path.dirname(os.path.realpath(__file__)))
    parent_dir = script_dir_path.parent.absolute()
    bin_dir = os.path.join(parent_dir, "binaries")

    if platform == "linux":
        args.ranger_path = os.path.join(bin_dir, "Ranger-DTL-Fast.linux")
        args.eccetera_path = os.path.join(bin_dir, "ecceTERA_linux64")
    elif platform == "darwin":
        args.ranger_path = os.path.join(bin_dir, "Ranger-DTL-Fast.mac")
        args.eccetera_path = os.path.join(bin_dir, "ecceTERA_mac")
    elif platform == "win32":
        args.ranger_path = os.path.join(bin_dir, "Ranger-DTL-Fast.win")
        args.eccetera_path = os.path.join(bin_dir, "")
    else:
        print(
            "Uknown platform. Please specify the absolute path to the Ranger-DTL or ecceTERA executable."
        )


ranger_args = {
    "path": args.ranger_path,
    "D": args.dup_cost,
    "T": args.transfer_cost,
    "L": args.loss_cost,
    "use_ecceTERA": int(args.use_ecceTERA),
}

eccetera_args = {
    "path": args.eccetera_path,
    "D": args.dup_cost,
    "T": args.transfer_cost,
    "L": args.loss_cost,
    "unrooted": int(args.unrooted),
    "amalgamate": int(args.amalgamate),
    "use_ecceTERA": int(args.use_ecceTERA),
}


def check_args(args):
    if not os.path.isfile(args.input_file):
        print("Specified gene tree file does not exist.")
        return False
    if args.start_tree_file:
        if not os.path.isfile(args.start_tree_file):
            print("Specified starting tree file does not exist.")
            return False
    if not os.path.isfile(args.ranger_path) and not args.use_ecceTERA:
        print("Ranger executable file not found.")
        return False
    if not os.path.isfile(args.eccetera_path) and args.use_ecceTERA:
        print("ecceTERA executable file not found.")
        return False
    return True


# endregion Arguments

# region Main


def main(args):
    start_time = time.time()
    f = open(log_file_path, "w+")
    f.close()
    out_handle = log_file_path

    if args.use_ecceTERA:
        print("Using ecceTERA for reconciliation.\n")
        dtl_args = eccetera_args
    else:
        print("Using Ranger-DTL-Fast for reconciliation.\n")
        dtl_args = ranger_args

    # region Input Trees

    # Read in gene trees. Each line should be a separate gene tree.
    # If using weights, prepend each tree string with [w], where w is the weight.
    gene_trees = []
    weights = []
    with open(args.input_file) as f:
        for line in f:
            line = line.strip()
            if len(line) > 0:
                if line[0] == "[":
                    end = line.index("]")
                    w = float(line[1:end])
                    weights.append(w)
                    gene_trees.append(line[end + 1 :])
                else:
                    gene_trees.append(line)
    if len(weights) > 0 and len(weights) != len(gene_trees):
        print("Must include weight for every tree or no trees.")
        return

    # if args.use_ecceTERA:
    # Check input gene trees for ecceTERA compatibility
    # Gene trees must have unique leaf names (no duplicates) and no internal nodes with the same name as a leaf.
    with open(gene_trees_path, "w+") as f:
        for gene_tree_str in gene_trees:
            # Create a MasterTree object for each gene tree and set the node names.
            # This is necessary for the Ranger-DTL-Fast program to work correctly.
            # The node names are set to the leaf names of the tree.
            # The MasterTree object is then written to a file in Newick format.
            # The format is set to 1, which means that the tree is written in Newick format with branch lengths.
            mt = utilities.MasterTree(tree=Tree(gene_tree_str, format=1))

            if args.use_ecceTERA and mt.check_duplicate_leaves():
                print(
                    "Gene trees must have unique leaf names to be compatible with ecceTERA."
                )
                print(
                    f"Gene tree {gene_tree_str} has duplicate leaves {mt.duplicate_leaves()}."
                )
                print(f"Please fix this and try again.")
                return

            # Remove the leaf name suffixes if they exist.
            # This is necessary for the Ranger-DTL-Fast program to work correctly.
            if not args.use_ecceTERA:
                for leaf in mt.tree:
                    if "_" in leaf.name:
                        leaf.name = leaf.name[: leaf.name.index("_")]

            # Write the tree to a file in Newick format with branch lengths.
            # The format is set to 9, which means that the tree is written in Newick format with branch lengths and node names.
            if args.unrooted and not args.use_ecceTERA:
                f.write(f"[&U]{mt.write(form=9)}\n")
            else:
                f.write(f"{mt.write(form=9)}\n")

    # endregion Input Trees

    # region Start Tree

    # check if starting species tree file exists, if so read tree
    print("**Initializing Start Tree**\n")
    if args.start_tree_file:
        current_tree = utilities.MasterTree(
            tree=Tree(args.start_tree_file, format=1)
        )
        current_tree.set_node_names()
        starting_score = utilities.getTotalRecCost(
            dtl_args,
            current_tree.write(),
            gene_trees_path,
            weights,
            temp_dir,
            args.keep_temp,
        )
    else:
        # get taxa names from gene trees
        if args.random_start:
            taxa_ids = set()
            for gene_tree_str in gene_trees:
                t = Tree(gene_tree_str)
                for leaf in t:
                    if "_" in leaf.name:
                        taxa_ids.add(leaf.name[: leaf.name.index("_")])
                    else:
                        taxa_ids.add(leaf.name)
            current_tree = utilities.MasterTree(
                tree=start.buildRandomStartingTree(list(taxa_ids))
            )
            current_tree.set_node_names()
            starting_score = utilities.getTotalRecCost(
                ranger_args,
                current_tree.write(),
                gene_trees_path,
                weights,
                temp_dir,
                args.keep_temp,
            )
        else:
            current_tree, starting_score = start.additive_start_tree(
                gene_trees_path,
                dtl_args,
                weights,
                args.unrooted,
                temp_dir,
                args.num_cores,
                args.keep_temp,
            )
    current_tree.save(os.path.join(part_tree_dir, "0tree.nwk"))
    print(f"Start Tree: {current_tree.write(9)}")
    current_tree.print_structure()
    print("")

    # endregion Start Tree

    # return

    # region Tree Search

    # Check if the starting tree is a valid tree
    print("**Tree search**\n")
    print(f"  {'Iteration' + '':<10}{'Score' + '':<10}Num optimal trees\n")
    print(f"  {'0' + '':<10}{starting_score:<10}1")
    with open(out_handle, "a+") as f:
        f.write("**Tree search**\n\n")
        f.write("\tIteration\tScore\tNum optimal trees\tTime\n\n")
        cur_time = time.strftime("%X %x %Z")
        f.write(f"\t0\t{starting_score}\t1\t{cur_time}\n")

    prev_val = starting_score
    condition_reached = False
    it = 1
    if args.max_iter == 0:
        condition_reached = True
    while not condition_reached:
        if args.full_SPR:
            optimal_trees, min_score = utilities.SPRsearch(
                current_tree,
                prev_val,
                dtl_args,
                gene_trees_path,
                args.num_cores,
                weights,
                temp_dir,
                args.keep_temp,
            )
        else:
            optimal_trees, min_score = utilities.greedySPRsearch(
                current_tree,
                prev_val,
                dtl_args,
                gene_trees_path,
                args.num_cores,
                weights,
                temp_dir,
                args.keep_temp,
            )

        # Print/log the results and current time
        print(f"  {it:<10}{min_score:<10}{len(optimal_trees)}")
        with open(out_handle, "a+") as f:
            cur_time = time.strftime("%X %x %Z")
            f.write(f"\t{it}\t{min_score}\t{len(optimal_trees)}\t{cur_time}\n")

        # Check if the stopping condition is reached
        # 1. The difference in scores is less than the threshold
        if abs(min_score - prev_val) <= args.threshold:
            condition_reached = True
        # 2. The maximum number of iterations is reached
        elif args.max_iter != -1 and it >= args.max_iter:
            condition_reached = True

        # Update the current tree to the new optimal tree at a random index
        # This is done to avoid getting stuck in a local minimum
        # The optimal trees are stored in a list, and a random index is selected from that list.
        # The current tree is then set to the tree at that index.
        # The format is set to 8, which means that the tree is written in Newick format with branch lengths and node names.
        index = np.random.randint(len(optimal_trees))
        current_tree = utilities.MasterTree(
            tree=Tree(optimal_trees[index], format=8)
        )

        # Save the current tree to a file in Newick format with branch lengths
        # The format is set to 9, which means that the tree is written in Newick format with branch lengths and node names.
        current_tree.save(os.path.join(part_tree_dir, f"{it}tree.nwk"))
        prev_val = min_score
        it += 1

    # endregion Tree Search

    # region Final Tree

    current_tree.print_newick()
    current_tree.print_structure()
    print("")
    current_tree.save(os.path.join(args.output_dir, "final_tree.nwk"))

    if not args.keep_temp:
        shutil.rmtree(temp_dir)
    end_time = time.time() - start_time
    print(f"Total runtime: {end_time}")
    with open(out_handle, "a+") as f:
        f.write(f"\nTotal runtime: {end_time}")

    # endregion Final Tree


if __name__ == "__main__":
    if check_args(args):
        main(args)

# endregion Main
