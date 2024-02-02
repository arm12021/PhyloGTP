import argparse
import os
from pathlib import Path
from sys import platform
import time
from ete3 import Tree
import shutil

from utilities import *
from start import additive_start_tree

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input-file', type=str, required=True, help='Path to input file containing input data of (optional) weights and gene trees. ')
parser.add_argument('-s', '--start-tree-file', type=str, default=None, help='Path to species tree used as the starting point in the search. If none is provided, will use an additive taxon procedure.')
parser.add_argument('-x', '--random-start', action='store_true', help='Toggle to start tree search from a random topology instead of the additive taxon procedure.')
parser.add_argument('-o', '--output-dir', type=str, default='out', help='Path to output dir.')
parser.add_argument('-r', '--ranger-path', type=str, default=None, help='Path to Ranger-DTL-Fast executable. Defaults to the platform-compatible executable located in the binaries directory.')
parser.add_argument('-D', '--dup-cost', type=str, default='2', help='Duplication cost used in Ranger.')
parser.add_argument('-T', '--transfer-cost', type=str, default='3', help='Transfer cost used in Ranger.')
parser.add_argument('-L', '--loss-cost', type=str, default='1', help='Loss cost used in Ranger.')
parser.add_argument('-t', '--threshold', type=float, default=0, help='Terminates loop after difference in scores falls below the given threshold.')
parser.add_argument('-m', '--max-iter', type=int, default=-1, help='Max number of iterations before termination. Runs until convergence by default.')
parser.add_argument('-f', '--full-SPR', action='store_true', help='Use the full SPR search.')
parser.add_argument('-c', '--num-cores', type=int, default=1, help='Number of CPUs to use in parallel.')
parser.add_argument('-u', '--unrooted', action='store_true', help='Input gene trees are unrooted.')
args = parser.parse_args()

if os.path.isdir(args.output_dir):
    shutil.rmtree(args.output_dir)
os.makedirs(args.output_dir)
log_file_path = os.path.join(args.output_dir, 'log.txt')

temp_dir = os.path.join(args.output_dir, '__TEMP__')
if not os.path.isdir(temp_dir):
    os.mkdir(temp_dir)
part_tree_dir = os.path.join(args.output_dir, 'partial_trees')
if not os.path.isdir(part_tree_dir):
    os.mkdir(part_tree_dir)
gene_trees_path = os.path.join(temp_dir, f'GENETREES-{os.getpid()}')

if args.ranger_path is None:
    script_dir_path = Path(os.path.dirname(os.path.realpath(__file__)))
    parent_dir = script_dir_path.parent.absolute()
    bin_dir = os.path.join(parent_dir, 'binaries')

    if platform == 'linux':
        args.ranger_path = os.path.join(bin_dir, 'Ranger-DTL-Fast.linux')
    elif platform == 'darwin':
        args.ranger_path = os.path.join(bin_dir, 'Ranger-DTL-Fast.mac')
    elif platform == 'win32':
        args.ranger_path = os.path.join(bin_dir, 'Ranger-DTL-Fast.win')
    else:
        print('Uknown platform. Please specify the absolute path to the Ranger-DTL executable.')




ranger_args = {
    'path': args.ranger_path,
    'D': args.dup_cost,
    'T': args.transfer_cost,
    'L': args.loss_cost
}

def check_args(args):
    if not os.path.isfile(args.input_file):
        print('Specified gene tree file does not exist.')
        return False
    if args.start_tree_file:
        if not os.path.isfile(args.start_tree_file):
            print('Specified starting tree file does not exist.')
            return False
    if not os.path.isfile(args.ranger_path):
        print('Ranger executable file not found.')
        return False
    return True

def main(args):
    start_time = time.time()
    f = open(log_file_path, 'w+')
    f.close()
    out_handle = log_file_path

    # Read in gene trees. Each line should be a separate gene tree.
    # If using weights, prepend each tree string with [w], where w is the weight. 
    gene_trees = []
    weights = []
    with open(args.input_file) as f:
        for line in f:
            line = line.strip()
            if len(line) > 0:
                if line[0] == '[':
                    end = line.index(']')
                    w = float(line[1:end])
                    weights.append(w)
                    gene_trees.append(line[end+1:])
                else:
                    gene_trees.append(line)
    if len(weights) > 0 and len(weights) != len(gene_trees):
        print('Must include weight for every tree or no trees.')
        return
    with open(gene_trees_path, 'w+') as f:
        for gene_tree_str in gene_trees:
            mt = MasterTree(tree=Tree(gene_tree_str, format=1))
            for leaf in mt.tree:
                if '_' in leaf.name:
                    leaf.name = leaf.name[:leaf.name.index('_')]
            if args.unrooted:
                f.write(f'[&U]{mt.write(form=9)}\n')
            else:
                f.write(f'{mt.write(form=9)}\n')

    # check if starting species tree file exists, if so read tree
    print('**Initializing Start Tree**\n')
    if args.start_tree_file:
        current_tree = MasterTree(tree=Tree(args.start_tree_file, format=1))
        starting_score = getTotalRecCost(ranger_args, current_tree.write(), gene_trees_path, weights, temp_dir)
    else:
        # get taxa names from gene trees
        if args.random_start:
            taxa_ids = set()
            for gene_tree_str in gene_trees:
                t = Tree(gene_tree_str)
                for leaf in t:
                    if '_' in leaf.name:
                        taxa_ids.add(leaf.name[:leaf.name.index('_')])
                    else:
                        taxa_ids.add(leaf.name)
            current_tree = MasterTree(tree=buildRandomStartingTree(list(taxa_ids)))
            current_tree.set_node_names()
            starting_score = getTotalRecCost(ranger_args, current_tree.write(), gene_trees_path, weights, temp_dir)
        else:
            current_tree, starting_score = additive_start_tree(gene_trees_path, ranger_args, weights, args.unrooted, temp_dir, args.num_cores)
    current_tree.save(os.path.join(part_tree_dir, '0tree.nwk'), format=8)

    print('**Tree search**\n')
    print('\tIteration\tScore\tNum optimal trees\n')
    print(f'\t0\t{starting_score}\t1')
    with open(out_handle, 'a+') as f:
        f.write('**Tree search**\n\n')
        f.write('\tIteration\tScore\tNum optimal trees\tTime\n\n')
        cur_time = time.strftime('%X %x %Z')
        f.write(f'\t0\t{starting_score}\t1\t{cur_time}\n')
    
    prev_val = starting_score
    condition_reached = False
    it = 1
    if args.max_iter == 0:
        condition_reached = True
    while not condition_reached:
        if args.full_SPR:
            optimal_trees, min_score = SPRsearch(current_tree, prev_val, ranger_args, gene_trees_path, args.num_cores, weights, temp_dir)
        else:
            optimal_trees, min_score = greedySPRsearch(current_tree, prev_val, ranger_args, gene_trees_path, args.num_cores, weights, temp_dir)

        index = random.randint(0, len(optimal_trees)-1)
        print(f'\t{it}\t{min_score}\t{len(optimal_trees)}')
        with open(out_handle, 'a+') as f:
            cur_time = time.strftime('%X %x %Z')
            f.write(f'\t{it}\t{min_score}\t{len(optimal_trees)}\t{cur_time}\n')
        if abs(min_score - prev_val) <= args.threshold:
            condition_reached = True
        elif args.max_iter != -1 and it >= args.max_iter:
            condition_reached = True
        current_tree = MasterTree(tree=Tree(optimal_trees[index], format=8))
        current_tree.save(os.path.join(part_tree_dir, f'{it}tree.nwk'), format=8)
        prev_val = min_score
        it += 1

    print('')
    current_tree.print_newick()
    current_tree.save(os.path.join(args.output_dir, 'final_tree.nwk'), format=8)
    #current_tree.save(log_file_path)

    shutil.rmtree(temp_dir)
    end_time = time.time() - start_time
    print(f'Total runtime: {end_time}')
    with open(out_handle, 'a+') as f:
        f.write(f'\nTotal runtime: {end_time}')

if __name__ == '__main__':
    if check_args(args):
        main(args)