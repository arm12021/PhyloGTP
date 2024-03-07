from ete3 import Tree
import subprocess
import shutil
import time
import multiprocessing
import os
import numpy as np

## Main tree class

class MasterTree():
    def __init__(self, tree=Tree()):
        self.tree = tree

    def print_structure(self):
        print(self.tree.get_ascii(show_internal=True))

    def print_newick(self, form=9):
        if form == 8:
            root_name = self.tree.name
            tree_str = self.tree.write(format=8)
            tree_str = tree_str[:-1] + root_name + ';'
            print(tree_str)
        else:
            print(self.tree.write(format=form))

    def write(self, form=8):
        if form == 8:
            root_name = self.tree.name
            tree_str = self.tree.write(format=8)
            tree_str = tree_str[:-1] + root_name + ';'
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
    #def prune(self, p_node_node):
        #p_node = self.tree&p_node_node
        assert not p_node.is_root(), 'pruned node cannot be the root'
        w_node = p_node.up
        w_node_name = w_node.name

        p_node.detach()
        sibling = w_node.children[0]

        if w_node.is_root():
            #sibling = w_node.get_children()[0].detach()
            sibling.detach()
            self.tree = sibling
        #else:
            #w_node.delete()
        w_node.delete()

        return w_node_name, sibling
        #return p_node, w_node_name

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
                node.name = 'n' + str(c)
                c += 1
    
## Helper functions
# Build a starting tree with a random topology 
def buildRandomStartingTree(names):
    t = Tree()
    t.populate(len(names), names_library=names)
    return t

# Makes system call to ranger and returns DTL score
def callRanger(ranger_args, temp_path):
    # make system call to ranger
    params = [ranger_args['path'], '-i', temp_path, '-D', ranger_args['D'], '-T', ranger_args['T'], '-L', ranger_args['L'], '-q', '-s']
    raw_out = subprocess.run(params, check=True, stdout=subprocess.PIPE, stderr = subprocess.DEVNULL, universal_newlines=True)
    #raw_out = subprocess.run(params, check=True, stdout=subprocess.PIPE, universal_newlines=True)
    ran_out = raw_out.stdout.strip()
    rec_cost = int(ran_out[ran_out.index(':')+2:ran_out.index('(')-1])
    return rec_cost

def callRangerWeighted(ranger_args, temp_path, weights):
    params = [ranger_args['path'], '-i', temp_path, '-D', ranger_args['D'], '-T', ranger_args['T'], '-L', ranger_args['L'], '-q']
    raw_out = subprocess.run(params, check=True, stdout=subprocess.PIPE, stderr = subprocess.DEVNULL, universal_newlines=True)
    ran_out = raw_out.stdout.strip().split('\n')
    costs = []
    for line in ran_out:
        if 'The minimum reconciliation cost' in line:
            c = int(line[line.index(':')+2:line.index('(')-1])
            costs.append(c)
    weighted_costs = [costs[i]*weights[i] for i in range(len(costs))]
    return sum(weighted_costs)

def getTotalRecCost(ranger_args, spec_tree, gene_trees_path, weights, temp_dir):
    temp_path = os.path.join(temp_dir, f'TEMPFILE-{os.getpid()}')
    # write species tree to temp file, then copy contents of gene trees file
    with open(temp_path, 'w+') as f1, open(gene_trees_path, 'r') as f2:
        f1.write(f'{spec_tree}\n')
        shutil.copyfileobj(f2, f1)
    if len(weights) == 0:
        score = callRanger(ranger_args, temp_path)
    else:
        score = callRangerWeighted(ranger_args, temp_path, weights)
    os.remove(temp_path)
    #print(f'Worker: {os.getpid()}, {score}')
    return score

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

# Single iteration of SPR search
def SPRsearch(template, start_score, ranger_args, gene_trees_path,  num_cores, weights, temp_dir):
    neighborhood = genFullNeighborhood(template)
    #print(neighborhood)

    global getTotalRecCost_helper
    def getTotalRecCost_helper(spec_tree):
        score = getTotalRecCost(ranger_args, spec_tree, gene_trees_path, weights, temp_dir)
        return score

    if num_cores > 1:
        p = multiprocessing.Pool(processes=num_cores)
        scores = p.map(getTotalRecCost_helper, neighborhood)
    else:
        scores = [getTotalRecCost(ranger_args, spec_tree, gene_trees_path, weights, temp_dir) for spec_tree in neighborhood]

    min_score = min(scores)
    optimal_trees = [neighborhood[i] for i in range(len(scores)) if scores[i] == min_score]
        
    return optimal_trees, min_score

# Choosest highest score among all regraft positions of a single subtree if improvement, otherwise iterates to next subtree
def greedySPRsearch(template, start_score, ranger_args, gene_trees_path, num_cores, weights, temp_dir):
    global getTotalRecCost_helper
    def getTotalRecCost_helper(spec_tree):
        score = getTotalRecCost(ranger_args, spec_tree, gene_trees_path, weights, temp_dir)
        return score
    
    n = len(template.tree)
    equiv_class = [template.write()]

    subtrees = [node.name for node in template.tree.iter_descendants() if len(node) < n-1]
    np.random.shuffle(subtrees)
    for p_node_name in subtrees:
        neighborhood = genRegraftNeighborhood(template, p_node_name)

        if num_cores > 1:
            p = multiprocessing.Pool(processes=num_cores)
            scores = p.map(getTotalRecCost_helper, neighborhood)
        else:
            scores = [getTotalRecCost(ranger_args, spec_tree, gene_trees_path, weights, temp_dir) for spec_tree in neighborhood]
        

        min_score = min(scores)
        optimal_trees = [neighborhood[i] for i in range(len(scores)) if scores[i] == min_score]
        if min_score < start_score:
            return optimal_trees, min_score
        elif min_score == start_score:
            equiv_class.extend(optimal_trees)
    
    return equiv_class, start_score


