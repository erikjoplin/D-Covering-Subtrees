import os
os.chdir('C:\\Users\\ERIK\\Documents\\Python\\Bioinformatics Project\\')

from treeswift import read_tree_newick

tree_file = "sugarbeets_10leaves.txt"
tree = read_tree_newick(tree_file)

leaves = []
for label in tree.labels(internal=False):
    leaves.append(label)
    
M = tree.distance_matrix()
map_leaves = tree.label_to_node()


import random

def Dcovering_random(D):
    random.shuffle(leaves)
    subtreeslist = []
    spent_nodes = set()
    for leaf in leaves:
        if leaf in spent_nodes:
            continue
        node1 = map_leaves[leaf]
        spent_nodes.add(leaf)
        subtree = [leaf]
        for node2, distance in M[node1].items():
            if distance > D:
                continue
            node2 = node2.get_label()
            if node2 not in spent_nodes:
                subtree.append(node2)
                spent_nodes.add(node2)
        subtreeslist.append(subtree)
    return subtreeslist

import time

def randomstarts(D, n):
    subtreeslist = Dcovering_random(D)
    for i in range(n):
        s = Dcovering_random(D)
        if len(s) < len(subtreeslist):
            subtreeslist = s
    return subtreeslist

D = 1.878
n = 60000

a = randomstarts(D, n)
print(len(a)) # should be length 3

def outputsubtrees(subtreeslist):
    subtrees = []
    for subtree in subtreeslist:
        st = tree.extract_tree_with(subtree)
        subtrees.append(st)
    return subtrees

b = outputsubtrees(a)

# with open('results.txt', 'w') as f:
#     for st in b:
#         print(st, file = f)