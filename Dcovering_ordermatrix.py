import os
os.chdir('C:\\Users\\ERIK\\Documents\\Python\\Bioinformatics Project\\')

from treeswift import read_tree_newick

tree_file = "Yeast_68leaves.txt"
tree = read_tree_newick(tree_file)

leaves = []
for label in tree.labels(internal=False):
    leaves.append(label)

M = tree.distance_matrix()
map_leaves = tree.label_to_node()

import random

def Dcovering_ordermatrix_random(D):
    random.shuffle(leaves)
    labels = [map_leaves[i] for i in leaves]
    subtreeslist = []
    spent_nodes = set()
    for i in range(len(labels)):
        if leaves[i] in spent_nodes:
            continue
        node1 = labels[i]
        spent_nodes.add(leaves[i])
        subtree = [leaves[i]]
        sorted_M = sorted(M[node1].items(), key=lambda x: x[1])
        for node2, distance in sorted_M:
            node2 = node2.get_label()
            if distance > D:
                break
            if node2 not in spent_nodes:
                subtree.append(node2)
                spent_nodes.add(node2)
        subtreeslist.append(subtree)
    return subtreeslist

def randomstarts(D, n):
    subtreeslist = Dcovering_ordermatrix_random(D)
    for i in range(n):
        s = Dcovering_ordermatrix_random(D)
        if len(s) < len(subtreeslist):
            subtreeslist = s
    return subtreeslist

import time

D = 0.5
n = 10000
start_time = time.time()
a =randomstarts(D, n)
print(len(a))

print("--- %s seconds ---" % (time.time() - start_time))