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

import itertools

def Dcovering_permutations(D):
    subtreeslist = [0]*len(leaves)
    for n in itertools.permutations(leaves):
        st = []
        spent_nodes = set()
        for leaf in n:
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
            st.append(subtree)
        if len(st) < len(subtreeslist):
            subtreeslist = st
    return subtreeslist


D = 1.878
a = Dcovering_permutations(D)
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