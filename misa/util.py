

def index_edges(tree):
    counter=0
    for node in tree.traverse_postorder():
        node.edge_index = counter
        counter += 1
        node.valid = True
    tree.root.valid = False
