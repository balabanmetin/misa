

def index_edges(tree):
    counter=0
    for node in tree.traverse_postorder():
        node.edge_index = counter
        node.desc = {node.edge_index}
        if not node.is_leaf():
            for c in node.children:
                node.desc.update(c.desc)
        counter += 1
        node.valid = True
    tree.root.valid = False
