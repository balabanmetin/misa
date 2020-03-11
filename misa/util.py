

def index_edges(tree):
    counter=0
    for node in tree.traverse_postorder():
        node.edge_index = counter
        counter += 1
        node.valid = True
    tree.root.valid = False

def mrca_matrix(tree):
    '''Return a dictionary storing all pairwise MRCAs. ``M[u][v]`` = MRCA of nodes ``u`` and ``v``. Excludes ``M[u][u]`` because MRCA of node and itself is itself

    Returns:
        ``dict``: ``M[u][v]`` = MRCA of nodes ``u`` and ``v``
    '''
    M = dict()
    leaves_below = dict()
    for node in tree.traverse_postorder():
        leaves_below[node] = list()
        leaves_below[node].append(node)
        M[node] = dict()
        if not node.is_leaf():
            for i in range(len(node.children) - 1):
                for l1 in leaves_below[node.children[i]]:
                    M[l1][node] = node
                    M[node][l1] = node
                    leaves_below[node].append(l1)
                    for j in range(i + 1, len(node.children)):
                        for l2 in leaves_below[node.children[j]]:
                            M[l1][l2] = node
                            M[l2][l1] = node
            if len(node.children) != 1:
                for l2 in leaves_below[node.children[-1]]:
                    M[l2][node] = node
                    M[node][l2] = node
                    leaves_below[node].append(l2)
    return M