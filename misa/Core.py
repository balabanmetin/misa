

class Core:
    def __init__(self, tree):
        self.tree = tree
        index=0
        for leaf in self.tree.traverse_postorder(leaves=True, internal=False):
            leaf.id=index
            index += 1

        self.create_S_table()
        self.create_R_table()



    # def validate_edges(self, obs_dist):
    #     for node in self.tree.traverse_postorder():
    #         if node.is_leaf():
    #             if obs_dist[node.label] == -1.0:
    #                 node.valid = False
    #             else:
    #                 node.valid = True
    #         else:
    #             if sum([i.valid for i in node.children]) == 0:
    #                 node.valid = False
    #             else:
    #                 node.valid = True
    #     for node in self.tree.traverse_preorder():
    #         if node == self.tree.root:
    #             node.valid = False
    #         elif sum([i.valid for i in filter(lambda x: x != node, node.parent.children)]) == 0 \
    #                 and not node.parent.valid:
    #             node.valid = False

    def create_R_table(self):
        for node in self.tree.traverse_preorder():
            if node == self.tree.root:
                continue
            node.Rd={}
            for sibling in node.parent.children:
                if sibling != node:
                    node.Rd.update({k: v + sibling.edge_length for k, v in sibling.Sd.items()})
            if node.parent != self.tree.root:
                node.Rd.update({k: v + node.parent.edge_length for k, v in node.parent.Rd.items()})
        return

    def create_S_table(self):
        for node in self.tree.traverse_postorder():
            if node.is_leaf():
                node.Sd={node.id:0}
            else:
                node.Sd={}
                for child in node.children:
                    node.Sd.update({k:v+child.edge_length for k,v in child.Sd.items()})
        return





