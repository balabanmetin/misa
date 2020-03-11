#!/usr/bin/env python3
import json
import sys
from optparse import OptionParser

import treeswift as ts
from misa.util import index_edges,mrca_matrix
from misa.estimate_d3 import estimate_d3

from misa.Core import Core
from misa.Optimize import optimize_for_two
from misa.Distance import jc69
from misa.jutil import extended_newick
import re
import multiprocessing as mp


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-t", "--tree", dest="tree_fp",
                      help="path to the reference tree", metavar="FILE")
    parser.add_option("-d", "--distances", dest="dist_fp",
                      help="path to the table of observed distances", metavar="FILE")
    parser.add_option("-o", "--output", dest="output_fp",
                      help="path for the output jplace file",
                      metavar="FILE")
    parser.add_option("-l", "--model", dest="model_name", default="HAR",
                      help="name of the model (LIN, HAR, JAC)", metavar="MODEL")
    parser.add_option("-m", "--method", dest="method_name", default="FM",
                      help="name of the weighted least squares method (OLS, FM)", metavar="METHOD")
    # parser.add_option("-p", "--protein", dest="protein_seqs", action='store_true',
    #                  help="input sequences are protein sequences")
    parser.add_option("-T", "--threads", dest="num_thread", default="0",
                      help="number of cores used in placement. 0 to use all cores in the running machine",
                      metavar="NUMBER")
    parser.add_option("-i", "--iterations", dest="num_iterations", default="5000",
                      help="maximum number of iterations run by the optimizer",
                      metavar="NUMBER")
    parser.add_option("-k", "--kmer", dest="kmer_size", default="31",
                      help="size of the kmers",
                      metavar="NUMBER")
    (options, args) = parser.parse_args()
    tree_fp = options.tree_fp
    dist_fp = options.dist_fp
    output_fp = options.output_fp
    model_name = options.model_name
    method_name = options.method_name
    num_thread = int(options.num_thread)
    kmer_size = int(options.kmer_size)
    num_iterations = int(options.num_iterations)
    if not num_thread:
        num_thread = mp.cpu_count()

    tree = ts.read_tree_newick(tree_fp) # "data/backbone.tree")

    #estimate_d3(tree)
    index_edges(tree)
    extended_newick_string = extended_newick(tree)

    f = open(dist_fp) # "data/dist.mat")

    def read_dismat(f):
        tags = list(re.split("\s+", f.readline().rstrip()))[1:]

        line = f.readline()
        dists = list(re.split("\s+", line.strip()))
        query_name = dists[0]
        obs_dist = dict(zip(tags, map(float, dists[1:])))
        return (query_name, None, obs_dist)


    query_name, _ , obs_dist = read_dismat(f)
    f.close()

    core = Core(tree)
    # each leaf has an 'id' and each node has 'edge_index'. don't confuse them
    ind_key_obs = {n.id:(2*obs_dist[n.label])/(obs_dist[n.label]+1) for n in core.tree.traverse_postorder(internal=False)}

    mrca = mrca_matrix(tree)
    # print(mrca.keys())

    from misa.Method import OLS
    import math
    # meta = ts.read_tree_newick("yeast/VIN7/meta_backbone.tree")
    #
    # dm = meta.distance_matrix(leaf_labels=True)
    #
    # s1 = "Saccharomyces_cerevisiae"
    # s2 = "Saccharomyces_kudriavzevii"
    #
    # a_vec = {n.id:dm[s1][n.label] for n in core.tree.traverse_postorder(internal=False)}
    # [print([n.label]) for n in core.tree.traverse_postorder(internal=False)]
    #
    # b_vec = {n.id:dm[s2][n.label] for n in core.tree.traverse_postorder(internal=False)}
    #
    # print(a_vec)
    # print(b_vec)
    # x_test=[0]*(2*len(a_vec)+4)
    # for i in range(len(a_vec)):
    #     x_test[i]= math.exp(-kmer_size*a_vec[i])
    #     x_test[i+len(a_vec)]=math.exp(-kmer_size*b_vec[i])
    # x_test[-4]=0.05117223
    # x_test[-3]=0.05683918
    # x_test[-2]=0.01680891
    # x_test[-1]=0.07946303
    #
    # print(x_test)

    # orig_err = 0
    # for k in range(12):
    #     x = ind_key_obs[k] + (a_vec[k]-b_vec[k])/2
    #     y = ind_key_obs[k] - (a_vec[k]-b_vec[k])/2
    #     orig_err += (x-a_vec[k])**2 + (y - b_vec[k])**2
    #
    # print("ORIG ERR: ", orig_err)

    alpha=0.50

    def prepare_edge_pairs():
        for i in core.tree.traverse_postorder():
            if i == core.tree.root:
                continue
            for j in core.tree.traverse_postorder():
                if j == core.tree.root:
                    continue
                #if i.edge_index in [0] and  j.edge_index in [4,7]:
                #if i.edge_index <=  j.edge_index:
                if i != j and mrca[i][j] not in [i,j]:
                #    print(i.edge_index,j.edge_index)
                #if i.label and k.label and i.label == "Drosophila_pseudoobscura" and k.label == "Drosophila_sechellia":
                #if i.label and k.label and i.label == "Drosophila_persimilis" and k.label == "Drosophila_mojavensis":
                #if i.label and k.label and i.label == "Drosophila_persimilis" and k.label == "Drosophila_persimilis":
                    #dlist = [tree.distance_between(i, j)]
                    # if i.parent != tree.root:
                    #     dlist += [tree.distance_between(i.parent, j)]
                    # if j.parent != tree.root:
                    #     dlist += [tree.distance_between(i, j.parent)]
                    # d = max(dlist)
                    #
                    # mvec = [ind_key_obs[k] for k in sorted(ind_key_obs)]
                    # tup=(i, j, alpha, kmer_size, d, mvec)
                    # print(tup)
                    # print(OLS.f(x_test,*tup))
                    #exit(1)
                    yield (i, j, tree, ind_key_obs, model_name, method_name,num_iterations,kmer_size, alpha)
    all_edge_pairs = prepare_edge_pairs()

    #print(len(all_edge_pairs))
    pool = mp.Pool(num_thread)
    results = pool.starmap(optimize_for_two, all_edge_pairs)
    results_no_error = list(filter(lambda x: x[0] != None, results))

    res, b1, b2 = min(results_no_error, key=lambda x: x[0].fun)
    print(res.fun)
    #print(res.x)

    jplace = dict()
    jplace["tree"] = extended_newick_string
    jplace["metadata"] = {"invocation": " ".join(sys.argv)}
    jplace["fields"] = ["edge_num", "likelihood", "like_weight_ratio", "distal_length", "pendant_length"]
    jplace["version"] = 3
    jplace["placements"] = [{"p": [[b1.edge_index, res.fun, 1, res.x[-4], res.x[-3]]], "n": [query_name+"_1"]},
                            {"p": [[b2.edge_index, res.fun, 1, res.x[-2], res.x[-1]]], "n": [query_name+ "_2"]}]


    if output_fp:
        f = open(output_fp, "w")
    else:
        f = sys.stdout
    f.write(json.dumps(jplace, sort_keys=True, indent=4))
    f.close()

    # result = join_jplace(results)
    # result["tree"] = extended_newick_string
    # result["metadata"] = {"invocation":" ".join(sys.argv)}
    # result["fields"] = ["edge_num", "likelihood", "like_weight_ratio", "distal_length", "pendant_length"]
    # result["version"] = 3




            # if i.label and k.label and i.label == "Drosophila_pseudoobscura" and k.label == "Drosophila_sechellia":
            #     import pdb; pdb.set_trace()
            #count += 1

