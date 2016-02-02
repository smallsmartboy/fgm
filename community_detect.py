#!/usr/bin/python
# -*- coding: utf-8 -*-
from igraph import *
import graph
from graph_file import *
import matplotlib.pyplot as plt


def extract_nodes_from_edges(egs=[]):
    sep = '\t'
    nods = set()
    for edge in egs:
        ips = edge.split(sep)
        nods.add(ips[0])
        nods.add(ips[1])
    #print '%d nodes in total' % len(nods)
    return list(nods)


def construct_graph(name=None,nods=None, egs=None):
    sep = '\t'
    g = Graph()
    g['name'] = name
    g['layers'] = {}
    [g.add_vertex(nods.index(node)) for node in nods]
    #print "add nodes done!"
    eg_list = [tuple(nods.index(item) for item in eg.split(sep)) for eg in egs]
    g.add_edges(eg_list)
    #print "add edges done!"
    for n in g.vs:
        n['label'] = nods[n.index]
        n['layers'] = 100
        n['linkto'] = None
        #print "node atttibutes", n.attributes()
        #print n,n['label'],g.incident(n)
    #print "num 2 node's layer: ", g.vs[2]['layers']
    #print "edge 2: ", g.es[2]
    return g


def main(filepath):
    edges = read_edges(filepath)
    if not len(edges):
        return
    nodes = extract_nodes_from_edges(edges)
    g = construct_graph(name=filepath, nods=nodes, egs=edges)
    stars = graph.find_star(g)
    print "find {0} stars".format(len(stars))
    graph.mark(g)


def test(filepath="E:\\data\\graph\\lan_noscan\\201511\\20151102.txt"):
    file = filepath.split("\\")[-1]
    edges = read_edges(filepath)
    if not len(edges):
        return
    nodes = extract_nodes_from_edges(edges)
    g = construct_graph(name=filepath, nods=nodes, egs=edges)
    graph.mark(g)
    print file, len(nodes), g['layers']
    #print g.vs['layers']


def test2(filepath="E:\\data\\graph\\lan_noscan\\201511\\20151102.txt"):
    sep = '\t'
    edges = read_edges(filepath)
    #edge2csv(csv_file, edges)
    if not len(edges):
        return
    nodes = extract_nodes_from_edges(edges)
    g = construct_graph(name=filepath, nods=nodes, egs=edges)
    #graph.find_star(g)
    graph.mark(g)
    #layer2file(filepath, g)
    layers2file(filepath, g, -4)
    # print g['layers']
    # print len([n['name'] for n in g.vs if n['layers'] == flag])
    #out_layer = [n['name'] for n in g.vs if n['layers'] == 1]
    s = [layer for layer in g['layers']]
    print filepath, len(nodes), g['layers']
    #print sep.join([n['label'] for n in g.vs if n['layers'] == s[-2]])
    return g['layers']

if __name__ == '__main__':
    test2()
    print "here"
    filepath = "E:\\data\\graph\\lan_noscan\\201511\\20151105.txt"
    dirpath = "E:\\data\\graph\\lan_noscan\\201511"
    filelist = get_dir_files(dirpath)
    """
    for fi in filelist:
        name = fi.split('.')[0].split('\\')[-1]
        # main(fi)
        test(fi)
        layers = test2(fi)
        if layers:
            vs = [layers[k] for k in layers]
            plt.plot(range(len(vs[:-1])), vs[:-1], url=name)
            plt.text(0.5, (vs[1]+vs[2])/2, name)
    plt.show()
    """




