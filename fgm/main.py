#!/usr/bin/python
# -*- encoding: utf-8 -*-
from graph import  graph
from gph import *
import collections

edges1 = read_gph_from_file('lan_20151102.txt',30000)
nodes1 = extract_nodes_from_edges(edges1)
m1, g1 = gph_eg2mat(node=nodes1, edge=edges1)

#edges2 = read_gph_from_file('20151103.txt', 60000)
#nodes2 = extract_nodes_from_edges(edges2)
#m2, g2 = gph_eg2mat(node=nodes2, edge=edges2)
ghp1 = graph('lan_20151102', nodes=nodes1, edges=edges1, g=g1, m=m1)
find_star_graph(nodes1, m1)