#!/usr/bin/python
# -*- encoding: utf-8 -*-
import os
import numpy as np


def get_dir_files(dir_name=None):
    """to gain the input file list

    :parameter
        dir_name:string
            directory name or file name

    :returns
        file_list: list
            file list contain all the absolute file path under the given dir
    """
    file_list = []
    if not os.path.exists(dir_name):
        return 0
    if os.path.isfile(dir_name):
        return file_list.append(dir_name)
    else:
        for strRoot, lsDir, lsFiles in os.walk(dir_name):
            for fn in lsFiles:
                file_list.append(os.path.join(strRoot, fn))
            if len(lsDir) > 0:
                for dn in lsDir:
                    file_list.append(os.path.join(strRoot, dn))
    return file_list


def read_gph_from_file(file_name=None, lines=1000):
    """
    :parameter
        file_name:string
    :returns

    """
    sep = '\t'
    edges = {}
    edge_num = 0
    i = 0
    with open(file_name, 'r') as f:
        for line in f:
            i += 1
            if i > lines:
                break
            line = line.strip().rstrip('\n')
            words = line.split(sep)
            if len(words) != 3:
                continue
            tmp = [words[0], words[1], words[2]]
            edges[edge_num] = tmp
            edge_num += 1
    print "read %d edges"% edge_num
    return edges


def extract_nodes_from_edges(edges):
    sep = '\t'
    nodes = set()
    [nodes.update([edge[0],edge[1]]) for edge in edges]
    nodes = list(nodes)
    print '%d nodes in total' % len(nodes)
    return nodes


def gph_eg2mat(node, edge):
    """
    :returns
        g:node-edge关联矩阵ni x mi
        h:node-edge增广矩阵ni x (mi + ni)
    """
    n = len(node)
    m = len(edge)

    g = np.zeros([n, m], dtype=np.int8)
    h = np.zeros([n, (n+m)], dtype=np.int8)
    for i in edge:
        sip = edge[i][0]
        dip = edge[i][1]
        node_index1 = node[sip]
        node_index2 = node[dip]
        h[node_index1][node_index2] = 1
        h[node_index2][node_index1] = 1
        g[node_index1][i] = 1
        g[node_index2][i] = 1
    return h, g


def gph_eg2feat(pt=np.array([]), eg=np.array([])):
    """calculate graph features
    :parameter
        pt      -  graph node, d x n
        eg      -  graph edge, 2 x m | []
    :returns
        Ptd     -  edge vector, d x m
        dists    -  distance, 1 x m
        angSs   -  angle (symmetric), 1 x m
        angAs   -  angle (asymmetric), 1 x m
    """
    if pt.size and eg.size:
        dists = np.sqrt(np.sum())
        ptd = np.array([])
        dists = np.array([])
        angss = np.array([])
        angas = np.array([])
    else:
        dists = np

# def calc_dists(nodes=np.array([]),edege=np.array([])):


def find_star_graph(nodes=np.array([]), g=np.array([])):
    """
    :parameter
        g:图的邻接矩阵
    :returns
        stars={center:[vertices]}
    """
    vsum = np.sum(g, 0)
    print g
    v = np.nonzero(vsum == 1)  # 找出所有星图的边缘顶点
    print v
    print len(v[0])
    esum = np.sum(g, 1)
    eg = np.nonzero(esum > 1)