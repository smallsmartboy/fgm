#/usr/bin/python
# -*- coding: utf-8 -*-
import igraph
import os

"""this module handles

"""


def read_edges(graph_file=None):
    """reading the original graph file
    each line contains a edge like: label label weight
    :returns
    edges: list, each item present a edge
    """
    sep = '\t'
    edges = set()
    with open(graph_file) as f:
        for line in f:
            line = line.rstrip('\n').strip()
            items = line.split('\t')
            if len(items) < 3:
                print 'items wrong', items
            tmp = sep.join([items[0], items[1]])
            edges.add(tmp)
    # print "read %d edges" % len(edges)
    return list(edges)


def nodes2file(g=igraph.Graph, filepath=None):
    sep='\t'
    if not isinstance(g, igraph.Graph):
        print "input should be a igraph.Graph"
        return None
    with open(filepath) as f:
        pass


def get_dir_files(dir_name=None):
    """to gain the input file list

    :parameter
        dir_name:string
            directory name or file name

    :returns
        file_list: list
            file list contain all the absolute file path under the given dir
    """
    import os

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


def edge2csv(filepath=None, edges=None):
    sep = ','
    root = filepath.split(".")[0]
    csv_file = root + '.csv'
    f = open(csv_file, 'w')
    for line in edges:
        line = line.rstrip('\n')
        items = line.split('\t')
        f.write(sep.join(items))
        f.write('\n')
    f.close()


def layer2file(filepath=None, g=igraph.Graph, layer=-1):
    sep = '\t'
    layers = len(g['layers'])
    root = filepath.split('.')[0]
    if layer < 0:
        layer += layers
        layer += 1
    outfile = root + "_layer_" + str(layer) + "_flag_" + '.txt'
    edges = [e for e in g.es if e['layers'] == layer]
    f = open(outfile, 'w')
    for e in edges:
        line = sep.join([g.vs[e.source]['label'], g.vs[e.target]['label']])
        f.write(line)
        f.write('\n')
    f.close()


def layers2file(filepath=None, g=igraph.Graph, layer=-1):
    sep = ','
    layers = len(g['layers'])
    name = filepath.split('\\')[-1].split('.')[0]
    root = "E:\\data\\graph\\csv\\"
    if layer < 0:
        layer += layers
    outfile = root + name + "_gt_layer_flag_" + str(layer) + ".csv"

    edges = [e for e in g.es if e['layers'] >= layer]
    f = open(outfile, 'w')
    for e in edges:
        line = sep.join([g.vs[e.source]['label'], g.vs[e.target]['label']])
        f.write(line)
        f.write('\n')
    f.close()





