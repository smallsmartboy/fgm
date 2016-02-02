#! /usr/bin/python
# -*-coding: utf-8 -*-
import os.path
from fgm import *

"""this file contains func reading graph from and writing graph to file"""


def read_eg_eattr(path, name):
    """read edge list and edge attribute from file"""
    file_path = os.path.join(path, name)
    with open(file_path) as f:
        lines = [line.rstrip('\n').strip() for line in f if line.rstrip('\n').strip()]
    edges = [get_edge(line) for line in lines]
    edge_attr = [get_edge_attr(line) for line in lines]
    return np.array(edges), np.array(edge_attr)


def read_node_nattr(path, name):
    """read nodes' attributes
    """
    file_path = os.path.join(path, name)
    with open(file_path) as f:
        lines = [line.rstrip('\n'.strip()) for line in f]
    lines = [line for line in lines if line]
    nods_attr = []
    [nods_attr.extend(get_node_attr(line)) for line in lines]
    nods_attr = [item.split(':') for item in nods_attr]
    nods_attr = zip(*nods_attr)
    na = []
    for line in nods_attr:
        na.append([int(item) for item in line])
    return np.array(na)


def get_node_attr(line):
    """get one node's attr
    :returns
        node_index,and a list of [feature_index, value] pairs
    """
    line = line.rstrip('\n').strip()
    items = line.split()
    node_index = items[0]
    feature_num = len(items[1:])
    return [node_index + ':' + item for item in items[1:]]


def get_edge(line):
    """get edge pair from a line"""
    line = line.rstrip('\n').strip()
    item = line.split()
    return [int(item[3].strip('"')), int(item[5].strip('"'))]


def get_edge_attr(line):
    """get edge attibute from a line"""
    line = line.rstrip('\n').strip()
    item = line.split()
    item =[it.strip('"') for it in item]
    if len(item) == 7:
        return [1, 0, 0]
    elif item[7] == "false":
        return [0, 1, 0]
    else:
        return [0, 0, 1]




