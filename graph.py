#!/usr/bin/python
# -*- coding: utf-8 -*-



class StarGraph:
    def __init__(self, g=igraph.Graph, cen=None, ver=None):
        if isinstance(g, igraph.Graph):
            self.g = g
        else:
            print "the first argument should be a igraph.Graph Object"
        self.center = cen
        self.vertices = ver


class ChainGraph:
    def __init__(self, g):
        self.g = g
        self.main_chain = []


def is_star(g=igraph.Graph()):
    if not g.is_connected():
        return False


def find_star(g=igraph.Graph):
    """find all the star graph in the given graph
    :parameter
        g: a Graph Object
    :returns
        stars:
    """
    v_candidate = set()
    cen_candidate = set()
    stars = []

    if not isinstance(g, igraph.Graph):
        print "input should be a igraph.Graph"
        return None

    for v in g.vs:
        v_adj = g.neighbors(v)
        if len(v_adj) <= 1:
            v_candidate.add(v['name'])
            cen_candidate.update(v_adj)
    print "v_candidate: ",len(v_candidate)

    def is_center(vj):
        for vi in g.neighbors(vj):
            if vi not in v_candidate:
                return False
        return True

    for vi in cen_candidate:
        if is_center(vi):
            star_g = igraph.Graph()
            star_g.add_vertex(vi)
            [star_g.add_vertex(vj) for vj in g.neighbors(vi)]
            star_g['name'] = vi
            [star_g.add_edge(0, star_g.vs['name'].index(vj)) for vj in g.neighbors(vi)]
            stars.append(StarGraph(star_g, cen=vi, ver=g.neighbors(vi)))
    return stars


def find_chain(g=igraph.Graph):
    chain_cen_candidate = set()
    pass


def mark(g=igraph.Graph, flag=1):
    n = 0
    candidate = [i for i in g.vs['name'] if g.vs[i]['layers'] >= flag]
    for v in candidate:
        v_adj = [j for j in g.neighbors(v) if g.vs[j]['layers'] >= flag]
        if len(v_adj) <= flag:  # 是否是边缘点
            g.vs[v]['layers'] = flag

            if v_adj:
                g.vs[v]['linkto'] = v_adj[0]

            n += 1

    if n:
        g['layers'][flag] = n
        mark(g, flag+1)
    else:
        g['layers'][flag] = len([n['name'] for n in g.vs if n['layers'] > flag])
        for n in g.vs:
            if n['layers'] > flag:
                n['layers'] = flag
    for e in g.es:
        e['layers'] = min(g.vs[e.source]['layers'], g.vs[e.target]['layers'])
        #print e




