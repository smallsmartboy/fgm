#!/usr/bin/python
# -*- encoding: utf-8 -*-


class GraphData:
    def __init__(self):
        __edge_x = []        #边一边的id
        __edge_y = []        #边另一边的id
        __node_labels = []   #节点标号
        __node_visibles = [] #节点是否可用
        __edge_labels = []   #边的集合标号
        __edge_visibles = [] #边是否可用

    def get_node_labels(self):
        return self.__node_labels

    def set_node_labels(self, node_labels):
        self.__node_labels = node_labels

    def get_node_visibles(self):
        return self.__node_visibles

    def set_node_visibles(self, node_visibles):
        self.__node_visibles = node_visibles

    def get_edge_labels(self):
        return self.__edge_labels

    def set_edge_labels(self, edge_labels):
        self.__edge_labels = edge_labels

    def get_edge_x(self):
        return self.__edge_x

    def set_edge_x(self, edge_x):
        self.__edge_x = edge_x

    def get_edge_y(self):
        return self.__edge_y

    def set_edge_y(self, edge_y):
        self.__edge_y = edge_y

    def get_edge_visibles(self):
        return self.__edge_visibles

    def set_edge_visibles(self, edge_visibles):
        self.__edge_visibles = edge_visibles

    def remove_not_fre_node_edge(self, freq_nodel_label=[], freq_edge_label=[], min_support=0):
        """
        根据点边的的频繁度一处图中不频繁的点和边
        """
        label = 0
        x = 0
        y = 0

        n = len(self.__node_labels)
        for i in xrange(n):
            label = self.__node_labels[i]
            if freq_nodel_label[label] < min_support:
                self.__node_visibles[i] = False

        n = len(self.__edge_labels)
        for i in xrange(n):
            label = self.__edge_labels[i]
            if freq_edge_label < min_support:
                self.__edge_visibles[i] = False
                continue

            x = self.__edge_x[i]
            y = self.__edge_y[i]
            if not (self.__node_visibles[x] and self.__node_visibles[y]):
                self.__edge_visibles[i] = False

    def relabel_by_rank(self, node_label2rank=[], edge_label2rank=[]):
        """
        根据标号排序重新对满足条件的点边重新编号
        """
        label = 0
        count = 0
        temp = 0
        oldId2New = []
        for i in xrange(len(self.__node_labels)):
            label = self.__node_labels[i]
            if self.__node_visibles[i]:
                self.__node_labels[i] = node_label2rank[label]
                oldId2New.append(count)
                count += 1

        for i in xrange(len(self.__edge_labels)):
            label = self.__edge_labels[i]
            if self.__edge_visibles[i]:
                self.__edge_labels[i] = edge_label2rank[label]
                temp = self.__edge_x[i]
                self.__edge_x[i] = oldId2New[temp]
                temp = self.__edge_y[i]
                self.__edge_y[i] = oldId2New[temp]
