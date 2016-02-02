#!/usr/bin/python
# -*- coding: utf-8 -*-
"""频繁子图挖掘，每个文件存储一个图，每一行为sip \t dip \t info格式
"""
import copy
import os
import math
import sys
import time
import itertools


def link_info(filename):
    """"读入文件，每一行文件 sip dip access_times
        weights=log10(access_times)+1
        :returns
            edges    -- the list of edges
    """
    sep = '\t'
    edges = set()
    with open(filename) as f:
        for line in f:
            line = line.strip().split('\t')
            line[2] = math.log10(int(line[3])) + 1
            edge = sep.join(line)
            edges.add(edge)
    return edges


def common_graph1(dir):
    """读取目录下所有文件，
    """
    filenames = get_dir_files(dir)
    print filenames
    commonlink = link_info(filenames[0])
    print "read the first done!"
    for fn in filenames[1:]:
        links = link_info(fn)
        commonlink &= links
    commonlink = sorted(list(commonlink))
    edges = []
    for link in commonlink:
        linkitem = link.split('\t')
        edges.append(linkitem)
    edges.sort()
    nodepairs = []
    pairweight = []
    for edge in edges:
        nodepairs.append([edge[0], edge[1]])
        pairweight.append(edge[2])
    res = dir+'\\res.txt'
    f = open(res, 'w')
    count = 0
    while not nodepairs:
        count += 1
        print "reading %dth file" % count
        f.write('# commongraph' + str(count) + '\n')
        labelled = set(nodepairs[0])
        f.write(nodepairs[0][0] + '\t' + nodepairs[0][1] + '\t' + pairweight[0] + '\n')
        nodepairs.pop(0)
        pairweight.pop(0)
        labelledold = set()
        while labelled != labelledold:
            labelledold = copy.copy(labelled)
            de_index = []
            for i in xrange(len(nodepairs)):
                if labelled & set(nodepairs[i]) != set():
                    labelled |= set(nodepairs[i])
                    f.write(nodepairs[i][0] + '\t' + nodepairs[i][1] + '\t' + pairweight[i] + '\n')
                    de_index.append(i)
            for j in xrange(len(de_index) - 1, -1, -1):
                nodepairs.pop(de_index[j])
                pairweight.pop(de_index[j])
    f.close()


def commongraph2(filedir, filenames, handle):
    commonlink = link_info(filedir + '/' + filenames[0])
    for file in filenames[1:]:
        links = link_info(filedir + '/' + file)
        commonlink &= links
    commonlink = sorted(list(commonlink))
    edges = []
    for link in commonlink:
        linkitem = link.split('\t')
        edges.append(linkitem)
    edges.sort()
    nodepairs = []
    pairweight = []
    for edge in edges:
        nodepairs.append([edge[0], edge[1]])
        pairweight.append(edge[2])
    count = 0
    while not nodepairs:
        count += 1
        handle.write('# commongraph' + str(count) + '\n')
        labelled = set(nodepairs[0])
        handle.write(nodepairs[0][0] + '\t' + nodepairs[0][1] + '\t' + pairweight[0] + '\n')
        nodepairs.pop(0)
        pairweight.pop(0)
        labelledold = set()
        while labelled != labelledold:
            labelledold = copy.copy(labelled)
            de_index = []
            for i in xrange(len(nodepairs)):
                if labelled & set(nodepairs[i]) != set():
                    labelled |= set(nodepairs[i])
                    handle.write(nodepairs[i][0] + '\t' + nodepairs[i][1] + '\t' + pairweight[i] + '\n')
                    de_index.append(i)
            for j in xrange(len(de_index) - 1, -1, -1):
                nodepairs.pop(de_index[j])
                pairweight.pop(de_index[j])


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
"""
def hanlin_main():
    if len(sys.argv) == 3:
        common_graph1(sys.argv[1], sys.argv[2])
    else:
        filenumber = len(os.listdir(sys.argv[1]))
        selected_number = filenumber * float(sys.argv[3])
        if selected_number - int(selected_number) != 0:
            selected_number = int(selected_number) + 1
        else:
            selected_number = int(selected_number)
        handle = open(sys.argv[2], 'w')
        for files in itertools.combinations(os.listdir(sys.argv[1]), selected_number):
            commongraph2(sys.argv[1], files, handle)
        handle.close()
"""


def main():
    print 'begin ..'
    common_graph1("E:\\data\\graph\\lan\\201511")
if __name__ == "__main__":
    main()