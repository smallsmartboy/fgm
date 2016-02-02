import copy
import os
import math
import sys
import itertools


def linkinfo(filename):
    handle = open(filename)
    info = set()
    for eachLine in handle:
        eachLine = eachLine.strip().split('\t')
        eachLine[2] = int(math.log10(int(eachLine[2])))+1
        eachLine = eachLine[0]+'\t'+eachLine[1]+'\t'+str(eachLine[2])
        info.add(eachLine)
    handle.close()
    return info


def commongraph1(filedir):
    filenames = os.listdir(filedir)
    commonlink = linkinfo(filedir+'/'+filenames[0])
    for file in filenames[1:]:
        links = linkinfo(filedir+'/'+file)
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
        nodepairs.append([edge[0],edge[1]])
        pairweight.append(edge[2])
    res = filedir+'\\res.txt'
    f = open(res, 'w')
    count = 0
    while nodepairs != []:
        count = count+1
        f.write('# commongraph'+str(count)+'\n')
        labelled = set(nodepairs[0])
        f.write(nodepairs[0][0]+'\t'+nodepairs[0][1]+'\t'+pairweight[0]+'\n')
        nodepairs.pop(0)
        pairweight.pop(0)
        labelledold = set()
        while labelled != labelledold:
            labelledold = copy.copy(labelled)
            de_index = []
            for i in xrange(len(nodepairs)):
                if labelled&set(nodepairs[i])!= set():
                    labelled |= set(nodepairs[i])
                    f.write(nodepairs[i][0]+'\t'+nodepairs[i][1]+'\t'+pairweight[i]+'\n')
                    de_index.append(i)
            for j in xrange(len(de_index)-1, -1, -1):
                 nodepairs.pop(de_index[j])
                 pairweight.pop(de_index[j])
    f.close()
 

def commongraph2(filedir,filenames,handle):
    commonlink = linkinfo(filedir+'/'+filenames[0])
    for file in filenames[1:]:
        links = linkinfo(filedir+'/'+file)
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
        nodepairs.append([edge[0],edge[1]])
        pairweight.append(edge[2])
    count = 0
    while nodepairs != []:
        count = count+1
        handle.write('# commongraph'+str(count)+'\n')
        labelled = set(nodepairs[0])
        handle.write(nodepairs[0][0]+'\t'+nodepairs[0][1]+'\t'+pairweight[0]+'\n')
        nodepairs.pop(0)
        pairweight.pop(0)
        labelledold = set()
        while labelled != labelledold:
            labelledold = copy.copy(labelled)
            de_index = []
            for i in xrange(len(nodepairs)):
                if labelled&set(nodepairs[i])!= set():
                    labelled |= set(nodepairs[i])
                    handle.write(nodepairs[i][0]+'\t'+nodepairs[i][1]+'\t'+pairweight[i]+'\n')
                    de_index.append(i)
            for j in xrange(len(de_index)-1,-1,-1):
                 nodepairs.pop(de_index[j])
                 pairweight.pop(de_index[j])


def hanlin_main():
    if len(sys.argv) ==3:
        commongraph1(sys.argv[1],sys.argv[2])
    else:
        filenumber = len(os.listdir(sys.argv[1]))
        selected_number = filenumber*float(sys.argv[3])
        if selected_number-int(selected_number)!=0:
            selected_number = int(selected_number)+1
        else:
            selected_number = int(selected_number)
        handle = open(sys.argv[2],'w')
        for files in itertools.combinations(os.listdir(sys.argv[1]),selected_number):
            commongraph2(sys.argv[1],files,handle)
        handle.close()

commongraph1("E:\\data\\graph\\lan_noscan\\201511")