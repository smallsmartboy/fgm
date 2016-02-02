#!/usr/bin/python
#-*- encoding: utf-8 -*-

无向图因式分解图匹配算法
使用格式：fgmU(KP, KQ, Ct, gphs, asgT, pars{8}{:})
输入：
    KP：两图顶点相似性矩阵，n1 x n2
    KQ：两图边相似性矩阵，  m1 x m2
    Ct：两图的顶点对应矩阵，n1 x n2
    gphs：用于匹配的图，
        G：点-边邻接矩阵，ni x mi
        H：增量点-边邻接矩阵，ni x (mi + ni)
        asgI：真实任务，可以为空
    par：参数,字典
        nAlp:alpha{100}
        nItMa：每一步的最大迭代次数
        nHst ：FW算法修改的历史节点
        ip   ：是否利用IPFP来提高算法y or n
        thAlp：alpha的阈值，决定从什么时候启动CCCP，默认为0
        deb  ：是否调试，y or n
        idxAlp：需要追踪细节的alphas索引

输出：
    asg：
    alg：算法名“fgmU”
    X  ：两图的顶点对应矩阵n1 x n2
    acc：匹配度
    obj：



