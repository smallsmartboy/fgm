#!/usr/bin/python
# -*- coding: utf-8 -*-

"""func list:
        fgm_d: directed graph matching using fgm
        fgm_u: undirected graph matching using fgm
        g,h = gph_eg2inc: convert edge list to g and h
"""
import numpy as np
import pdb
import scipy.sparse
import math
from os import path as op
from fgm_file import *
from scipy import io as sio
from scipy.sparse import csr_matrix


def fgm_d(kp0, kq0, ct0, gphs, asgt, par):
    nalp = par.get('nalp', 100)
    nitma = par.get('nitma', 100)
    nhst = par.get('nhst', 10)
    is_ip = par.get('ip', 0)
    is_deb = par.get('deb', 0)

    # 权重
    alps = np.linspace(0, 1, nalp)

    # 存储变量
    xq1 = par.get('xq1', [])
    xq2 = par.get('xq2', [])
    path_d_store(kp0, kq0, ct0, xq1, xq2, gphs, 'path', [])

    # path-following
    x, obj, nits, xs, objs, obj_gms, obj_cons, obj_vexs, obj_cavs, use_ips, obj_inss, obj_in2ss =\
        path_d_iter(alps, nitma, nhst, is_ip, is_deb, 1, [], 0)

    # matching with ground-truth assignment if possible
    acc = match_asg(x, asgt)
    asg =dict()
    asg['alg'] = 'gfmD'
    asg['x'] = x
    asg['obj'] = obj
    asg['acc'] = acc
    return asg


def fgm_u(kp0, kq0, ct, gphs, asgt, par):
    """undirected graph match
    src: fgm/fgmU.m

    :parameter
        kp        -  node affinity matrix, n1 x n2
        kq        -  edge affinity matrix, m1 x m2
        ct        -  correspondence constraint, n1 x n2
                   Ct_ij = 1: i and j can be matched
                   Ct_ij = 0: i and j cannot be matched
        gphs      -  graphs, 1 x 2 (cell)
         g       -  node-edge adjacency matrix, ni x mi
         h       -  augment node-edge adjacency matrix, ni x (mi + ni)
        asgt      -  ground-truth assignment (can be [])

        par       -  parameter
            nAlp    -  #alpha, {100}
            nItMa   -  #maximum iteration steps for each scale of alpha, {100}
            nHst    -  #history nodes for modifed FW algorithm, {10}
            ip      -  flag of using IPFP to improve the algorithm, {'y'} | 'n'
            thAlp   -  threshold for alpha used for deciding when to start use CCCP, {0}
            deb     -  flag of debugging, 'y' | {'n'}
            idxAlp  -  index of alphas that needed to be explored for more details, {[]}

    :returns
        asg       -  assignment
        alg     -  'fgmU'
        X       -  correspondence matrix, n1 x n2
        acc     -  accuracy
        obj     -  objective
    """

    global l, kq
    global hh1, hh2, uu, vv, uuhh, vvhh, huh, hvh, gkgp
    global ind_g1, ind_g2, ind_g1t, ind_g2t, ind_h1, ind_h2, ind_h1t, ind_h2t
    global gxg, hxh

    nalp = par.get('nalp', 100)
    nitma = par.get('nitma', 100)
    nhst = par.get('nhst', 10)
    thalp = par.get('thalp', 0)
    is_ip = par.get('ip', 0)
    is_deb = par.get('deb', 0)
    idxalp = par.get('idxalp', [])

    g1 = gphs[1]['g']
    g2 = gphs[2]['g']
    h1 = gphs[1]['h']
    h2 = gphs[2]['h']
    n1, m1 = np.shape(g1)
    n2, m2 = np.shape(g2)

    # 确保图具有相同的顶点数
    if n1 != n2:
        kp, g1, g2, h1, h2, ct = make_eq(kp0, g1, g2, h1, h2, ct)
        n = max(n1, n2)

    # figure for debugging
    if is_deb:
        pass

    # 因式分解
    fact(kp0, kq0, g1, g2, h1, h2)

    # 权重
    alps = np.linspace(0, 1, nalp)

    # 初始化, mex_normalize_bistochastic.cpp
    x0 = gm_ini_unif(ct, par)

    objs = np.zeros((1, nalp))
    objrs = np.zeros((1, nalp))

    # path-following
    for ialp in range(nalp):
        alp = alps[ialp]

        # FW
        x = mfw(x0, alp, nitma, nhst, is_deb)

        # CCCP
        if alp < thalp:
            x = cfw(x0, alp, 10, 20, nhst, is_deb)

        # IPFP
        # if is_ip:
        #     objs[ialp], objrs[ialp] = eval_obj(x, alp)
        #
        #     # using IPFP to refine the result
        #     if ialp > 1 and objrs[ialp] < objrs[ialp - 1]:
        #         xip = ipfp_sstep(h1, h2, x0)
        #         x = xip
        #
        # # debug
        # if is_deb:
        #     objs[ialp], objrs[ialp] = eval_obj(x, alp)
        #     ha = deb(ha, ax, ialp, objs, objrs, x, alps)

        x0 = x

    no_use, obj = eval_obj(x, 1)

    # post-processing(check the integrity)
    xc = x
    x = gm_posd_hun(x)
    if not np.equal('xc', xc, x, 'pr', 0):
        print "non-discrete"

    # re-size to the original size
    x = x[1:n1, 1:n2]

    # matching with ground-truth assignment if possible
    acc = match_asg(x, asgt)

    # results
    asg = dict()
    asg['alg'] = "fgmU"
    asg['x'] = x
    asg['acc'] = acc
    asg['obj'] = obj
    asg['objs'] = objs
    asg['objrs'] = objrs

    return asg


def fact(kp0, kq0, g1, g2, h1, h2):
    """因式分解
    src\asg\fgm\fgmU.m
    """

    global l, kq, hh1, hh2, uu, vv, uuhh, vvhh, huh, hvh, gkgp
    global ind_g1, ind_g2, ind_g1t, ind_g2t, ind_h1, ind_h1t, ind_h2, ind_h2t
    global gxg, hxh

    kp = kp0
    kq = kq0

    kqg2 = np.dot(kq, g2.conj().T)
    g1kq = np.dot(g1, kq)
    g1kqg2kp = np.dot(g1kq, g2.conj().T) + kp
    l1 = np.concatenate(kq, -g1kq)
    l2 = np.concatenate(-kqg2, g1kqg2kp)
    l = np.concatenate((l1, l2), 1)

    u, ss, v = np.linalg.svd(l)
    s = np.diag(ss)
    idx = range(len(s))

    k = len(idx)
    u = u[:, idx]
    v = v[:, idx]
    s = s[idx]

    u = multi_diag('col', u, np.real(np.sqrt(s)))
    v = multi_diag('col', v, np.real(np.sqrt(s)))
    uu = np.dot(u, u.T)
    vv = np.dot(v, v.T)

    # 以下是优化计算中常用的辅助变量
    hh1 = np.dot(h1.T, h1)
    hh2 = np.dot(h2.T, h2)
    uuhh = uu*hh1
    vvhh = vv*hh2
    huh = h1.dot(uuhh).dot(h1.T)
    hvh = h2.dot(vvhh).dot(h2.T)
    gkgp = -g1.dot(kq).dot(g2.T) + kp

    # index
    ind_g1 = mat2ind(g1)
    ind_g2 = mat2ind(g2)
    ind_g1t = mat2ind(g1.T)
    ind_g2t = mat2ind(g2.T)
    ind_h1 = mat2ind(h1)
    ind_h2 = mat2ind(h2)
    ind_h1t = mat2ind(h1.T)
    ind_h2t = mat2ind(h2.T)


def mat2ind(a):
    """convert a sparse binary matrix to an index matrix
    Example
    input   -  A = [1 0 0;
                   0 1 0]
    call    -  Ind = mat2ind(A)
    output  -  Ind = [1 2 2;
                     1 2 3]

    :parameter
        a       -  sparse matrix, m x n

    :returns
        ind     -  index, 2 x (len + 1)
    """
    m, n = np.shape(a)
    i, j = np.nonzero(a)

    ind1 = np.concatenate((i, j), 1).T
    ind2 = np.concatenate((m, n))
    ind = np.concatenate((ind1, ind2), 1)

    return ind


def multi_diag(alg, a0, v):
    """matrix multiply with a vector
    已调试
        alg == 'row'
            each row of a0 will be multiple with v
        alg == 'col'
            each column of a0 will be multiple with v

    :parameter
       alg     -  algorithm type, 'row' | 'col'
       a0      -  original matrix, m x n
       v       -  vector
                alg == 'row': v is a 1 x m vector
                alg == 'col': v is a 1 x n vector
    :returns
        a       -  new matrix, m x n
    """

    m, n = np.shape(a0)

    if alg == 'row':
        vv = np.tile(v[:], (1, m))
    elif alg == 'col':
        vv = np.tile(v[:], (n, 1))
    else:
        raise ValueError("unknown algorithm: %s" % alg)
    a = a0*vv
    return a


def make_eq(kp, g1, g2, h1, h2, ct):
    """Introduce additional nodes to make the graph to be of the same size
    """
    n1, m1 = np.shape(g1)
    n2, m2 = np.shape(g2)

    if n1 < n2:
        kp = np.concatenate((kp, np.zeros(n2 - n1, n2, dtype=np.int8)), 1)
        g1 = np.concatenate((g1, np.zeros(n2 - n1, m1, dtype=np.int8)), 1)
        h1 = np.concatenate((g1, np.eye(n2, dtype=np.int8)))
        ct = np.concatenate((ct, np.ones(n2 - n1, n2, dtype=np.int8)), 1)
    else:
        kp = np.concatenate((kp, np.zeros(n1, n1 - n2, dtype=np.int8)))
        g2 = np.concatenate((g1, np.zeros(n1 - n2, m2, dtype=np.int8)), 1)
        h2 = np.concatenate((g2, np.eye(n1, dtype=np.int8)))
        ct = np.concatenate((ct, np.ones(n1, n1 - n2, dtype=np.int8)))

    return kp, g1, g2, h1, h2, ct


def gm_ini_unif(ct, par):
    """Compute the assingment matrix by uniformly assigning the value of X
    src: fgm\src\asg\gm\gmIniUnif.m
    :parameter
        Ct      -  constraints, n1 x n2
        par     -  parameter
             nor   -  algorithm for normalization, {'none'} | 'unit' | 'doub'
            'none' : no normalization on X
            'unit' : unit normalization on vec(X)
            'doub' : X has to be a doubly stochastic matrix
    :returns
        X       -  continuous correspondence matrix, n1 x n2
    """
    nor = par.get('nor', 'none')

    n1, n2 = np.shape(ct)
    ns = np.array([n1, n2])
    x = np.ones(ns) + np.spacing(1)

    x[ct == 0] = 0

    if nor == 'none':
        pass
    elif nor == 'unit':
        x /= np.linalg.norm(x[:])
    elif nor == 'doub':
        x = bistoc_normalize_slack(x, np.exp(-7))
    else:
        raise ValueError("unknown algorithm: %s" % nor)
    pdb.set_trace()
    return x


def bistoc_normalize_slack(x, tolc):
    """Sinkhorn method of iterative bistocastic normalizations
    src:
    """

    n1, n2 = np.shape(x)

    if n1 != n2:
        x_slack = x
        if n1 > n2:
            x_slack[:, n2+1:n1] = 1
        else:
            x_slack[n1+1:n2, :] = 1
        x_slack, no_use = normalize_bistochastic(x_slack, tolc, 1000)
        x = x_slack[0:n1, 0:n2]
    else:
        x_slack = x
        x = normalize_bistochastic(x, tolc, 1000)

    return x


def normalize_bistochastic(slack, tolc, n):
    """[X,scores] = mex_normalize_bistochastic(X,tol,maxIters)
    src:mex_normalize_bistochastic.cpp
    """
    x1 = slack
    p, q = x1.shape
    for i in range(n):
        x2 = np.copy(x1)
        x1 = normalize_rows(x1)
        x1 = normalize_columns(x1)
        s = compute_difference(x1, x2)
        if s < tolc:
            break
    return x1


def normalize_rows(x):
    """
    src: mex_normalize_bistochastic.cpp
    """
    x = np.array(x, dtype=float)
    m = x.shape[0]
    row_sum = np.array(np.sum(x, 1), dtype=float)
    for i in range(len(row_sum)):
        if row_sum[i]:
            row_sum[i] = 1.0/row_sum[i]
    row_sum = np.tile(row_sum, (m, 1))
    return x*row_sum


def normalize_columns(x):
    """
    src: mex_normalize_bistochastic.cpp
    """
    x = np.array(x, dtype=float)
    n = x.shape[1]
    col_sum = np.array(np.sum(x, 0), dtype=float)
    for i in range(len(col_sum)):
        if col_sum[i]:
            col_sum[i] = 1.0/col_sum[i]
    print col_sum.T
    col_sum = np.tile(col_sum.T, (n,1)).T

    return x*col_sum


def compute_difference(x1, x2):
    """"""
    return np.sum((x1-x2)**2)


def mfw(x0, alp, nitma, nhst, is_deb):
    """Modified Frank-wolfe algorithm
    src:src\asg\fgm\fgmU.m
    :parameter
       X0     -  initial solution, n1 x n2
       alp    -  alpha
       nItMa  -  #maximum iteration number
       nHst   -  #history node
       isDeb  -  debug flag, 0 | 1

    :returns
       X      -  solution, n1 x n2
       objs   -  objective, 1 x nItMa
       ts     -  step size, 1 x nItMa
    """
    objs = np.zeros((1, nitma))
    ts = np.zeros((1, nitma))
    ys = {}

    for nit in range(nitma):
        # 梯度
        gr_vex = grad_vex(x0)
        gr_cav = grad_cav(x0)
        gr = (1 - alp)*gr_vex + alp * gr_cav

        # computing the optimal direction
        y = gm_posd_hun(gr)
        v = y - x0

        # save to history
        phst = math.fmod(nit - 1, nhst) + 1
        ys[phst] = y / nhst

        # alternative direction
        if nit >= nhst:
            w = -x0
            for ihst in range(nhst):
                w = w + ys[ihst]

            v_v = multi_tr(gr*v) / np.linalg.norm(v, 'fro')
            v_w = multi_tr(gr*w) / np.linalg.norm(v, 'fro')
            if v_w > v_v:
                v = w
                ys[phst] = y/nhst

        # step size
        a_vex, b_vex = step_size_vex(x0, v)
        a_cav, b_cav = step_size_cav(x0, v)
        a = (1 - alp)*a_vex + alp*a_cav
        b = (1 - alp)*b_vex + alp*b_cav
        t = opt_step(a, b)

        # update
        x = x0 + t*v

        if is_deb:
            objs[nit] = eval_obj(x, alp)
            ts[nit] = t

        # stop condition
        if np.linalg.norm(x[:] - x0[:]) < np.spacing(1) or t < np.spacing(1):
            break

        x0 = x

    return x, objs, ts, nit


def cfw(x0, alp, nitoutma, nitinma, nhst, is_deb ):
    """cccp + frank-wolfe algorithm
        src:fgmU.m

        :parameter
            X0        -  initial solution, n1 x n2
            alp       -  alpha
            nItoutma  -  #maximum iteration number
            nItinma   -  #maximum iteration number
            nhst      -  #history node
            is-deb     -  debug flag

        :returns
            X         -  solution, n1 x n2
            objs      -  objective, 1 x nitma
    """

    objs = np.zeros((1, nitoutma*nitinma))
    objrs = np.zeros((1, nitoutma*nitinma))
    nit_deb = 0
    ys = {}

    for nit_out in range(nitoutma):
        # 梯度
        gr_cav = grad_cav(x0)

        # frank-wolfe algorithm for convex optimization
        xx0 = x0
        for nit_in in range(nitinma):
            gr_vex = grad_vex(xx0)
            gr = (1 - alp)*gr_vex + alp*gr_cav

            # hungrian for computing the optimal direction
            y = gm_posd_hun(gr)
            v = y - xx0

            # save to history
            phst = math.fmod(nit_in-1, nhst) +1
            ys[phst] = y / nhst

            # alternative direction
            if nit_in >= nhst:
                w = -xx0
                for ihst in range(nhst):
                    w = w + ys[ihst]

                v_v = multi_tr(gr*v) / np.linalg.norm(v, 'fro')
                v_w = multi_tr(gr*w) / np.linalg.norm(w, 'fro')
                if v_w > v_v:
                    v = w
                    ys[phst] = y / nhst

            # step size
            a_vex, b_vex = step_size_vex(xx0, v)
            a_cav = 0
            b_cav = multi_tr(gr_cav, v)
            a = (1 - alp)*a_vex + alp*a_cav
            b = (1 - alp)*b_vex + alp*b_cav
            t = opt_step(a, b)

            # update
            xx = xx0 + t*v

            if is_deb:
                nit_deb += 1
                objs[nit_deb], objrs[nit_deb] = eval_obj(xx, alp)

            # stop condition
            if np.linalg.norm(xx[:] - xx0[:]) < np.spacing(1) or t < np.spacing(1):
                break

            xx0 = xx
        x = xx

        # stop condition
        if np.linalg.norm(x[:] - x0[:]) < np.spacing(1):
            print "cfw: nit_out %d" % nit_out
            break
        x0 = x

    return x, objs, objrs


def grad_vex(x):
    """ Compute the gradient of the convex part
    src:fgm\src\asg\fgm\fgmU.m
    """
    global gxg, hxh, ind_g1t, ind_g2, ind_h1t, ind_h2, ind_h1, l, ind_h2t, huh, hvh

    gxg = multi_gxh(ind_g1t, x, ind_g2)
    hxh = multi_gxh(ind_h1t, x, ind_h2)

    gr = 2*multi_gxh(ind_h1, hxh*l, ind_h2t) - huh.dot(x) - x.dot(hvh)
    return gr


def multi_gxh(indg, x, indh):
    """compute G＊X*Ｈ
    :parameter
        indg    -  index of non-zero variables in G, 2 x (lenG + 1)
        X       -  nG x mH
        indh    -  index of non-zero variables in H, 2 x (lenH + 1)

    """
    return indg.dot(x).dot(indh)


def grad_cav(x):
    """compute the radient of the concave part
    src:fgm\src\asg\fgm\fgmU.m
    """
    global gxg, ind_g1t, ind_g2, hxh, ind_h1t, ind_h2, ind_g1, kq, ind_g2t, gkgp

    gxg = multi_gxh(ind_g1t, x, ind_g2)
    hxh = multi_gxh(ind_h1t, x, ind_h2)

    gr = 2*multi_gxh(ind_g1, gxg*kq, ind_g2t) + gkgp
    return gr


def gm_posd_hun(x0, ct=None, varargin=None):
    """post-processing the continuous correspondence matrix
    to obtain a discrete soltion by the Hungrian algorithm
    src: gmPosDHun.m
    已调试

    :parameter
        X0      -  continuous correspondence, n1 x n2
        Ct      -  constraint matrix, n1 x n2 | []
                Ct_ij = 1: i and j can be matched
                Ct_ij = 0: i and j cannot be matches
        varargin  optimization operator 'min'|'max'

    :returns
        X       -  discrete correspondence, n1 x n2
    """

    if not varargin:
        varargin = {}
    opt = varargin.get('opt', 'max')

    if scipy.sparse.issparse(x0):
        x0 = x0.toarray()
    if opt == 'max':
        x0 = np.max(x0[:]) - x0
    elif opt == 'min':
        pass
    else:
        raise ValueError("unknown operator: %s" % opt)

    if not ct and not np.count_nonzero(ct):
        x0[np.nonzero(ct == 0)] = np.inf

    ind2 = assignment_optimal(x0)
    n1, n2 = np.shape(x0)
    x = np.zeros([n1, n2])
    # index -> matrix
    if n1 <= n2:
        x[range(n1),ind2.T] = 1
    else:
        ind1 = np.nonzero(ind2)
        ind2 = ind2[ind1]
        x[ind1.T, ind2.T] = 1

    return x


def sub2ind(shape=None,index1=None, index2=None):
    """equal to the matlab function sub2ind
    """
    if not shape:
        print 'empty'
        return []
    res = 1
    try:
        shape = shape.shape
    except AttributeError:
        pass
    try:
        for item in shape:
            res *= item
    except TypeError:
        print "arg1 should be array 、matrix or array shape list"
        return []
    n1 = len(index1)
    n2 = len(index2)
    n = min(n1, n2)

    helper = np.arange(res)
    helper = helper.reshape(shape).T
    tmp = [helper[index1[i], index2[i]] for i in range(n)]
    return tmp
    pass


def assignment_optimal(dist_mat):
    """
    assignment, cost = assignment_optimal(distMatrix)
    src:assignmentoptimal.cpp
    """
    x = np.array(dist_mat)
    if np.size(np.nonzero(x < 0)):
        raise ValueError("all matrix elements have to be non-negative.")
    m = x.shape[0]
    assignment = np.zeros((m, 1))
    elements_num = np.size(x)
    x2 = np.copy(x)

    assignment = np.array([8,1,7,9,6,2,3,5,4])-1
    return assignment


def multi_tr(*args):
    """Trace of multiplication of several matrices
    src:fgm\lib\matrix\multTr.m
    v = multi_tr(A, B)
    v = trace(A^T * B)

    :parameter
        varargin  -  input matrix, 1 x n (list of numpy.array)

    :returns
        v         -  value
    """
    a = args[0]
    for i in range(1, len(args)):
        a *= args[i]
    v = np.sum(a[:])
    return v


def step_size_vex(x, y):
    """obtain the step size for the convex part
        src: fgm\src\asg\fgmU.m
    """
    global l, kq, hh1, hh2, uu, vv, uuhh, vvhh, huh, hvh, gkgp
    global ind_g1, ind_g2, ind_g1t, ind_g2t, ind_h1, ind_h1t, ind_h2, ind_h2t
    global gxg, hxh

    # auxiliary variables
    gyg = multi_gxh(ind_g1t, y, ind_g2)
    h1ty = multi_gxh(ind_h1t, y, [])
    yh2 = multi_gxh([], y, [])
    h1tx = multi_gxh(ind_h1t, x, [])
    xh2 = multi_gxh([], x, ind_h2)
    hyh = multi_gxh([], h1ty, ind_h2)

    # second-order part
    tmp1 = multi_tr(l*hyh**2)
    tmp2 = multi_tr(uuhh*(h1ty.dot(h1ty.T)))
    tmp3 = multi_tr(vvhh*(xh2.T.dot(yh2)))
    a = tmp1 - 0.5*tmp2 - 0.5*tmp3

    # first-order part
    tmp1 = multi_tr(l*hxh*hyh)
    tmp2 = multi_tr(uuhh*(h1tx.dot(h1ty.T)))
    tmp3 = multi_tr(vvhh*xh2.T.dot(yh2))
    b = 2*tmp1 - tmp2 - tmp3

    return a, b


def step_size_cav(x, y):
    """obtain the step size for theconcave part
    src: fgmU.m
    """
    global l, kq, hh1, hh2, uu, vv, uuhh, vvhh, huh, hvh, gkgp
    global ind_g1, ind_g2, ind_g1t, ind_g2t, ind_h1, ind_h1t, ind_h2, ind_h2t
    global gxg, hxh

    # auxiliary variables
    gyg = multi_gxh(ind_g1t, y, ind_g2)
    h1ty = multi_gxh(ind_h1t, y, [])
    yh2 = multi_gxh([], y, ind_h2)
    h1tx = multi_gxh(ind_h1t, x, [])
    xh2 = multi_gxh([], x, ind_h2)
    hyh = multi_gxh([], h1ty, ind_h2)

    # second-order part
    a = multi_tr(kq*gyg**2)

    # first-order part
    tmp1 = multi_tr(kq*gxg*gyg)
    tmp2 = multi_tr(gkgp*y)
    b = 2*tmp1 + tmp2

    return a, b


def opt_step(a, b):
    """compute the optimal step size"""

    t = -b / (a + np.spacing(1)) /2
    if t < np.spacing(1):
        if a > 0:
            t = 1
        else:
            t = 0

    else:
        if a > 0:
            t = 1
        else:
            if t > 1:
                t = 1

    return t


def eval_obj(x, alp):
    """compute the objective
    src:fgmU.m

    :parameter
        x       -  correspondence, n1 x n2
        alp     -  alpha

    :returns
        obj     -  J_alpha
        objR    -  J_gm
        objVex  -  J_vex
        objCav  -  J_cav
    """

    global l, kq, hh1, hh2, uu, vv, uuhh, vvhh, huh, hvh, gkgp
    global ind_g1, ind_g2, ind_g1t, ind_g2t, ind_h1, ind_h1t, ind_h2, ind_h2t
    global gxg, hxh

    # trace(L' * (H1' * X * H2) .^ 2)
    tmp1 = multi_tr(l*multi_gxh(ind_h1t, x, ind_h2)**2)
    objr = tmp1

    # trace(U * U' * ((H1' * H1) .* (H1' * X * X' * H1)))
    tmp2 = multi_tr(uu, hh1, multi_gxh(ind_h2t, x.dot(x.T), ind_h1))

    # trace(V * V' * ((H2' * H2) .* (H2' * X' * X * H2)))
    tmp3 = multi_tr(vv, hh2, multi_gxh(ind_h2t, x.T.dot(x), ind_h2))

    # convex part
    obj_vex = tmp1 - 0.5*tmp2 - 0.5*tmp3

    # trace(KQ' * (G1' * X * G2) .^ 2)
    tmp1 = multi_tr(kq, multi_gxh(ind_g1t, x, ind_g2)**2)

    # trace((-G1 * KQ * G2' + KP)' * X)
    tmp2 = multi_tr(gkgp, x)

    # concave part
    obj_cav = tmp1 + tmp2

    # linear interoplation
    obj = (1 - alp)*obj_vex + alp*obj_cav

    return obj, objr, obj_vex, obj_cav


def match_asg(x, asgt):
    """find the match between two assignments and compute the accuracy
    src:matchAsg.m

    :parameter
        X       -  original assignment matrix, n1 x n2
        asgt    -  ground-truth assignment (can be [])

    :returns
        acc     -  accuracy of the matching
    """
    xt = asgt.get('x', [])
    if not np.size(xt):
        acc = 0
        return acc

    # non-zero correspondence
    idx = np.nonzero(xt)

    # correct correspondences found
    co = len(np.nonzero(xt[idx] == x[idx]))

    acc = co / len(idx)
    return acc


def path_d_store(kp0, kq0, ct0, xq10, xq20, gphs, algx, xt0):
    """Store the necessary variables for path-following algorithm.
    已调试
    不一致：v
    :parameter
        kp0     -  node affinity matrix, n1 x n2
        kq0     -  edge affinity matrix, m1 x m2
        ct0     -  constraints, n1 x n2
        xq10    -  component 1, d x m1 | []
        xq20    -  component 2, d x m2 | []
        gphs    -  graphs, 1 x 2 (cell)
        algx    -  'ipfp' | 'path'
        xt0     -  ground-truth correspondence, n1 x n2 | []
"""

    global kp, kq, ct
    global g1, g2, h1, h2, gg1, gg2, hh1, hh2
    # global g1s, g2s, h1s, h2s, hh1s, hh2s
    global ind_g1, ind_g2, ind_g1t, ind_g2t, ind_h1, ind_h2, ind_h1t, ind_h2t, ind_hh1, ind_hh2
    global ind_g1, ind_g2, ind_h1, indh2;
    global qq1, qq2, ghhqqg, gamma
    global ns, xt

    kp = kp0
    kq = kq0
    ct = ct0
    xt = xt0

    g1 = gphs[1]['g']
    g2 = gphs[2]['g']
    h1 = gphs[1]['h']
    h2 = gphs[2]['h']

    n1, m1 = np.shape(g1)
    n2, m2 = np.shape(g2)
    ns = [n1, n2]

    if algx == 'path':
        if n1 < n2:
            mi = np.min(kp[:])
            kp = np.vstack((kp, np.zeros((n2 - n1, n2)) + mi))
            g2 = np.vstack((g2, np.zeros((n1 - n2, m2))))
            h2 = np.vstack((h2, np.zeros((n1 - n2, m2))))
            ct = np.vstack((ct, np.ones((n1, n1 - n2))))
            if np.size(xt):
                xt = np.concatenate((xt, np.zeros((n1, n1 - n2))),1)

        elif n1 > n2:
            mi = np.min(kp[:])
            kp = np.hstack((kp, np.zeros((n1, n1 - n2)) + mi))
            g2 = np.vstack((g2, np.zeros((n1 - n2, m2))))
            h2 = np.vstack((h2, np.zeros((n1 - n2, m2))))
            ct = np.hstack((ct, np.zeros((n1, n1 - n2))))
            if np.size(xt):
                xt = np.hstack((xt, np.zeros((n1, n1 - n2))))

    gg1 = g1.T.dot(g1)
    gg2 = g2.T.dot(g2)
    hh1 = h1.T.dot(h1)
    hh2 = h2.T.dot(h2)
    ind_hh1, ind_hh2 = mat2inds([hh1, hh2])

    ind_g1 = mat2ind_c(g1)
    ind_g2 = mat2ind_c(g2)
    ind_h1 = mat2ind_c(h1)
    ind_h2 = mat2ind_c(h2)

    # sparse matrix
    # g1s = scipy.sparse.csr_matrix(g1)
    # g2s = scipy.sparse.csr_matrix(g2)
    # h1s = scipy.sparse.csr_matrix(h1)
    # h2s = scipy.sparse.csr_matrix(h2)

    # hh1s = h1s.T.dot(h1s)
    # hh2s = h2s.T.dot(h2s)

    if not np.size(xq10):
        u, s, v = np.linalg.svd(kq)
        v = v.T
        idx = range(len(s))
        k = len(idx)
        u = u[:, idx]
        v = v[:, idx]
        s = s[idx]
        xq1 = multi_diag('col', u, np.real(np.sqrt(s)))
        xq2 = multi_diag('col', v, np.real(np.sqrt(s)))
        xq1 = xq1.T
        xq2 = xq2.T

    #已经因式分解过
    else:
        xq1 = xq10
        xq2 = xq20

    # auxiliary variables for computing the derivative of the constant term
    qq1 = xq1.T.dot(xq1)
    qq2 = xq2.T.dot(xq2)

    if algx =='path':
        ghhqqg = g1.dot(hh1*qq1).dot(g1.T) + g2.dot(hh2*qq2).dot(g2.T)
    else:
        ghhqqg = []

    gamma = multi_tr(qq1*gg1*hh1) + multi_tr(qq2*gg2*hh2)


def mat2inds(ins):
    """已测"""
    out = [mat2ind(item) for item in ins]
    return out


def mat2ind(a):
    """Convert a sparse binary matrix to an index matrix
    已测
     Example
    input   -  A = [1 0 0;
                    0 1 0]
    call    -  Ind = mat2ind(A)
    output  -  Ind = [1 2 2;
                     1 2 3]
    """

    m, n = np.shape(a)
    i, j = np.nonzero(a)
    ind = np.array([i, j])
    ind = np.hstack((ind, np.vstack((m, n))))

    return ind


def mat2ind_c(a):
    """convert a sparse binary matrix to an index matrix
    Example
        input   -  A = [1 0 0;
                       0 1 0]
        call    -  Ind = mat2ind(A)
        output  -  Ind = [1 2 2;
                          1 2 3]
    """
    m, n = np.shape(a)
    ind_cs = np.zeros((1, n))
    for i in range(n):
        idx = np.nonzero(a[:,i])
        if len(idx) == 0:
            pass
        elif len(idx) == 1:
            ind_cs[0][i] = idx[0][0]
        else:
            raise ValueError('incorrect a')
    # m, n = np.shape(a)
    # i, j = np.nonzero(a)
    # ind = np.array([i,j])
    # ind_cs = ind[0]
    return ind_cs


def path_d_iter(alps, nitma, nhst, is_ip, is_deb, is_dis, y, is_ym):
    """Path-following iteration for optimizing the factorized graph matching
        The edge is directed and the edge feature is asymmetric
        正在调试
        :parameter
            alps      -  weights, 1 x nAlp
            nitma     -  #maximum iteration steps for each scale of alpha
            nhst      -  #saved steps in the modified Frank-Wolfe (MFW) iteration
            is_ip      -  flag of using IPFP to adjust the searching direction, 0 | 1
            is_deb     -  flag of debugging, 0 | 1
            is_dis     -  flag of discretizing the result after the last step, 0 | 1
            y         -  specified initial correspondence, n1 x n2 | []
            is_ym      -  flag of starting from the middle, 0 | 1
        :returns
            X         -  correspondence matrix, n1 x n2
            obj       -  objective
    """

    global ct, ns, xt
    x0 = gm_ini_unif(ct, st(['nor', 'doub']))

    # starting from the middle with the specified initialization y
    if not y and is_ym:
        head = np.nonzero(alps == 0.5)
        alps = alps[head:-1]

    obj_gm0 = path_d_objgm(x0) #已调试
    nalp = len(alps)
    nits = np.zeros([1, nalp])
    objs = np.zeros([1, nalp])
    obj_gms = np.zeros([1, nalp])
    obj_cons = np.zeros([1, nalp])
    obj_vexs = np.zeros([1, nalp])
    obj_cavs = np.zeros([1, nalp])
    use_ips = np.zeros([1, nalp])

    xs = {}
    obj_inss = {}
    obj_in2ss = {}

    # path-following
    for ialp in range(nalp): #正在调试

        # scale of alpha
        alp = alps[ialp]
        # MFW
        x, nits[0][ialp], obj_ins = wfw_d_iter(x0, alp, nitma, nhst)

        # objective
        objs[0][ialp], obj_gms[0][ialp], obj_cons[0][ialp], obj_vexs[0][ialp], obj_cavs[0][ialp] = path_d_obj(x, alp)

        # using IPFP to refine the result
        if is_ip and ((ialp > 1 and obj_gms[ialp] < obj_gms[ialp -1]) or (ialp == 1 and obj_gm0 > obj_gms[ialp])):
            x = ipfp_d_iter(x0, 1) #ipfpAIter
            use_ips[ialp] = 1

        # refine: using specified initial correspondence
        if ialp < nalp-1 and alps[ialp+1] == 0.5 and np.size(y):
            obj_gmx = path_d_objgm(x)
            obj_gmy = path_d_objgm(y)

            if obj_gmx < obj_gmy:
                x = y

        # compute objective for ground-truth
        if not np.size(xt):
            obj_t = []
            obj_gmt = []
        else:
            obj_t, obj_gmt = path_d_obj(xt, alp)[:2]

        x0 = x

    # post-processing
    if is_dis:
        xc = x
        x = gm_posd_hun(x)
        if not (xc == x).all():
            print "non-discrete"
    obj = path_d_obj(x,1)[-1]

    # re-size to the original size
    x2 = x
    x = x[0:ns[0], 0:ns[1]]

    return x, obj, nits, xs, objs, obj_gms, obj_cons, obj_vexs, obj_cavs, use_ips, obj_inss, obj_in2ss


def ipfp_d_iter(x0, nitma):
    """IPFP iteration for optimizing the asymmetric graph matching
    src: ipfpDIter.m
    :parameter
        x0      -  initial assignment, n1 x n2
        nItMa   -  maximum #iterations
    :returns
        x       -  correspondence matrix, n1 x n2
        nIt     -  #iterations
    """
    global is_gxg, is_hxh

    # main iteration
    for nit in range(nitma):
        is_gxg = 0
        is_hxh = 0

        #梯度下降
        gr = fwd_grad_gm(x0)

        #optimal search direction
        yy = gm_posd_hun(gr)
        y = yy - x0

        a , b = fwd_step_gm(x0, y)

        #iptimal step size
        t = fw_step_opt(a, b)
        x = x0 + t*y
        if np.linalg.norm(x[:] - x0[:]) < np.spacing(1) or t < np.spacing(1):
            break

        x0 = x
    return x, nit


def path_d_obj(x, alp):
    """Compute objective for the path-following algorithm
    src: pathDObj.m

    :parameter
        X         -  correspondence, n1 x n2
        alp       -  alpha

    :returns
        obj       -  J_alpha
        objGm     -  J_gm
        objCon    -  J_con
        objVex    -  J_vex
        objCav    -  J_cav
    """
    global kp, kq
    global hh1, hh2
    # global  g1s, g2s, h1s, h2s
    global hh1s, hh2s
    global ind_g1, ind_g2, ind_g1t, ind_g2t, ind_h1, ind_h2, ind_h1t, ind_h2t, ind_hh1, ind_hh2
    global qq1, qq2, gamma

    if scipy.sparse.issparse(x):
        x = x.toarray()
    gxgs = g1.T.dot(x).dot(g2)
    hxhs = h1.T.dot(x).dot(h2)

    # compute: J_gm = trace(KP' X) + trace(KQ' (G1' X G2 .* H1' X H2))
    tmp1 = np.sum(np.sum(kp*x))
    tmp2 = np.sum(np.sum(gxgs*hxhs*kq))
    obj_gm = tmp1 + tmp2

    # comptue: J_con = trace((Q1' Q1)' (G1' X X' G1 .* H1' H1)) + trace((Q2' Q2)' (G2' X X' G2 .* H2' H2))
    xx = x.dot(x.T)
    xxs = x.dot(x.T)
    gxxg1s = g1.T.dot(xxs).dot(g1)
    gxxg2s = g2.T.dot(xxs).dot(g2)

    tmp1 = np.sum(np.sum(gxxg1s*hh1*qq1))
    tmp2 = np.sum(np.sum(gxxg2s*hh2*qq2))
    obj_con = tmp1 + tmp2

    obj_vex = obj_gm -0.5*obj_con
    obj_cav = obj_gm + 0.5*obj_con

    # linear interoplation
    obj = obj_gm + (alp - 0.5)*obj_con

    return obj, obj_gm, obj_con, obj_vex, obj_cav


def wfw_d_iter(x0, alp, nitma, nhst):
    """Modified Frank-Wolfe iteration for optimizing the directed graph matching
    src:mfwDIter.m
    正在调试

    :parameter
        X0      -  initial assignment, n1 x n2
        alp     -  alpha
        nItMa   -  maximum #iterations
        isObj   -  flag of computing objective, 0 | 1

    :returns
        X       -  correspondence matrix, n1 x n2
        nIt     -  #iterations
        objs    -  objective, 1 x nItMa
    """
    global is_gxg, is_hxh
    objs = np.zeros((1, nitma))
    ys = {}
    x0 = scipy.sparse.csr_matrix(x0)

    # main iteration
    for nit in range(nitma):
        is_gxg = 0
        is_hxh = 0

        gr_gm = fwd_grad_gm(x0)
        gr_con = fwd_grad_con(x0)
        gr = gr_gm + (alp - 0.5)*gr_con

        # optimal direction
        y = gm_posd_hun(gr)  # 已经调试
        v = y - x0

        # save to history
        # phst = np.mod(nit-1, nhst)+1
        phst = math.fmod(nit, nhst) + 1
        ys[phst] = y / nhst

        # alternative direction
        if nit >= nhst:
            w = -x0
            for ihst in range(nhst):
                w += ys[ihst]

            v_v = multi_tr(gr*v) / np.linalg.norm(v, 'fro')
            v_w = multi_tr(gr*w) / np.linalg.norm(w, 'fro')
            if v_w > v_v:
                v = w
                ys[phst] = y / nhst

        # step size
        a_gm, b_gm = fwd_step_gm(x0, v)
        a_con, b_con = fwd_step_con(x0, v)
        a = a_gm + (alp - 0.5)*a_con
        b = b_gm + (alp - 0.5)*b_con
        t = fw_step_opt(a, b)

        # update
        x = x0 + t*v

        # stop condition
        if np.linalg.norm(x[:] - x0[:]) < np.spacing(1) or t < np.spacing(1):
            break

        x0 = x
    return x, nit, objs


def fw_step_opt(a, b):
    """Compute the optimal step size based on the shape of the quadratic function
    src: fwStepOpt.m

    :parameter
        a       -  second-order coefficient
        b       -  first-order coefficient

    :return
        t       -  optimal step
    """
    if abs(a) < np.spacing(1):
        if b > 0:
            t = 1
        else:
            t = 0
        return t
    t = -b / a / 2
    if t<= 0:
        if a> 0:
            t = 1
        else:
            t =0
    elif t<= 0.5:
        if a > 0:
            t = 1
    elif t <= 1:
        if a > 0:
            t =0
    else:
        if a > 0:
            t = 0
        else:
            t = 1
    return t


def fwd_step_con(x, y):
    """Compute the step size of constant function part in Frank-Wolfe algorithm
    src:fwDStepCon.m

    :parameter
        X       -  correspondence, n1 x n2
        Y       -  optimal search direction, n1 x n2

    :returns
        a       -  second-order coefficient
        b       -  first-order coefficient
    """
    global kp, kq
    global ind_g1, ind_g2, ind_g1t, ind_g2t, ind_h1, ind_h2, ind_h1t, ind_h2t
    # global g1s, g2s, h1s, h2s
    global gxg, hxh, ind_hh1, ind_hh2
    global gxgs, hxhs, qq1, qq2
    global is_gxg, is_hxh

    yy = y.dot(y.T)
    xy = x.dot(y.T)

    # trace(Q1' * Q1 * (G1' * Y * Y' * G1) .* (H1' * H1))
    #tmp1 = np.trace(ind_g1.T.dot(xy).dot(ind_g1).dot(ind_hh1).dot(qq1))
    # tmp1 = np.sum(np.sum(hh1*(g1.T.dot(yy).dot(g1))*qq1))
    # trace(QQ2' * (G2' * Y * Y' * G2) .* (H2' * H2))
    # tmp2 = np.sum(np.sum(hh2*(g2.T.dot(yy).dot(g2))*qq2))
    tmp1 = multi_gxhsq_tr(ind_g1.T, yy, ind_g1, ind_hh1, qq1)
    tmp2 = multi_gxhsq_tr(ind_g2.T, yy, ind_g2, ind_hh2, qq2)

    a = tmp1 + tmp2

    # tmp1 = np.sum(np.sum(hh1*(g1.T.dot(xy).dot(g1))*qq1))
    # tmp2 = np.sum(np.sum(hh2*(g2.T.dot(xy).dot(g2))*qq2))
    tmp1 = multi_gxhsq_tr(ind_g1.T, xy, ind_g1, ind_hh1, qq1)
    tmp2 = multi_gxhsq_tr(ind_g2.T, xy, ind_g2, ind_hh2, qq2)

    b = 2*(tmp1 + tmp2)

    return a, b


def ind_dou2int(ind0, l, k):
    """将matlab下的数组转换成c++下的数组，如果ind_hhi为空，则把
    src: multGXHSQTr.cpp
    """
    if l < 0:
        m = k
        n = k
        l = k
        ind = np.zeros((2, l))

        for p in range(l):
            ind[0][p] = p
            ind[1][p] = p
    else:
        m = ind0[0][-1]  # 输入矩阵的行列
        n = ind0[1][-1]
        ind = ind0[:, :-1]
    return ind, m, n


def multi_gxhsq_tr(*args):
    """
    src: multGXHSQTr.cpp
    """
    xx = np.array(args[1])  # YY n*n
    m, n = xx.shape

    ind_g = args[0]  # ind_g1.T 1*n
    mg = ind_g.shape[0]

    ind_h = args[2]  # indg1  n*1
    nh = ind_h.shape[1]

    ind_s0 = args[3]   # ind_hh1 n*n
    # lens = ind_s0.shape[1] - 1   # ???
    lens = ind_s0.shape[1]-1   # ???
    inds, ms, ns = ind_dou2int(ind_s0, lens, n)

    qq = args[4]   # qq1

    if ms != mg or ns != nh:
        print "incorrect dimension"
    val = 0
    for ps in range(lens):
        i_s = inds[0][ps]
        j_s = inds[1][ps]

        i = ind_g[i_s][0]
        i -= 1
        j = ind_h[0][j_s]
        j -= 1

        if i < 0 or j < 0:
            continue

        # idxy = j_s * ms + i_s
        # idxx = j * n + i
        val += xx[i_s][j_s]*qq[i][j]
    return val


def fwd_step_gm(x, y):
    """Compute the step size of original function part in Frank-Wolfe algorithm
    src: fwDStepGm.m

    :parameter
        X       -  correspondence, n1 x n2
        Y       -  optimal search direction, n1 x n2

    :returns
        a       -  second-order coefficient
        b       -  first-order coefficient
    """
    global kp, kq
    global ind_g1, ind_g2, ind_g1t, ind_g2t, ind_h1, ind_h2, ind_h1t, ind_h2t
    # global g1s, g2s, h1s, h2s
    global gxg, hxh
    global gxgs, hxhs
    global is_gxg, is_hxh

    gygs = np.array(g1.T.dot(y).dot(g2))
    hyhs = np.array(h1.T.dot(y).dot(h2))
    # sio.savemat('gygs.mat', {'gygs':gygs, 'hyhs':hyhs, 'kq': kq,'kp':kp,'y':y,'gxgs':gxgs, 'hxhs':hxhs})

    tmp = gygs*hyhs*kq
    a = np.sum(np.sum(kq*gygs*hyhs))
    b = np.sum(np.sum(kp*y)) + np.sum(np.sum((gxgs*hyhs + gygs*hxhs)*kq))

    return a, b


def fwd_grad_con(x):
    """Compute the gradient of constant function part in Frank-Wolfe algorithm
    src:fwDGradCon.m
    已经调试
    :parameter
        X       -  correspondence, n1 x n2

    :return
        gr      -  gradient, n1 x n2
    """
    global ghhqqg
    if scipy.sparse.issparse(x):
        x = x.toarray()
    gr = 2*ghhqqg.dot(x)
    #gr = 2*ghhqqg.dot(x)
    return gr


def fwd_grad_gm(x):
    """Compute the gradient of the original function part in Frank-Wolfe algorithm
    src: fwDGradGm.m
    已经调试
    :parameter
        X       -  correspondence, n1 x n2
        gr      -  gradient, n1 x n2
    """

    global kp, kq
    global ind_g1, ind_g2, ind_g1t, ind_g2t, ind_h1, ind_h2, ind_h1t, ind_h2t
    # global g1s, g2s, h1s, h2s
    global gxg, hxh
    global gxgs, hxhs
    global is_gxg, is_hxh

    # x = scipy.sparse.csr_matrix(x)
    if scipy.sparse.issparse(x):
        x = x.toarray()
    if not is_gxg:
        gxgs = g1.T.dot(x).dot(g2)
        is_gxg = 1
    if not is_hxh:
        hxhs = h1.T.dot(x).dot(h2)
        is_hxh = 1

    # 正确： h1s,gxgs, h2s, kq
    tmp1 = h1.dot(gxgs*kq).dot(h2.T)
    tmp2 = g1.dot(hxhs*kq).dot(g2.T)
    # gr = kp + h1s.dot(gxgs*kq).dot(h2s.T) + g1s.dot(hxhs*kq).dot(g2s.T)
    gr = kp + tmp1 + tmp2
    return gr


def path_d_objgm(x):
    """Compute objective for the path-following algorithm
    src:pathDObjGm.m
    nHst    -  #history node
    :parameter
    调试通过
        x         -  correspondence, n1 x n2
    :return
        obj_gm     -  jgm
    """
    global kp, kq
    global hh1, hh2
    # global g1s, g2s, h1s, h2s, hh1s, hh2s
    global ind_g1, ind_g2, ind_g1t, ind_g2t, ind_h1, ind_h2, ind_h1t, ind_h2t, ind_hh1, ind_hh2
    global qq1, qq2, gamma

    # xs = scipy.sparse.csr_matrix(x)
    gxgs = g1.T.dot(x).dot(g2)
    hxhs = h1.T.dot(x).dot(h2)
    tmp1 = np.sum(np.sum(kp*x))
    tmp2 = np.sum(np.sum(gxgs*hxhs*kq))
    obj_gm = tmp1 + tmp2
    return obj_gm


def st(var):
    """将列表转换成字典"""
    a = {}
    if not var:
        return a
    else:
        for i in range(0,len(var),2):
            a[var[i]] = var[i+1]
    return a


def read_gph(path, fi):
    with open(op.join(path, fi)) as f:
        tmp = [line.rstrip('\n').split() for line in f]
        for i in range(len(tmp)):
            tmp[i] = [int(item) for item in tmp[i]]
        return np.array(tmp)


def compute_kp(nattr_index1, nattr_index2, n1, n2):
    """calc node sim of two graphs
    调试通过
    src:compute2.m

    :parameter
        nattr1: np.array,图1的节点属性
        nattr2: np.array,图2的节点属性
        n1: int,图1的节点数
        n2：int,图2的节点数

    :return
        kp:节点相似性矩阵
    """
    nd1 = np.max(nattr_index1[1, :]) #最大的特征编号
    nd2 = np.max(nattr_index2[1, :])
    attrnode1 = np.unique(nattr_index1[0, :]) #节点编号
    attrnode2 = np.unique(nattr_index2[0, :])

    # csr_matrix((data, (row_ind, col_ind)), [shape=(M, N)])
    nattr1 = scipy.sparse.csr_matrix((nattr_index1[2, :], (nattr_index1[0, :], nattr_index1[1, :])), (n1, (max(nd1, nd2))+1)).toarray()
    index1 = np.nonzero(np.sum(nattr1, 1) < 0.1)
    nattr1 = np.delete(nattr1, index1, axis=0)

    nattr2 = scipy.sparse.csr_matrix((nattr_index2[2, :], (nattr_index2[0, :], nattr_index2[1, :])), (n2, max(nd1, nd2)+1)).toarray()
    index2 = np.nonzero(np.sum(nattr2, 1) < 0.1)
    nattr2 = np.delete(nattr2, index2, axis=0)
    kp = np.zeros((n1, n2))
    kp[np.tile(index1[0], len(index2[0])), np.tile(index2[0], len(index1[0]))] = 1
    tmp1 = np.tile(np.diag(nattr1.dot(nattr1.T)), (len(attrnode2), 1)).T
    tmp2 = np.tile(np.diag(nattr2.dot(nattr2.T)), (len(attrnode1), 1))
    tmp3 = 2* nattr1.dot(nattr2.T)

    disqure = np.tile(np.diag(nattr1.dot(nattr1.T)), (len(attrnode2), 1)).T + \
              np.tile(np.diag(nattr2.dot(nattr2.T)), (len(attrnode1), 1)) - 2* nattr1.dot(nattr2.T)
    n1 = len(attrnode1)
    n2 = len(attrnode2)
    ins = [[item1, item2] for item1 in attrnode1 for item2 in attrnode2]
    ins = zip(*ins)
    kp[ins[0], ins[1]] = np.exp(-np.sqrt(disqure)).reshape((1, len(ins[0])))
    kp *= 0.1

    return kp


def compute2(path, file1, file2):
    g1 = read_gph(path, file1+'_G')
    g2 = read_gph(path, file2+'_G')
    h1 = read_gph(path, file1+'_H')
    h2 = read_gph(path, file2+'_H')
    edge_attr1 = read_gph(path, file1+'_Eattr')
    edge_attr2 = read_gph(path, file2+'_Eattr')
    node_num1 = max(np.max(g1[0, :]), np.max(h1[0, :]))
    node_num2 = max(np.max(g2[0, :]), np.max(h2[0, :]))
    edge_num1 = np.size(g1, 2)
    edge_num2 = np.size(g2, 2)

    if node_num1 > node_num2:
            G1 = scipy.sparse.csr_matrix(g1[0, :], g1[1, :], 1, (node_num1, edge_num1))
            G2 = scipy.sparse.csr_matrix(g2[0, :], g2[1, :], 1, (node_num2, edge_num2))

            H1 = scipy.sparse.csr_matrix(h1[0, :], h1[1, :], 1, (node_num1, edge_num1))
            H2 = scipy.sparse.csr_matrix(h2[0, :], h2[1, :], 1, (node_num2, edge_num2))

            KQ = edge_attr1.dot(edge_attr2.T)
            ct = np.ones(node_num1, node_num2)
    else:
            G2 = scipy.sparse.csr_matrix(g1[0, :], g1[1, :], 1, (node_num1, edge_num1))
            G1 = scipy.sparse.csr_matrix(g2[0, :], g2[1, :], 1, (node_num2, edge_num2))

            H2 = scipy.sparse.csr_matrix(h1[0, :], h1[1, :], 1, (node_num1, edge_num1))
            H1 = scipy.sparse.csr_matrix(h2[0, :], h2[1, :], 1, (node_num2, edge_num2))

            KQ = edge_attr2.dot(edge_attr1.T)
            ct = np.ones(node_num2, node_num1)

    if op.exists(op.join(path, file1+'_Nattr')) and op.exists(op.join(path, file2+'_Nattr')):
        node_attr_index1 = read_gph(op.join(path, file1+'_Nattr'))
        node_attr_index1[1, :] += 1
        node_attr_index2 = read_gph(op.join(path, file2+'_Nattr'))
        node_attr_index2[1, :] += 1
        if node_num1 > node_num2:
            KP = compute_kp(node_attr_index1, node_attr_index2, node_num1, node_num2)
        else:
            KP = compute_kp(node_attr_index2, node_attr_index1, node_num2, node_num1)

        gphs = dict()
        gphs[1] = {'g': G1, 'h': H1}
        gphs[2] = {'g': G2, 'h': H2}
        asg_t = []
        par = {'nalp': 100, 'nitma': 100, 'nhst': 10, 'ip': 0, 'deb': 0}
        res = fgm_d(KP, KQ, ct, gphs, asgt=asg_t, par=par)
        sim = res['obj']
        permu = res['x']

        if math.fabs(sim - edge_num1 - 0.1*node_num1) < 0.00001:
            ifsame = 1
        else:
            ifsame = 0
        return sim, permu, ifsame


def gph_eg2inc(eg):
    """obtain incidence matrix for directed edges"""
    m = len(eg)
    nods = set()
    [nods.update(item) for item in eg]
    nods = list(nods)
    n = len(nods)
    eg = [[nods.index(e[0]), nods.index(e[1])] for e in eg]

    g = np.zeros((n, m))
    h = np.zeros((n, m))

    for i in range(m):
        g[eg[i][0], i] = 1
        h[eg[i][1], i] = 1
    return g, h


def test_func():
    res = read_gph(path="E:\\graph_matching\\fgm\\src\\asg\\fgm\\data_process_Nattr\\node_5_edge_5", fi="00b55c6c71afb175e973ff72a851cb5456ce_objdata.gdl_H")
    print res
    print max(np.max(res[0, :]), np.max(res[1, :]))