#from fgm_file import read_eg_eattr
from fgm import *
from fgm_file import *
from scipy import io as sio
import numpy as np
import pdb





root = "E:\\data\\shellcode"
name1 = "d8479b994c3a4f84c57389c7a7415159a1de_sv.gdl"
name2 = "dd5464d5aacbed5bf5ec63dc7d4b33022a45_sv.gdl"
# name2 = "d8479b994c3a4f84c57389c7a7415159a1de_sv.gdl"
e1, ea1 = read_eg_eattr(root, name1)
e2, ea2 = read_eg_eattr(root, name2)
g1, h1 = gph_eg2inc(e1)
g2, h2 = gph_eg2inc(e2)
node_attr1 = read_node_nattr(root, name1+'_node')
node_attr2 = read_node_nattr(root, name2 + '_node')
node_num1 = np.shape(g1)[0]
node_num2 = np.shape(g2)[0]

ct = np.ones([node_num1, node_num2], dtype=np.int8)
kp = compute_kp(node_attr1, node_attr2, node_num1, node_num2)
kq = ea1.dot(ea2.T)

gphs = dict()
gphs[1] = {'g': g1, 'h': h1}
gphs[2] = {'g': g2, 'h': h2}
asgt = {}
par = {'nalp': 100, 'nitma': 100, 'nhst': 10, 'ip': 0, 'deb': 0}
res = fgm_d(kp, kq, ct, gphs, asgt, par)
print res['obj']
# print res['x']