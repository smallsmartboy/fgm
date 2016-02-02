import shutil
import os.path
from see_mat import *

root = "E:/graph_matching/fgm/src/asg/fgm/data_process_Nattr/"
dir_list = []
f = open(os.path.join(root, 'dirlog2.txt'))
n = 0
for line in f:
    line = line.rstrip('\n')
    shutil.copy(os.path.join(root, line, line+"_result"), os.path.join(root,'result'))
    #n += len(open(os.path.join(root, line, "file_list")).readlines())
f.close()
print n

