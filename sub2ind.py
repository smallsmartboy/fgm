import numpy as np

a = np.random.random([2,3,4])
# print "a\n",a
# b = a.T.reshape(a.size)
# print "b:\n",b
#
# print "b1:\n", b[1]
res = 1
for item in a.shape:
    res *= item

helper = np.arange(res)
helper = helper.reshape(a.shape)
print "a\n", a
print "helper", helper.shape
print helper