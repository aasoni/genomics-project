#!/usr/bin/env python
#Authors: Alessandro L. Asoni <ale.luca.asoni@gmail.com>
#         Kyle Wong           <wong.kyle5@gmail.com>
#         Guannan Ren         <gren3@jhu.edu>

import sys

k = sys.argv[1]
t = sys.argv[2]
d = sys.argv[3]

in_pos = "pos_svm_light_" + str(k) + "_" + str(t) + "_" + str(d)
in_neg = "neg_svm_light_" + str(k) + "_" + str(t) + "_" + str(d)
out_train = str(k) + "_" + str(t) +  "_" + str(d) + ".train"
out_test = str(k) + "_" + str(t) +  "_" + str(d)+ ".test"

inp = open(in_pos, 'r') 
inn = open(in_neg, 'r')
tr = open(out_train, 'w')
ts = open(out_test,  'w') 

for i in xrange(0,1750):
    ln = inn.readline()
    tr.write(ln)
    ln = inp.readline()
    tr.write(ln)
    ln = inn.readline()
    tr.write(ln)

tr.close()

for i in xrange(0,500):
    ln = inn.readline()
    ts.write(ln)
    ln = inp.readline()
    ts.write(ln)

while True:
    ln = inp.readline()
    if len(ln) == 0: break
    ts.write(ln)

inn.close()
inp.close()
ts.close()
