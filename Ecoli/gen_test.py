#!/usr/bin/env python 
import os
import sys
import random

# denotation - 	ecoli 0
#				staph 1
#				difficile 2

indir = '.'
rand_seq = ""
marker = ""
for root, dirs, filenames in os.walk(indir):
  for f in filenames:
    if 'fna' in f:
      log = open(os.path.join(root, f), 'r')
      input_str = log.readline()
      temp_seq = ""
      while True:
        input_str = log.readline()
        if len(input_str) == 0: break
        input_str = input_str.strip().upper()
        temp_seq += input_str

      for i in xrange(1):
        start = random.randint(1,10000)  #getting starting position
        length = 100 #random.randint(100,1000) #getting length of random sequence, ~read length of 100 to 1000
        rand_seq = temp_seq[start:start+length]
        print ">"
        print rand_seq
# print "seq: "
# print rand_seq
# print "marker: "
# print marker