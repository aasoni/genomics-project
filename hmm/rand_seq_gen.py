#!/usr/bin/env python
import os
import sys
import random

# denotation - 	ecoli 0
#				staph 1
#				difficile 2

indir = '../species_test/'
rand_seq = ""
marker = ""
for root, dirs, filenames in os.walk(indir):
    for f in filenames:
    	# print f
    	if 'ecoli' in f:
    		mark = '0'
    	elif 'staph' in f:
    		mark = '1'
    	elif 'difficile' in f:
    		mark = '2'
    	# print mark
        log = open(os.path.join(root, f), 'r')
        input_str = log.readline()
       	if '>' in input_str:
       		start = random.randint(1,2000)  #getting starting position
       		length = random.randint(100,1000) #getting length of random sequence, ~read length of 100 to 1000
       		temp_seq = ""
       		while True:
       			input_str = log.readline()
       			if len(input_str) == 0: break
       			input_str = input_str.strip().upper()
       			temp_seq += input_str
		rand_seq += temp_seq[start:start+length]
		marker += length*mark
print "seq: "
print rand_seq
print "marker: "
print marker