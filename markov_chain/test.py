#!/usr/bin/env python
import re
import sys
import string
import numpy
import math
import argparse
import itertools

ln =1
bits = 0
bp = ['A','C','G','T']
dna = ln*[bp]
kmers = list(itertools.product(*dna))
# print kmers
for i in xrange(len(kmers)):
	kmers[i] = ''.join(kmers[i])
	print kmers[i]
count = 0
nucmap = {}
for k in kmers:
	nucmap[k] = count
	count += 1
# for n in sorted(nucmap.iterkeys()):
# 	print n + str(nucmap[n])
seq = "ACGTA"
for dinuc in [ seq[i:i+ln+1] for i in xrange(0, len(seq)-ln) ]:
	print dinuc
	i, j = nucmap[dinuc[0:ln]], nucmap[dinuc[1:ln+1]]
	print i
	print j
    # bits += lrTab[i, j]
# return bits