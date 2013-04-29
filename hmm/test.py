import sys
import itertools
from ghmm import *

ln = 2
bp = ['A','C','G','T']
dna = ln*[bp]
kmers = list(itertools.product(*dna))
kmer_dict = {}
for i in xrange(len(kmers)):
	kmer_dict[''.join(kmers[i])] = 1
for j in sorted(kmer_dict.iterkeys()):
	print j