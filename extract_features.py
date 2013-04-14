#!/usr/bin/env python
#Authors: Alessandro L. Asoni <ale.luca.asoni@gmail.com>
#         Kyle Wong           <wong.kyle5@gmail.com>
#         Guannan Ren         <gren3@jhu.edu>

""" Feature-extraction from DNA sequences in FASTA format

This script uses kmer counts to convert DNA sequences into feature vectors in
SVM-light format.

It reads two sequence files in FASTA format. For each file it must be
specified if the data belongs to the positive or negative set.

All kmers of a specified length 'k' are extracted from both files. The set of
extracted kmers is compressed so that the edit (or hamming) distance between any
two remaining kmers is greater than a specified threshold 't'.

Each remaining kmer corresponds to a feature. For every read in the two files
the number of times each kmer appears with up to 't' edits correspends to the value of
the feature associated with it. 

The following is a usage example:

    using edit distance:
        feature_extraction -k kmer_size -t threshold -d edit epositive.fa negative.fa 
    using hamming distance:
        feature_extraction -k kmer_size -t threshold -d hamm positive.fa negative.fa

"""

import sys
import string
import os
import argparse
import numpy

parser = argparse.ArgumentParser(description="Create feature vectors from kmer counts")
parser.add_argument("-k", "--k", type=int, help="kmer length", default=20)
parser.add_argument("-t", "--threshold", type=int, help="kmer neighborhood threshold", default=3)
parser.add_argument("-d", "--distance", choices=["edit", "hamm"], help="distance scheme")
parser.add_argument("pos", help="positive data set")
parser.add_argument("neg", help="negative data set")
args = parser.parse_args()

def kedDistDp(x, y, k):
    """ Calculate edit distance between sequences x and y using
        matrix dynamic programming.  Return distance. """
    D = numpy.zeros((len(x)+1, len(y)+1), dtype=int)
    D[:] = k+1
    D[0, 1:] = range(1, len(y)+1)
    D[1:, 0] = range(1, len(x)+1)
    for i in xrange(1, len(x)+1):
        for j in xrange(max(1, i-k), min(i+k+1,len(y)+1)):
            delt = 1 if x[i-1] != y[j-1] else 0
            D[i, j] = min(D[i-1, j-1]+delt, D[i-1, j]+1, D[i, j-1]+1)
    return D[len(x), len(y)] <= k

def khammDist(x, y, k):
    """ Calculuate hamming distance between sequences x and y """
    assert len(x) == len(y)
    nmm = 0
    for i in xrange(0,	len(x)):
        if x[i]	!= y[i]:
            nmm += 1
            if nmm > k:
                return False
    return True

def partition(p, parts=2):
    """ Divide p into non-overlapping partitions. If there are excess
    characters, distribute them round-robin starting with 1st. """
    base, mod = len(p) / parts, len(p) % parts
    idx = 0
    ps = []

    modAdjust = 1
    for i in xrange(0, parts):
        if i >= mod:
            modAdjust = 0
        newIdx = idx + base + modAdjust
        ps.append((p[idx:newIdx], idx))
        idx = newIdx
    return ps

class kmerIndex(object):
    """ Substring index using hash table map. Maps substrings to kmers in which
    they occur """

    def __init__(self, ln=12, t=0, gaps=True):
        """ Create index with substrings of length 'ln' """
        self.ln = ln
        self.index = {}
        self.gaps = gaps

        # number of edits/mismatches 
        self.t = t  

    def query(self, p):
        """ Return candidate kmers for p """
        return self.index.get(p[:self.ln]) or []

    def insert(self, kmer):
        """ insert kmer in index if it doesn't appear already with up to t
        edits/mismatches """
        for part, off in partition(kmer, self.t+1):
            for hit in self.query(part):
                d = 0
                if self.gaps and kedDistDp(kmer, hit, self.t):
                    return
                elif khammDist(kmer, hit, self.t):
                    return
                
        l = self.ln
        for i in xrange(0, len(kmer)-l+1):
            substr = kmer[i:i+l]
            if substr not in self.index:
                self.index[substr] = []
            self.index[substr].append(kmer)

def readFasta(fn):
    """ Read a fasta file. Return list of sequences """
    reads = []
    if not os.path.exists(fn):
        raise RuntimeError("Fasta file doesn't exist: '%s'" % fn)
    with open(fn, 'r') as fafh:
        for ln in fafh:
            ln = ln.rstrip().upper()
            if ln.startswith('>'):
                continue
            else:
                seq = []
                seq.append(ln)
                reads.append( ''.join(seq) )

    return reads

def kmerInSeq(seq, k):
    """ Read a sequence. Return list of kmers """
    kmerList = []
    for i in xrange(0, len(seq)-k+1):
        substr = seq[i:i+k]
        kmerList.append(substr)

    return kmerList

posSeqList = readFasta(args.pos)
negSeqList = readFasta(args.neg)

klist = []
for seq in posSeqList:
    klist.extend( kmerInSeq(seq, args.k) )

for seq in negSeqList:
    klist.extend( kmerInSeq(seq, args.k) )

kset = list(set(klist))
kmerIndex = kmerIndex(args.k/(args.threshold+1),args.threshold,True)

count = 0
for k in kset:
    if count == 1000:
        break
    kmerIndex.insert(k)
    count = count+1

print len(kmerIndex.index)

