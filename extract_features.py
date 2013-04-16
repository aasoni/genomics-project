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
from bitarray import bitarray
from bitarray import bitdiff

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
    D[0,0] = 0
    D[0, 1:] = range(1, len(y)+1)
    D[1:, 0] = range(1, len(x)+1)
    for i in xrange(1, len(x)+1):
        for j in xrange(max(1, i-k), min(i+k+1,len(y)+1)):
            delt = 1 if x[i-1] != y[j-1] else 0
            D[i, j] = min(D[i-1, j-1]+delt, D[i-1, j]+1, D[i, j-1]+1)
    return D[len(x), len(y)] <= k

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

klist = list(set(klist))

print "Done extracting all substrings of length " + repr(args.k) + " appearing in all sequences ... "
print repr(len(klist)) + " distinct kmers extracted ... "

encoding = {'A':bitarray('0001'), 'G':bitarray('0010'), 'C':bitarray('0100'), 'T':bitarray('1000')}
encodedklist = []
for kmer in klist:
    a = bitarray()
    a.encode(encoding,kmer)
    encodedklist.append(a)

print "Done encoding all kmers as bit arrays ... "
print "Starting compression of kmer set using hamming distance " + repr(args.threshold) + " ... "

#variables to keep track of progress
progress_count = 0
total = len(encodedklist)

#bitmask to keep track of kmers that can be eliminated
junk = bitarray(total)

for i in xrange(0, len(encodedklist)):
    if progress_count % 5000 == 0:
        p = progress_count/float(total)
        sys.stdout.write("\r%.2f%%" %p)
        sys.stdout.flush()
    progress_count += 1
    print progress_count

    if junk[i]: continue
    for j in xrange(0, len(encodedklist)):
        if bitdiff(encodedklist[i],encodedklist[j])/2 <= args.threshold:
            junk[j] = 1

print
print junk.count()
