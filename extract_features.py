#!/usr/bin/env python
#Authors: Alessandro L. Asoni <ale.luca.asoni@gmail.com>
#         Kyle Wong           <wong.kyle5@gmail.com>
#         Guannan Ren         <gren3@jhu.edu>

""" Feature-extraction from DNA sequences in FASTA format

This script uses kmer counts to convert DNA sequences into feature vectors in
SVM-light format.

It reads two sequence files in FASTA format. For each file it must be
specified if the data belongs to the positive or negative set.

All kmers of a specified length 'k' are extracted from both files. 
If the option 'neighbors' is not set, the set of extracted kmers is 
compressed so that the hamming distance between any two remaining kmers is 
greater than a specified threshold 't'.
Each remaining kmer corresponds to a feature. For every read in the two files
the number of times each kmer appears with up to 't' mismatches correspends to 
the value of the feature associated with it. 

If the option neighbors is set then for each extracted kmer all string
neighbors within edit (or hamming) distance t are generated, and each one is
associated with the same feature number in a Map. The input DNA sequences are
then parsed and the feature counts are obtained with a lookup in this Map

If the option neighbors is not set then each read in the FASTFA files is
processed by looking at each kmer in the read and incrementing the count of all the
kmers in the feature set whose hamming distance from that kmer is smaller than 
the threshold 't'.

If 'edit' distance is chosen, the option neighbors defaults to true. This is for
performance reasons. Using the dynamic programming approach to calculate edit
distance is very slow.

When the option --skip is set, only every t-th kmer is read in from the FASTFA
files. This option deafaults to true if 'edit' distance is chosen, because any
kmer that is skipped this way is necessarily in the t-neighborhood of the kmers
that are included.
It is recommended to always use this option for long kmers (k > 15) even when
using hamming distance. This is for performance reasons.

The following is a usage example:

    using edit distance (with neighbors, with skipping):
        feature_extraction -k kmer_size -t threshold -d edit --neighbors --skip epositive.fa negative.fa 
    using hamming distance (without neighbors, without skipping):
        feature_extraction -k kmer_size -t threshold -d hamm positive.fa negative.fa

Using neighbors (and therefore edit distance) is only racommended for small k.
For long kmers (k > 10) using neighbors becomes intractable. 

IMPORTANT: This script makes use of the python library 'bitarray'
"""


import sys
import string
import os
import argparse
import numpy
import math
from bitarray import bitarray
from bitarray import bitdiff

parser = argparse.ArgumentParser(description="Create feature vectors from kmer counts")
parser.add_argument("-k", "--k", type=int, help="kmer length", default=20)
parser.add_argument("-t", "--threshold", type=int, help="kmer neighborhood threshold", default=3)
parser.add_argument("-d", "--distance", choices=["edit", "hamm"], help="distance scheme")
parser.add_argument("--neighbors", help="use neighbors to create feature map", action="store_true")
parser.add_argument("--skip", help="extract every t-th kmer", action="store_true")
parser.add_argument("pos", help="positive data set")
parser.add_argument("neg", help="negative data set")
args = parser.parse_args()
print args

if args.distance == "edit" and not args.neighbors:
    print "Warning. 'edit' distance only supported using neighbors. Option --neighbors is set here!"
    args.neighbors = True

if args.distance == "edit" and not args.skip:
    print "Using 'edit' distance. Setting option --skip"
    args.skip = True

def khammDist(x, y, k):
    assert len(x) == len(y)
    nmm = 0
    for i in xrange(0, len(x)):
        if x[i] != y[i]:
            nmm += 1
            if nmm > k:
                return False
    return True

def stringNeighbors(st, alph, edits=1, gaps=True):
    """ Given a string, an alphabet, and a maximum edit or Hamming
        distance, return all strings within that distance. """
    ret = []
    def __editNeighborsHelp(st, edits, ii):
        for i in xrange(ii, len(st)):
            if edits > 0:
                if gaps:
                    # Insertion just before position i
                    for a in alph:
                        newst = st[:i] + a + st[i:]
                        __editNeighborsHelp(newst, edits - 1, ii)
                    # Deletion of position i
                    newst = st[:i] + st[i+1:]
                    __editNeighborsHelp(newst, edits - 1, ii+1)
                # Mismatch at position i
                for a in alph:
                    if a != st[i]:
                        newst = st[:i] + a + st[i+1:]
                        __editNeighborsHelp(newst, edits - 1, ii+1)
        if gaps and edits > 0:
            # Insertion just after last position
            for a in alph: ret.append(st + a)
        ret.append(st)
    __editNeighborsHelp(st, edits, 0)
    return ret

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

def kmerInSeq(seq, k, t, skip):
    """ Read a sequence. Return list of kmers extracted """
    kmerList = []
    if skip: 
        for i in xrange(0, len(seq)-k+1, t+1):
            substr = seq[i:i+k]
            kmerList.append(substr)
    else:
        for i in xrange(0, len(seq)-k+1):
            substr = seq[i:i+k]
            kmerList.append(substr)

    return kmerList


def createFeatureVector(pos_neg, seq, features, k, t, edit, neighb):
    """ Create SVM-light format feature vector """
    feature_vector = {}
    if neighb:
            if edit:
                for i in xrange(0,len(seq)-k-t+1):
                    for j in xrange(-t,t+1):
                        feature = features[seq[i:i+k+j]]
                        if feature not in feature_vector:
                            feature_vector[feature] = 0
                        feature_vector[feature] += 1
            else:
                for i in xrange(0,len(seq)-k+1):
                    feature = features[seq[i:i+k]]
                    if feature not in feature_vector:
                        feature_vector[feature] = 0
                    feature_vector[feature] += 1

    else:
        for (kmer,feature) in features.items():
            for i in xrange(0,len(seq)-k+1):
                if khammDist(kmer, substr, t):
                    if feature not in feature_vector:
                        feature_vector[feature] = 0
                    feature_vector[feature] += 1
        
    feature_list = []
    for (feature, count) in feature_vector.items():
        feature_list.append((feature,count))
    feature_list.sort(key=lambda tup: tup[0])
    vector = ''.join(pos_neg)
    for (feature,count) in feature_list:
        vector += " " + repr(feature) + ":" + repr(count)

    return vector

posSeqList = readFasta(args.pos)
negSeqList = readFasta(args.neg)


kset = set() 
for seq in posSeqList:
    for kmer in kmerInSeq(seq, args.k, args.threshold, args.skip):
        kset.add( kmer )

for seq in negSeqList:
    for kmer in kmerInSeq(seq, args.k, args.threshold, args.skip):
        kset.add( kmer )

klist = list(kset)


print "Done extracting all substrings of length " + repr(args.k) + " appearing in all sequences ... "
print repr(len(klist)) + " distinct kmers extracted ... "
    
featureMap = {}
f = 1

#variables to keep track of progress
progress_count = 0
total = len(klist)

edit = args.distance == "edit"
#using neighbors
if args.neighbors:
    print "Using neighbors to create string->feature mapping ... "
    for k in xrange(0,len(klist)):
        if progress_count % math.ceil(total/1000.0) == 0:
            p = (float(progress_count)/total)*100
            sys.stdout.write("\r%.2f%% progress. " % p)
            sys.stdout.flush()
        progress_count += 1

        for n in stringNeighbors(klist[k],"ACGT", edits=args.threshold, gaps=edit):
            if n in featureMap: continue
            featureMap[n] = f
        f += 1

#using hamming distance to reduce kmer set
else:
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
    junk.setall(False)
    for i in xrange(0, len(encodedklist)):
        if progress_count % math.ceil(total/10000.0) == 0:
            p = (float(progress_count)/total)*100
            reduction = total - junk.count()
            sys.stdout.write("\r%.2f%% progress. " % p + "Feature set reduced to %d kmers ... " % reduction)
            sys.stdout.flush()
        progress_count += 1
    
        if junk[i]: continue
        for j in xrange(i+1, len(encodedklist)):
            if junk[j]: continue
            if bitdiff(encodedklist[i],encodedklist[j])/2 <= args.threshold:
                junk[j] = 1
    print
    print "Done reducing kmer set ... "

    for i in xrange(0,len(encodedklist)):
        if junk[i]: continue
        featureMap[''.join(encodedklist[i].decode(encoding))] = f 
        f += 1



print "Creating positive feature vectors in svm-light format ... "
output = "pos_svm_light_" + repr(args.k) + "_" + repr(args.threshold) + "_" + args.distance

print "Feature vectors will be saved in " + output + " ... "

f = open(output,'w')

progress_count = 0
total = len(posSeqList)
for seq in posSeqList:
    if progress_count % math.ceil(total/100.0) == 0:
        p = (float(progress_count)/total)*100
        sys.stdout.write("\r%.2f%% progress ..." % p)
        sys.stdout.flush()
    progress_count += 1

    feature_vector = createFeatureVector('1', seq, featureMap, args.k, args.threshold, edit, args.neighbors)
    f.write(feature_vector)
    f.write("\n")
print
f.close()

print "Creating negative feature vectors in svm-light format ... "
output = "neg_svm_light_" + repr(args.k) + "_" + repr(args.threshold) + "_" + args.distance

print "Feature vectors will be saved in " + output + " ... "
f = open(output,'w')

progress_count = 0
total = len(negSeqList)
for seq in negSeqList:
    if progress_count % math.ceil(total/100.0) == 0:
        p = (float(progress_count)/total)*100
        sys.stdout.write("\r%.2f%% progress ..." % p)
        sys.stdout.flush()
    progress_count += 1

    feature_vector = createFeatureVector('-1', seq, featureMap, args.k, args.threshold, edit, args.neighbors)
    f.write(feature_vector)
    f.write("\n")

print
print "Done ..."
f.close()
