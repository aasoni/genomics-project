import re
import sys
import string
import numpy
import math

def fastaKmerParser(fns, k):
    """ Parse all (overlapping) k-mers of length k from all the FASTA
        filenames provided. """
    for fn in fns:
        with open(fn, 'r') as fh:
            kmer = []
            off = 0
            for ln in fh:
                if ln[0] == '>':
                    kmer = []
                    off = 0
                else:
                    for c in ln:
                        if c.isalpha():
                            if len(kmer) == k:
                                kmer.pop(0)
                            kmer.append(c.upper())
                            off += 1
                            if len(kmer) == k:
                                yield (''.join(kmer), off - k)

def fastaKmerParserIslands(fns, k, isles):
    """ Yield kmer, boolean pairs where boolean indicates whether the
        kmer lies entirely within an island or not """
    curIsland = 0
    for kmer, off in fastaKmerParser(fns, k):
        while curIsland < len(isles) and off >= isles[curIsland][1]:
            curIsland += 1
        if curIsland < len(isles) and off >= isles[curIsland][0]:
            if off + k <= isles[curIsland][1]:
                yield kmer, "island"
        else:
            yield kmer, "non-island"

def parseIslands(fn):
    """ Parse a file with island annotations, per the output from Hao's
        model-based approach: http://rafalab.jhsph.edu/CGI/ """
    isles = []
    with open(fn, 'r') as fh:
        for ln in fh:
            chr, st, en, _ = string.split(ln, '\t', 3)
            st, en = int(st), int(en)
            isles.append((st-1, en))
    return isles

def islandDinucs(fns, ifn, stopAfter=None):
    """ Compile dinucleotide count tables for dinucleotides inside or
        outside CpG islands """
    hist = {}
    n = 0
    isles = parseIslands(ifn)
    # fastaKmerParserIslands parses all dinucleotides from the FASTA
    # file; returns the k-mer along with a label stating whether kmer
    # came from inside or outside an island 
    for kmer, lab in fastaKmerParserIslands(fns, 2, isles):
        if lab not in hist: hist[lab] = {}
        hist[lab][kmer] = hist[lab].get(kmer, 0) + 1
        n += 1
        if stopAfter is not None and n >= stopAfter:
            break # stop early; helpful for testing
    return hist

def islandTransitionTables(fns, ifn, stopAfter=None):
    """ Given FASTA filenames and island annotation filename, obtain
        the island / non-island dinucleotide frequency tables.  Turn
        them into transition probabilities for the island & non-island
        markov chains.  Also, make log ratio table. """
    hist = islandDinucs(fns, ifn, stopAfter)
    iTab = numpy.zeros(shape=(4, 4), dtype=float)
    nTab = numpy.zeros(shape=(4, 4), dtype=float)
    for i in xrange(0, 4):
        srcc = "ACGT"[i]
        itot, ntot = 0, 0
        for j in xrange(0, 4):
            dstc = "ACGT"[j]
            itot += hist["island"].get(srcc + dstc, 0)
            ntot += hist["non-island"].get(srcc + dstc, 0)
        assert itot > 0
        assert ntot > 0
        for j in xrange(0, 4):
            dstc = "ACGT"[j]
            icnt = hist["island"].get(srcc + dstc, 0)
            ncnt = hist["non-island"].get(srcc + dstc, 0)
            iTab[i, j] = icnt / float(itot)
            nTab[i, j] = ncnt / float(ntot)
    return iTab, nTab

def classify(seq, lrTab):
    """ Classify seq using the already-computed table of dinucleotide
        log ratios """
    bits = 0
    nucmap = { 'A':0, 'C':1, 'G':2, 'T':3 }
    for dinuc in [ seq[i:i+2] for i in xrange(0, len(seq)-1) ]:
        i, j = nucmap[dinuc[0]], nucmap[dinuc[1]]
        bits += lrTab[i, j]
    return bits
