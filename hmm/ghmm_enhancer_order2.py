#!/usr/bin/env python
from ghmm import *
import sys
import os
import argparse
import itertools
import random

def extract_kmer(s,ln):
	index = []
	for i in xrange (0,len(s)-ln+1):
		index.append(s[i:i+ln])
	return index

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Create feature vectors from kmer counts")
	parser.add_argument("pos", help="positive data set")
	parser.add_argument("neg", help="negative data set")
	# parser.add_argument("order", help="order of HMM")
	args = parser.parse_args()

	bp = ['A','C','T','G']

	dna = [bp,bp,bp]
	kmers = list(itertools.product(*dna))
	for i in xrange(len(kmers)):
		kmers[i] = ''.join(kmers[i])
	# print kmers
	sigma = Alphabet(kmers)

	normal = []
	enhancer = []
	for i in xrange(len(kmers)):
		normal.append(0.5)
		enhancer.append(0.5)
	E = [normal,enhancer]
	# print len(E)
	# print len(E[0])

	T = [[0.5, 0.5], [0.5, 0.5]]
	pi = [0.5] * 2
	m=HMMFromMatrices(sigma,DiscreteDistribution(sigma),T,E,pi)
	# print m

	temp_seq = ""
	# trained = False
	kmer_arr = []
	with open(args.pos, 'r') as fafh:
		input_str = fafh.readline()  #get rid of initial >
		# while trained == False:
		while True:
			input_str = fafh.readline()
			if len(input_str) == 0: break
			input_str = input_str.strip().upper()

			if '>' in input_str:
				# temp_seq = ""
				trained = True
			else:
				temp_seq += input_str
		train_arr = extract_kmer(temp_seq,3)

		# for i in kmer_arr:
			# print i
		# print train_arr
		train_seq = EmissionSequence(sigma, train_arr)
		m.baumWelch(train_seq)

		# print m

	test_arr = []
	with open(args.neg, 'r') as fafh:
		input_str = fafh.readline()  #get rid of initial >
		for input_str in fafh:
			input_str = input_str.strip().upper()
			if '>' in input_str:
				test_arr = extract_kmer(temp_seq[:20],3)
				test_seq = EmissionSequence(sigma, test_arr)
				hidden_state,log = m.viterbi(test_seq)
				print log
				temp_seq = ""
			else:
				temp_seq += input_str
		test_arr = extract_kmer(temp_seq[:20],3)
		test_seq = EmissionSequence(sigma, test_arr)
		hidden_state,log =  m.viterbi(test_seq)
		print log

