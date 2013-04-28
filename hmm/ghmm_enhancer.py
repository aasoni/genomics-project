#!/usr/bin/env python
from ghmm import *
import sys
import os
import argparse
import itertools

def extract_kmer(s,ln,dict):
	tot_kmer_count = 0
	for i in xrange (0,len(s)-ln+1):
		dict[(s[i:i+ln])] += 1
		tot_kmer_count += 1
	return tot_kmer_count

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Create feature vectors from kmer counts")
	parser.add_argument("pos", help="positive data set")
	# parser.add_argument("neg", help="negative data set")
	args = parser.parse_args()

	ln = 1
	bp = ['A','C','T','G']
	dna = ln*[bp]
	kmers = list(itertools.product(*dna))
	kmer_dict = {}
	for i in xrange(len(kmers)):
		kmer_dict[''.join(kmers[i])] = 0
	# print kmer_dict

	temp_seq = ""
	tot_kmer_count = 0
	with open(args.pos, 'r') as fafh:
		input_str = fafh.readline()  #get rid of initial >
		while True:
			input_str = fafh.readline()
			if len(input_str) == 0: break
			input_str = input_str.strip().upper()

			if '>' in input_str:
				tot_kmer_count += extract_kmer(temp_seq,ln,kmer_dict)
				temp_seq = ""
			else:
				temp_seq += input_str
		tot_kmer_count += extract_kmer(temp_seq,ln,kmer_dict)
	tot_percent = 0
	for k in kmer_dict:
		kmer_dict[k] = float(kmer_dict[k])/float(tot_kmer_count)
		tot_percent += kmer_dict[k]

	# print tot_kmer_count
	print kmer_dict
	# print tot_percent




	# dna = ['A','C','T','G']
	# sigma = Alphabet(dna)

	# A = [[0.9, 0.1], [0.3, 0.7]]
	# normal = [.25,.15,.35,.25]
	# enhancer = [.25,.25,.25,.25]
	# B=[normal,enhancer]
	# pi = [0.5] * 2
	# m=HMMFromMatrices(sigma,DiscreteDistribution(sigma),A,B,pi)
	# # print m

	# temp_seq = ""
	# trained = False

	# with open(args.pos, 'r') as fafh:
	# 	input_str = fafh.readline()  #get rid of initial >
	# 	# while trained == False:
	# 	while True:
	# 		input_str = fafh.readline()
	# 		if len(input_str) == 0: break
	# 		input_str = input_str.strip().upper()

	# 		if '>' in input_str:
	# 			temp_seq = ""
	# 			trained = True
	# 		else:
	# 			temp_seq += input_str
	# 	train_seq = EmissionSequence(sigma, list(temp_seq))
	# 	m.baumWelch(train_seq)

	# with open(args.neg, 'r') as fafh:
	# 	input_str = fafh.readline()  #get rid of initial >
	# 	for input_str in fafh:
	# 		input_str = input_str.strip().upper()
	# 		if '>' in input_str:
	# 			#running 20mer test seq
	# 			test_seq = EmissionSequence(sigma, list(temp_seq[:20]))
	# 			hidden_state,log = m.viterbi(test_seq)
	# 			print log
	# 			temp_seq = ""
	# 		else:
	# 			temp_seq += input_str
	# 	test_seq = EmissionSequence(sigma, list(temp_seq[:20]))
	# 	hidden_state,log =  m.viterbi(test_seq)
	# 	print log