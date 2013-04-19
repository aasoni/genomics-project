#!/usr/bin/env python
from ghmm import *
import sys
import os
import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Create feature vectors from kmer counts")
	parser.add_argument("pos", help="positive data set")
	parser.add_argument("neg", help="negative data set")
	args = parser.parse_args()

	dna = ['A','C','T','G']
	sigma = Alphabet(dna)

	A = [[0.9, 0.1], [0.3, 0.7]]
	normal = [.25,.15,.35,.25]
	enhancer = [.25,.25,.25,.25]
	B=[normal,enhancer]
	pi = [0.5] * 2
	m=HMMFromMatrices(sigma,DiscreteDistribution(sigma),A,B,pi)
	# print m

	temp_seq = ""
	trained = False

	with open(args.pos, 'r') as fafh:
		input_str = fafh.readline()  #get rid of initial >
		# while trained == False:
		while True:
			input_str = fafh.readline()
			if len(input_str) == 0: break
			input_str = input_str.strip().upper()

			if '>' in input_str:
				temp_seq = ""
				trained = True
			else:
				temp_seq += input_str
		train_seq = EmissionSequence(sigma, list(temp_seq))
		m.baumWelch(train_seq)

	with open(args.neg, 'r') as fafh:
		input_str = fafh.readline()  #get rid of initial >
		for input_str in fafh:
			input_str = input_str.strip().upper()
			if '>' in input_str:
				#running 20mer test seq
				test_seq = EmissionSequence(sigma, list(temp_seq[:20]))
				hidden_state,log = m.viterbi(test_seq)
				print log
				temp_seq = ""
			else:
				temp_seq += input_str
		test_seq = EmissionSequence(sigma, list(temp_seq[:20]))
		hidden_state,log =  m.viterbi(test_seq)
		print log