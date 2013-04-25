#!/usr/bin/env python
from ghmm import *
from math import *
import sys
import os
import argparse
import itertools

###############################################################
## Program compute the sum of square difference of 2 DNA seq
## the length of kmer used for comparison can be changed in program
## command line : ./prob_two.py __file1__ __file2__
###############################################################

def extract_kmer(s,ln,dict):
	tot_kmer_count = 0
	for i in xrange (0,len(s)-ln+1):
		dict[(s[i:i+ln])] += 1
		tot_kmer_count += 1
	return tot_kmer_count

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Create feature vectors from kmer counts")
	parser.add_argument("pos", help="positive data set")
	parser.add_argument("neg", help="negative data set")
	args = parser.parse_args()

	ln = 1
	bp = ['A','C','T','G']
	dna = ln*[bp]
	kmers = list(itertools.product(*dna))
	kmer_dict1 = {}
	kmer_dict2 = {}
	for i in xrange(len(kmers)):
		kmer_dict1[''.join(kmers[i])] = 0
		kmer_dict2[''.join(kmers[i])] = 0
	# print kmer_dict

	temp_seq = ""
	with open(args.pos, 'r') as fafh:
		input_str = fafh.readline()  #get rid of initial >
		while True:
			input_str = fafh.readline()
			if len(input_str) == 0: break
			input_str = input_str.strip().upper()

			if '>' in input_str:
				# temp_seq = ""
				trained = True
			else:
				temp_seq += input_str
	tot_kmer_count1 = extract_kmer(temp_seq,ln,kmer_dict1)
	tot_percent = 0
	for k in kmer_dict1:
		kmer_dict1[k] = float(kmer_dict1[k])/float(tot_kmer_count1)
		tot_percent += kmer_dict1[k]

	# print tot_kmer_count
	print kmer_dict1
	# print tot_percent

	temp_seq = ""
	with open(args.neg, 'r') as fafh:
		input_str = fafh.readline()  #get rid of initial >
		while True:
			input_str = fafh.readline()
			if len(input_str) == 0: break
			input_str = input_str.strip().upper()

			if '>' in input_str:
				# temp_seq = ""
				trained = True
			else:
				temp_seq += input_str
	tot_kmer_count2 = extract_kmer(temp_seq,ln,kmer_dict2)
	tot_percent = 0
	for k in kmer_dict2:
		kmer_dict2[k] = float(kmer_dict2[k])/float(tot_kmer_count2)
		tot_percent += kmer_dict2[k]

	print kmer_dict2

	sq_diff_sum = 0
	for k in kmer_dict1:
		sq_diff_sum += math.sqrt(math.pow((float(kmer_dict1[k])-float(kmer_dict2[k])),2))
	print sq_diff_sum
