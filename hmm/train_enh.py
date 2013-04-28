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

	ln = 5
	bp = ['A','C','T','G']
	dna = ln*[bp]
	kmers = list(itertools.product(*dna))
	kmer_dict = {}
	for i in xrange(len(kmers)):
		kmer_dict[''.join(kmers[i])] = 0
	# print kmer_dict

	temp_seq = ""
	tot_kmer_count = 0
	read_count = 1
	with open(args.pos, 'r') as fafh:
		input_str = fafh.readline()  #get rid of initial >
		while True:
			input_str = fafh.readline()
			if len(input_str) == 0: break
			input_str = input_str.strip().upper()

			if '>' in input_str:
				tot_kmer_count += extract_kmer(temp_seq,ln,kmer_dict)
				temp_seq = ""
				read_count += 1
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
	print read_count