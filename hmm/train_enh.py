#!/usr/bin/env python
from ghmm import *
import sys
import os
import argparse
import itertools

def extract_kmer(s,ln,transition_dict,tot_kmer_count):
	for i in xrange (0,len(s)-ln):  #took away +1
		if s[i:i+ln] not in transition_dict:
			bp = ['A','C','T','G']
			dna = ln*[bp]
			kmers = list(itertools.product(*dna))
			kmer_dict = {}
			for i in xrange(len(kmers)):
				kmer_dict[''.join(kmers[i])] = 0
			transition_dict[(s[i:i+ln])] = kmer_dict
		transition_arr = transition_dict[(s[i:i+ln])]
		if s[i+1:i+ln+1] not in transition_arr:
			transition_arr[s[i+1:i+ln+1]] = 0			
		transition_arr[s[i+1:i+ln+1]] += 1
		transition_dict[(s[i:i+ln])] = transition_arr
		if s[i:i+ln] not in tot_kmer_count:
			tot_kmer_count[s[i:i+ln]] = 0			
		tot_kmer_count[s[i:i+ln]] += 1

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Create feature vectors from kmer counts")
	parser.add_argument("seq", help="input sequence data set")
	parser.add_argument("which", help="input which (enh or null) data set")
	args = parser.parse_args()

	ln = 2
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
	kmer_dict = {}
	tot_kmer_count = {}
	with open(args.seq, 'r') as fafh:
		input_str = fafh.readline()  #get rid of initial >
		while True:
			input_str = fafh.readline()
			if len(input_str) == 0: break
			input_str = input_str.strip().upper()

			if '>' in input_str:
				extract_kmer(temp_seq,ln,kmer_dict,tot_kmer_count)
				temp_seq = ""
				read_count += 1
			else:
				temp_seq += input_str
		extract_kmer(temp_seq,ln,kmer_dict,tot_kmer_count)
	# tot_percent = 0
	for k in sorted(kmer_dict.iterkeys()):
	# for k in kmer_dict:
		k_dict_dict = kmer_dict[k]
		for j in k_dict_dict:
			k_dict_dict[j] = float(k_dict_dict[j])/float(tot_kmer_count[k])
		if args.which == "enhancer":
			dummy_arr = []
			for i in sorted(k_dict_dict.iterkeys()):
				dummy_arr.append(k_dict_dict[i])
			# dummy_arr = k_dict_dict.values()
			dummy_arr.extend([0.0]*(4**ln))
			print k + str(dummy_arr)
		elif args.which == "null":
			dummy_arr = [0.0]*(4**ln)
			for i in sorted(k_dict_dict.iterkeys()):
				dummy_arr.append(k_dict_dict[i])
			# dummy_arr.extend(k_dict_dict.values())
			print k + str(dummy_arr)

	print "-------------------transition------------------------"
	bp = ['A','C','G','T']
	dna = ln*[bp]
	kmers = list(itertools.product(*dna))
	kmer_dict = {}
	for i in xrange(len(kmers)):
		kmer_dict[''.join(kmers[i])] = [0]*4

	for j in sorted(kmer_dict.iterkeys()):
		temp_arr = kmer_dict[j]
		if j[-1:] == 'A':
			temp_arr[0] = 1;
		elif j[-1:] == 'C':
			temp_arr[1] = 1;
		elif j[-1:] == 'T':
			temp_arr[3] = 1;
		elif j[-1:] == 'G':
			temp_arr[2] = 1;
		kmer_dict[j] = temp_arr
		print j + str(kmer_dict[j])