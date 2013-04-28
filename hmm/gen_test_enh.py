#!/usr/bin/env python 
import os
import sys
import random
import argparse

# denotation - 	ecoli 0
#				staph 1
#				difficile 2

indir = '../'
rand_seq = ""
marker = ""
enh_ind = []
null_ind = []

# parser = argparse.ArgumentParser(description="Create test files")
# parser.add_argument("test_input", help="test data set")
# parser.add_argument("test_check", help="actual result set")
# args = parser.parse_args()

for root, dirs, filenames in os.walk(indir):
  for f in filenames:
    if 'enh_fb' in f:
      mark = '0'
      # print f
      log_enh = open(os.path.join(root, f), 'r')
      for i in xrange(2453/3):
        enh_ind.append(random.randint(1,2453))  #getting random seq from batch
      # for j in enh_ind:
      #   print j
      enh_ind_count = 1
      input_str = log_enh.readline()  #get rid of initial >
      temp_seq = ""
      while True:
        input_str = log_enh.readline()
        if len(input_str) == 0: break
        input_str = input_str.strip().upper()
        if '>' in input_str:
          if enh_ind_count in enh_ind:
            print input_str
            print temp_seq
          temp_seq = ""
          enh_ind_count += 1
        else:
          temp_seq += input_str
      if enh_ind_count in enh_ind:
        print input_str
        print temp_seq
      log_enh.close()
    elif 'nullseqsi' in f:
      mark = '1'
      # print f   
      log_null = open(os.path.join(root, f), 'r')
      for i in xrange(4000/3):
        null_ind.append(random.randint(1,4000))  #getting random seq from batch
      null_ind_count = 1
      input_str = log_null.readline()  #get rid of initial >
      temp_seq = ""
      while True:
        input_str = log_null.readline()
        if len(input_str) == 0: break
        input_str = input_str.strip().upper()
        if '>' in input_str:
          if null_ind_count in null_ind:
            print input_str
            print temp_seq
          temp_seq = ""
          null_ind_count += 1
        else:
          temp_seq += input_str
      if null_ind_count in null_ind:
        print input_str
        print temp_seq
      log_null.close()