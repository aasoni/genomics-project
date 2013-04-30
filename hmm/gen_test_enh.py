#!/usr/bin/env python 
import os
import sys
import random
import argparse

indir = '../'
rand_seq = ""
marker = ""
enh_ind = []
null_ind = []

num_enh = 1
num_null = 1
for root, dirs, filenames in os.walk(indir):
  for f in filenames:
    if 'enh_fb' in f:
      mark = '0'
      # print f
      log_enh = open(os.path.join(root, f), 'r')
      enh_ind_count = 1
      input_str = log_enh.readline()  #get rid of initial >
      temp_seq = ""
      while True:
        input_str = log_enh.readline()
        if len(input_str) == 0: break
        input_str = input_str.strip().upper()
        if '>' in input_str:
          if enh_ind_count > 2*2453/3:
            num_enh += 1
            print input_str
            print temp_seq
          temp_seq = ""
          enh_ind_count += 1
        else:
          temp_seq += input_str
      print input_str
      print temp_seq
      log_enh.close()
    elif 'nullseqsi' in f:
      mark = '1'
      # print f   
      log_null = open(os.path.join(root, f), 'r')
      null_ind_count = 1
      input_str = log_null.readline()  #get rid of initial >
      temp_seq = ""
      while True:
        input_str = log_null.readline()
        if len(input_str) == 0: break
        input_str = input_str.strip().upper()
        if '>' in input_str:
          if null_ind_count > 2*4000/3:
            num_null += 1
            print input_str
            print temp_seq
          temp_seq = ""
          null_ind_count += 1
        else:
          temp_seq += input_str
      print input_str
      print temp_seq
      log_null.close()
print "Number of enhancer sequences at the start: " + str(num_enh)
print "Number of non-enhancer sequences at the end: " + str(num_null)
    