#!/usr/bin/env python
from ghmm import *

dna = ['a','c','t','g']
sigma = Alphabet(dna)

A = [[0.9, 0.1], [0.3, 0.7]]
normal = [.25,.15,.35,.25]
island = [.25,.25,.25,.25]
B=[normal,island]
pi = [0.5] * 2
m=HMMFromMatrices(sigma,DiscreteDistribution(sigma),A,B,pi)
# print m

# obs_seq = m.sampleSingle(20)
# print obs_seq
# obs = map(sigma.external, obs_seq)
# print obs

train_seq = EmissionSequence(sigma, ['c'] * 20 + ['t'] * 10 + ['c'] * 40)
m.baumWelch(train_seq)
# v = m.viterbi(train_seq)
# print v
obs_seq = m.sampleSingle(70)
print obs_seq
v = m.viterbi(obs_seq)
print v
