#!/usr/bin/env python
from ghmm import *
from UnfairCasino import *

sigma = IntegerRange(1,7)

A = [[0.5, 0.5], [0.5, 0.5]]
efair = [1.0 / 6] * 6
#print efair
eloaded = [3.0 / 13, 3.0 / 13, 2.0 / 13, 2.0 / 13, 2.0 / 13, 1.0 / 13]
B = [efair, eloaded]
pi = [0.5] * 2
m = HMMFromMatrices(sigma, DiscreteDistribution(sigma), A, B, pi)
#print m

obs_seq = m.sampleSingle(20)
#print obs_seq
obs = map(sigma.external, obs_seq)
#print obs
#[1, 5, 4, 1, 3, 4, 1, 1, 3, 4, 3, 5, 5, 1, 5, 2, 1, 5, 3, 5]

# print train_seq
m.baumWelch(train_seq)
# print train_seq
print m
v,log = m.viterbi(test_seq)
# print v
# my_seq = EmissionSequence(sigma, [1] * 20 + [6] * 10 + [1] * 40)
# print my_seq
# print m.viterbi(my_seq)

