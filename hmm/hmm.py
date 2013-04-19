#!/usr/bin/env python
import numpy

class HMM(object):
    """ Simple Hidden Markov Model implementation.  User provides
        transition, emission and initial probabilities as dictionaries.
        Transition and emission probabilities are expressed as a dict
        mapping a 2-character code onto the floating-point probability
        for that table entry.  States and emissions are represented
        with single characters.  """

    def __init__(self, A, E, I):
        """ Initialize the HMM given transition, emission and initial
            probability tables. """
        self.Q = set()  #transitions
        self.S = set()  #emissions

        #adding all variable letters
        for a, prob in A.iteritems():
            assert len(a) == 2
            asrc, adst = a[0], a[1]
            self.Q.add(asrc)
            self.Q.add(adst)
        for e, prob in E.iteritems():
            assert len(e) == 2
            eq, es = e[0], e[1]
            self.Q.add(eq)
            self.S.add(es)
        self.Q = sorted(list(self.Q))
        self.S = sorted(list(self.S))

        #create dictionaries for Q and S
        qmap, smap = {}, {}
        for i in xrange(0, len(self.Q)):
            qmap[self.Q[i]] = i
        for i in xrange(0, len(self.S)):
            smap[self.S[i]] = i


        lenq = len(self.Q)
        self.A = numpy.zeros(shape=(lenq, lenq), dtype=float)
        self.E = numpy.zeros(shape=(lenq, len(self.S)), dtype=float)
        self.I = [ 0.0 ] * len(self.Q)


        #create probability matrices
        for a, prob in A.iteritems():
            asrc, adst = a[0], a[1]
            assert asrc in qmap
            assert adst in qmap
            self.A[qmap[asrc], qmap[adst]] = prob
        for e, prob in E.iteritems():
            eq, es = e[0], e[1]
            assert eq in qmap
            assert es in smap
            self.E[qmap[eq], smap[es]] = prob
        for a, prob in I.iteritems():
            self.I[qmap[a]] = prob
        self.qmap, self.smap = qmap, smap
        # Make log versions for log-space functions
        self.Alog = numpy.log2(self.A)
        self.Elog = numpy.log2(self.E)
        self.Ilog = numpy.log2(self.I)
    
    def jointProb(self, p, x):
        """ Return joint probability of path p and emission string x """
        p = map(self.qmap.get, p) # turn state characters into ids
        x = map(self.smap.get, x) # turn emission characters into ids
        tot = self.I[p[0]]
        for i in xrange(1, len(p)):
            tot *= self.A[p[i-1], p[i]]
        for i in xrange(0, len(p)):
            tot *= self.E[p[i], x[i]]
        return tot
    
    def viterbi(self, x):
        """ Given sequence of emissions, return the most probable path
            along with its probability. """
        x = map(self.smap.get, x) # turn emission characters into ids
        nrow, ncol = len(self.Q), len(x)
        mat   = numpy.zeros(shape=(nrow, ncol), dtype=float) # prob
        matTb = numpy.zeros(shape=(nrow, ncol), dtype=int)   # traceback
        # Fill in first column
        for i in xrange(0, nrow):
            mat[i, 0] = self.E[i, x[0]] * self.I[i]
            # print mat[i,0]   #0.5*0.5, 0.2*0.5
            # print self.E[i, x[0]]
        # Fill in rest of prob and Tb tables
        for j in xrange(1, ncol):
            for i in xrange(0, nrow):
                ep = self.E[i, x[j]]
                # print ep
                mx, mxi = mat[0, j-1] * self.A[0, i] * ep, 0
                for i2 in xrange(1, nrow):
                    pr = mat[i2, j-1] * self.A[i2, i] * ep
                    if pr > mx:
                        mx, mxi = pr, i2
                mat[i, j], matTb[i, j] = mx, mxi

        # Find final state with maximal probability
        omx, omxi = mat[0, ncol-1], 0
        for i in xrange(1, nrow):
            if mat[i, ncol-1] > omx:
                omx, omxi = mat[i, ncol-1], i
                
        # Traceback
        i, p = omxi, [omxi]
        for j in xrange(ncol-1, 0, -1):
            i = matTb[i, j]
            p.append(i)
        p = map(lambda x: self.Q[x], p[::-1])
        return omx, p # Return probability and path
    
    def viterbiL(self, x):
        """ Given sequence of emissions, return the most probable path
            along with its log base 2 probability.  Do all calculations in log
            space to avoid underflow. """
        x = map(self.smap.get, x) # turn emission characters into ids
        nrow, ncol = len(self.Q), len(x)
        mat   = numpy.zeros(shape=(nrow, ncol), dtype=float) # prob
        matTb = numpy.zeros(shape=(nrow, ncol), dtype=int)   # traceback
        # Fill in first column
        for i in xrange(0, nrow):
            mat[i, 0] = self.Elog[i, x[0]] + self.Ilog[i]
            # print mat[i,0]
        # Fill in rest of log prob and Tb tables
        for j in xrange(1, ncol):
            for i in xrange(0, nrow):
                ep = self.Elog[i, x[j]]
                mx, mxi = mat[0, j-1] + self.Alog[0, i] + ep, 0
                for i2 in xrange(1, nrow):
                    pr = mat[i2, j-1] + self.Alog[i2, i] + ep
                    if pr > mx:
                        mx, mxi = pr, i2
                mat[i, j], matTb[i, j] = mx, mxi
        # Find final state with maximal log probability
        omx, omxi = mat[0, ncol-1], 0
        for i in xrange(1, nrow):
            if mat[i, ncol-1] > omx:
                omx, omxi = mat[i, ncol-1], i
        # Traceback
        i, p = omxi, [omxi]
        for j in xrange(ncol-1, 0, -1):
            i = matTb[i, j]
            p.append(i)
        p = map(lambda x: self.Q[x], p[::-1])
        return omx, p # Return log probability and path


""" Example input for dishonest casino: """
hmm = HMM({"FF":0.6,"FL":0.4,"LF":0.4,"LL":0.6},{"FH":0.5,"FT":0.5,"LH":0.8,"LT":0.2},{"F":0.5,"L":0.5})
# prob,_ = hmm.viterbiL("THTHHHTHTTH")  #-18
# prob,_ = hmm.viterbiL("HHHHHHHHHHH")  #-11
# prob,_ = hmm.viterbiL("TTTTTTTTTTT")  #-19
prob,_ = hmm.viterbi("TTTTTTTTTTH")
# print prob
