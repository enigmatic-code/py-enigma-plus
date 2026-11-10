#! python
# -*- mode: Python; py-indent-offset: 2; -*-

from __future__ import print_function

# a class for solving Penney's game

from enigma import (Rational, Matrix, subsets, join)

__author__ = "Jim Randell <jim.randell@gmail.com>"
__version__ = "2024-04-01"

class Penney(object):

  def __init__(self, values='HT', probabilities=None):
    n = len(values)
    if probabilities is None:
      Q = Rational()
      probabilities = [Q(1, n)] * n
    self.values = values
    # look up the probabilities by the values
    self.probabilities = dict(zip(values, probabilities))
    # checks
    assert sum(self.probabilities.values()) == 1

  # the probability of sequence X occuring before sequence Y
  def P(self, X, Y):
    n = len(X)
    assert n == len(Y)

    # prefix sequences we are interested in
    ss = list(join(x) for x in subsets(self.values, max_size=n - 1, select="M"))
    ss.extend([X, Y])

    # generate the matrices to solve
    (mA, mB) = ([], [])
    for s in ss:
      eq = dict((k, 0) for k in ss)
      eq[s] = 1
      if s == X:
        k = 1
      elif s == Y:
        k = 0
      else:
        k = 0
        for x in self.values:
          sx = s + x
          if len(sx) == n and sx != X and sx != Y:
            sx = sx[1:]
          eq[sx] -= self.probabilities[x]
      mA.append(list(eq[x] for x in ss))
      mB.append(k)

    # solve the equations
    V = Matrix.linear(mA, mB)

    # we're interested in the probability for the empty sequence
    return V[ss.index('')]
