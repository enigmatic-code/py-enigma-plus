#!/usr/bin/env python3 -t
# -*- mode: Python; py-indent-offset: 2; -*-

from __future__ import print_function

# 

from enigma import enigma, basestring, chunk, join, printf

# shapes have 8 possible orientations

# R0 = rotated 0 degrees
# R1 = rotated 90 degrees
# R2 = rotated 180 degrees
# R3 = rotated 270 degrees
# M0 = mirrored and rotated 0 degrees
# M1 = mirrored and rotated 90 degrees
# M2 = mirrored and rotated 180 degrees
# M3 = mirrored and rotated 270 degrees
(R0, R1, R2, R3, M0, M1, M2, M3) = (1, 2, 4, 8, 16, 32, 64, 128)

# pre-computed polyominoes:
#   1-ominoes = O1
#   2-ominoes = I2
#   3-ominoes = I3 V3
#   4-ominoes = O4 I4 T4 Z4/Z4' L4/L4'
#   5-ominoes = I5 U5 T5 V5 X5 W5 L5/L5' Y5/Y5' P5/P5' N5/N5' F5/F5' Z5/Z5'

polyominoes = {

  # 1 square (1-omino, monomino) [1 achiral]
  #
  #  @ = O1

  # O1: 1 square - 1x1 block (achiral)
  "O1": (([(0, 0)], R0 | R1 | R2 | R3 | M0 | M1 | M2 | M3),),


  # 2 squares (2-ominos, dominos) [1 achiral]
  #
  #  @@ = I2

  # I2: 2 squares - 2x1 block (achiral)
  "I2": (
    ([(0, 0), (1, 0)], R0 | R2 | M0 | M2),
    ([(0, 0), (0, 1)], R1 | R3 | M1 | M3),    
  ),


  # 3 squares (3-ominos, trominos) [2 achiral]
  #
  #  @@@ = I3
  #
  #  @
  #  @@  = V3, L3

  # I3: 3 squares - linear, 3x1 block (achiral)
  "I3": (
    ([(0, 0), (1, 0), (2, 0)], R0 | R2 | M0 | M2),
    ([(0, 0), (0, 1), (0, 2)], R1 | R3 | M1 | M3),
  ),

  # V3: 3 squares - L shape (achiral)
  "V3": (
    ([(0, 0), (0, 1), (1, 0)], R0 | M1),
    ([(0, 0), (0, 1), (1, 1)], R1 | M2),
    ([(0, 1), (1, 0), (1, 1)], R2 | M3),
    ([(0, 0), (1, 0), (1, 1)], R3 | M0),
  ),


  # 4 squares (4-ominos, tetrominos) [3 achiral, 2 chiral]
  #
  #  @@
  #  @@   = O4
  #
  #  @@@@ = I4
  #
  #   @
  #  @@@  = T4
  #
  #  @@         |   @@
  #   @@  = Z4  |  @@   = Z4', S4
  #
  #    @        |  @
  #  @@@  = L4  |  @@@  = L4', R4

  # O4: 4 squares - square, 2x2 block (achiral)
  "O4": (
    ([(0, 0), (0, 1), (1, 0), (1, 1)], R0 | R1 | R2 | R3 | M0 | M1 | M2 | M3),
  ),

  # I4: 4 squares - 4x1 block (achiral)
  "I4": (
    ([(0, 0), (1, 0), (2, 0), (3, 0)], R0 | R2 | M0 | M2),
    ([(0, 0), (0, 1), (0, 2), (0, 3)], R1 | R3 | M1 | M3),
  ),

  # T4: 4 squares - T shape (achiral)
  "T4": (
    ([(0, 0), (1, 0), (1, 1), (2, 0)], R0 | M0),
    ([(0, 0), (0, 1), (0, 2), (1, 1)], R1 | M1),
    ([(0, 1), (1, 0), (1, 1), (2, 1)], R2 | M2),
    ([(0, 1), (1, 0), (1, 1), (1, 2)], R3 | M3),
  ),

  # Z4: 4 squares - Z shape (chiral)
  "Z4": (
    ([(0, 1), (1, 0), (1, 1), (2, 0)], R0 | R2),
    ([(0, 0), (0, 1), (1, 1), (1, 2)], R1 | R3),
    ([(0, 0), (1, 0), (1, 1), (2, 1)], M0 | M2),
    ([(0, 1), (0, 2), (1, 0), (1, 1)], M1 | M3),
  ),  

  # Z4': 4 squares - S shape (chiral)
  "Z4'": (
    ([(0, 0), (1, 0), (1, 1), (2, 1)], R0 | R2),
    ([(0, 1), (0, 2), (1, 0), (1, 1)], R1 | R3),
    ([(0, 1), (1, 0), (1, 1), (2, 0)], M0 | M2),
    ([(0, 0), (0, 1), (1, 1), (1, 2)], M1 | M3),
  ),  

  # L4: 4 squares - L shape (chiral)
  "L4": (
    ([(0, 0), (0, 1), (0, 2), (1, 0)], R0),
    ([(0, 0), (0, 1), (1, 1), (2, 1)], R1),
    ([(0, 2), (1, 0), (1, 1), (1, 2)], R2),
    ([(0, 0), (1, 0), (2, 0), (2, 1)], R3),
    ([(0, 0), (1, 0), (1, 1), (1, 2)], M0),
    ([(0, 0), (0, 1), (1, 0), (2, 0)], M1),
    ([(0, 0), (0, 1), (0, 2), (1, 2)], M2),
    ([(0, 1), (1, 1), (2, 0), (2, 1)], M3),
  ),
 
  # 4 squares - r shape (chiral)
  "L4'": (
    ([(0, 0), (1, 0), (1, 1), (1, 2)], R0),
    ([(0, 0), (0, 1), (1, 0), (2, 0)], R1),
    ([(0, 0), (0, 1), (0, 2), (1, 2)], R2),
    ([(0, 1), (1, 1), (2, 0), (2, 1)], R3),
    ([(0, 0), (0, 1), (0, 2), (1, 0)], M0),
    ([(0, 0), (0, 1), (1, 1), (2, 1)], M1),
    ([(0, 2), (1, 0), (1, 1), (1, 2)], M2),
    ([(0, 0), (1, 0), (2, 0), (2, 1)], M3),
  ),

  # 5 squares (5-ominoes, pentominoes) [6 achiral, 6 chiral]
  #
  #  @@@@@ = I5
  #
  #  @ @
  #  @@@   = U5
  #
  #  @@@
  #   @
  #   @    = T5
  #
  #  @
  #  @
  #  @@@   = V5
  #
  #   @
  #  @@@
  #   @    = X5
  #
  #    @
  #   @@
  #  @@    = W5
  #
  #     @        |  @
  #  @@@@  = L5  |  @@@@  = L5', R5
  #
  #    @         |   @
  #  @@@@  = Y5  |  @@@@  = Y5'
  #
  #  @@          |   @@
  #  @@@   = P5  |  @@@   = P5', Q5
  #
  #  @@          |    @@
  #   @@@  = N5  |  @@@   = N5'
  #
  #   @@         |  @@
  #  @@          |   @@
  #   @    = F5  |   @    = F5'
  #
  #  @@          |   @@
  #   @          |   @
  #   @@   = Z5  |  @@    = Z5'
  #

  # 5 squares - 5x1 block (achiral)
  "I5": (
    ([(0, 0), (1, 0), (2, 0), (3, 0), (4, 0)], R0 | R2 | M0 | M2),
    ([(0, 0), (0, 1), (0, 2), (0, 3), (0, 4)], R1 | R3 | M1 | M3),
  ),
  
  # 5 squares - U shape (achiral)
  "U5": (
    ([(0, 0), (0, 1), (1, 0), (2, 0), (2, 1)], R0 | M0),
    ([(0, 0), (0, 1), (0, 2), (1, 0), (1, 2)], R1 | M1),
    ([(0, 0), (0, 1), (1, 1), (2, 0), (2, 1)], R2 | M2),
    ([(0, 0), (0, 2), (1, 0), (1, 1), (1, 2)], R3 | M3),
  ),

  # 5 squares - T shape (achiral)
  "T5": (
    ([(0, 2), (1, 0), (1, 1), (1, 2), (2, 2)], R0 | M0),
    ([(0, 1), (1, 1), (2, 0), (2, 1), (2, 2)], R1 | M1),
    ([(0, 0), (1, 0), (1, 1), (1, 2), (2, 0)], R2 | M2),
    ([(0, 0), (0, 1), (0, 2), (1, 1), (2, 1)], R3 | M3),
  ),

  # 5 squares - V shape (achiral)
  "V5": (
    ([(0, 0), (0, 1), (0, 2), (1, 0), (2, 0)], R0 | M1),
    ([(0, 0), (0, 1), (0, 2), (1, 2), (2, 2)], R1 | M2),
    ([(0, 2), (1, 2), (2, 0), (2, 1), (2, 2)], R2 | M3),
    ([(0, 0), (1, 0), (2, 0), (2, 1), (2, 2)], R3 | M0),
  ),

  # 5 squares - X shape (achiral)
  "X5": (
    ([(0, 1), (1, 0), (1, 1), (1, 2), (2, 1)], R0 | R1 | R2 | R3 | M0 | M1 | M2 | M3),
  ),

  # 5 squares - W shape (achiral)
  "W5": (
    ([(0, 0), (1, 0), (1, 1), (2, 1), (2, 2)], R0 | M3),
    ([(0, 1), (0, 2), (1, 0), (1, 1), (2, 0)], R1 | M0),
    ([(0, 0), (0, 1), (1, 1), (1, 2), (2, 2)], R2 | M1),
    ([(0, 2), (1, 1), (1, 2), (2, 0), (2, 1)], R3 | M2),
  ),

  # 5 squares - L shape (chiral)
  "L5": (
    ([(0, 0), (1, 0), (2, 0), (3, 0), (3, 1)], R0),
    ([(0, 0), (0, 1), (0, 2), (0, 3), (1, 0)], R1),
    ([(0, 0), (0, 1), (1, 1), (2, 1), (3, 1)], R2),
    ([(0, 3), (1, 0), (1, 1), (1, 2), (1, 3)], R3),
    ([(0, 0), (0, 1), (1, 0), (2, 0), (3, 0)], M0),
    ([(0, 0), (0, 1), (0, 2), (0, 3), (1, 3)], M1),
    ([(0, 1), (1, 1), (2, 1), (3, 0), (3, 1)], M2),
    ([(0, 0), (1, 0), (1, 1), (1, 2), (1, 3)], M3),
  ),

  # 5 squares - r shape (chiral)
  "L5'": (
    ([(0, 0), (0, 1), (1, 0), (2, 0), (3, 0)], R0),
    ([(0, 0), (0, 1), (0, 2), (0, 3), (1, 3)], R1),
    ([(0, 1), (1, 1), (2, 1), (3, 0), (3, 1)], R2),
    ([(0, 0), (1, 0), (1, 1), (1, 2), (1, 3)], R3),
    ([(0, 0), (1, 0), (2, 0), (3, 0), (3, 1)], M0),
    ([(0, 0), (0, 1), (0, 2), (0, 3), (1, 0)], M1),
    ([(0, 0), (0, 1), (1, 1), (2, 1), (3, 1)], M2),
    ([(0, 3), (1, 0), (1, 1), (1, 2), (1, 3)], M3),
  ),

  # 5 squares - y shape (chiral)
  "Y5": (
    ([(0, 0), (1, 0), (2, 0), (2, 1), (3, 0)], R0),
    ([(0, 0), (0, 1), (0, 2), (0, 3), (1, 1)], R1),
    ([(0, 1), (1, 0), (1, 1), (2, 1), (3, 1)], R2),
    ([(0, 2), (1, 0), (1, 1), (1, 2), (1, 3)], R3),
    ([(0, 0), (1, 0), (1, 1), (2, 0), (3, 0)], M0),
    ([(0, 0), (0, 1), (0, 2), (0, 3), (1, 2)], M1),
    ([(0, 1), (1, 1), (2, 0), (2, 1), (3, 1)], M2),
    ([(0, 1), (1, 0), (1, 1), (1, 2), (1, 3)], M3),
  ),

  # 5 squares - y' shape (chiral)
  "Y5'": (
    ([(0, 0), (1, 0), (1, 1), (2, 0), (3, 0)], R0),
    ([(0, 0), (0, 1), (0, 2), (0, 3), (1, 2)], R1),
    ([(0, 1), (1, 1), (2, 0), (2, 1), (3, 1)], R2),
    ([(0, 1), (1, 0), (1, 1), (1, 2), (1, 3)], R3),
    ([(0, 0), (1, 0), (2, 0), (2, 1), (3, 0)], M0),
    ([(0, 0), (0, 1), (0, 2), (0, 3), (1, 1)], M1),
    ([(0, 1), (1, 0), (1, 1), (2, 1), (3, 1)], M2),
    ([(0, 2), (1, 0), (1, 1), (1, 2), (1, 3)], M3),
  ),

  # 5 squares - P shape (chiral)
  "P5": (
    ([(0, 0), (0, 1), (1, 0), (1, 1), (2, 0)], R0),
    ([(0, 0), (0, 1), (0, 2), (1, 1), (1, 2)], R1),
    ([(0, 1), (1, 0), (1, 1), (2, 0), (2, 1)], R2),
    ([(0, 0), (0, 1), (1, 0), (1, 1), (1, 2)], R3),
    ([(0, 0), (1, 0), (1, 1), (2, 0), (2, 1)], M0),
    ([(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)], M1),
    ([(0, 0), (0, 1), (1, 0), (1, 1), (2, 1)], M2),
    ([(0, 1), (0, 2), (1, 0), (1, 1), (1, 2)], M3),
  ),

  # 5 squares - q shape (chiral)
  "P5'": (
    ([(0, 0), (1, 0), (1, 1), (2, 0), (2, 1)], R0),
    ([(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)], R1),
    ([(0, 0), (0, 1), (1, 0), (1, 1), (2, 1)], R2),
    ([(0, 1), (0, 2), (1, 0), (1, 1), (1, 2)], R3),
    ([(0, 0), (0, 1), (1, 0), (1, 1), (2, 0)], M0),
    ([(0, 0), (0, 1), (0, 2), (1, 1), (1, 2)], M1),
    ([(0, 1), (1, 0), (1, 1), (2, 0), (2, 1)], M2),
    ([(0, 0), (0, 1), (1, 0), (1, 1), (1, 2)], M3),
  ),

  # 5 squares - N shape (chiral)
  "N5": (
    ([(0, 1), (1, 0), (1, 1), (2, 0), (3, 0)], R0),
    ([(0, 0), (0, 1), (0, 2), (1, 2), (1, 3)], R1),
    ([(0, 1), (1, 1), (2, 0), (2, 1), (3, 0)], R2),
    ([(0, 0), (0, 1), (1, 1), (1, 2), (1, 3)], R3),
    ([(0, 0), (1, 0), (2, 0), (2, 1), (3, 1)], M0),
    ([(0, 1), (0, 2), (0, 3), (1, 0), (1, 1)], M1),
    ([(0, 0), (1, 0), (1, 1), (2, 1), (3, 1)], M2),
    ([(0, 2), (0, 3), (1, 0), (1, 1), (1, 2)], M3),
  ),

  # 5 squares - N' shape (chiral)
  "N5'": (
    ([(0, 0), (1, 0), (2, 0), (2, 1), (3, 1)], R0),
    ([(0, 1), (0, 2), (0, 3), (1, 0), (1, 1)], R1),
    ([(0, 0), (1, 0), (1, 1), (2, 1), (3, 1)], R2),
    ([(0, 2), (0, 3), (1, 0), (1, 1), (1, 2)], R3),
    ([(0, 1), (1, 0), (1, 1), (2, 0), (3, 0)], M0),
    ([(0, 0), (0, 1), (0, 2), (1, 2), (1, 3)], M1),
    ([(0, 1), (1, 1), (2, 0), (2, 1), (3, 0)], M2),
    ([(0, 0), (0, 1), (1, 1), (1, 2), (1, 3)], M3),
  ),

  # 5 squares - F shape (chiral)
  "F5": (
    ([(0, 1), (1, 0), (1, 1), (1, 2), (2, 2)], R0),
    ([(0, 1), (1, 1), (1, 2), (2, 0), (2, 1)], R1),
    ([(0, 0), (1, 0), (1, 1), (1, 2), (2, 1)], R2),
    ([(0, 1), (0, 2), (1, 0), (1, 1), (2, 1)], R3),
    ([(0, 2), (1, 0), (1, 1), (1, 2), (2, 1)], M0),
    ([(0, 1), (1, 0), (1, 1), (2, 1), (2, 2)], M1),
    ([(0, 1), (1, 0), (1, 1), (1, 2), (2, 0)], M2),
    ([(0, 0), (0, 1), (1, 1), (1, 2), (2, 1)], M3),
  ),

  # 5 squares - F' shape (chiral)
  "F5'": (
    ([(0, 2), (1, 0), (1, 1), (1, 2), (2, 1)], R0),
    ([(0, 1), (1, 0), (1, 1), (2, 1), (2, 2)], R1),
    ([(0, 1), (1, 0), (1, 1), (1, 2), (2, 0)], R2),
    ([(0, 0), (0, 1), (1, 1), (1, 2), (2, 1)], R3),
    ([(0, 1), (1, 0), (1, 1), (1, 2), (2, 2)], M0),
    ([(0, 1), (1, 1), (1, 2), (2, 0), (2, 1)], M1),
    ([(0, 0), (1, 0), (1, 1), (1, 2), (2, 1)], M2),
    ([(0, 1), (0, 2), (1, 0), (1, 1), (2, 1)], M3),
  ),

  # 5 squares - Z shape (chiral)
  "Z5": (
    ([(0, 2), (1, 0), (1, 1), (1, 2), (2, 0)], R0 | R2),
    ([(0, 0), (0, 1), (1, 1), (2, 1), (2, 2)], R1 | R3),
    ([(0, 0), (1, 0), (1, 1), (1, 2), (2, 2)], M0 | M2),
    ([(0, 1), (0, 2), (1, 1), (2, 0), (2, 1)], M1 | M3),
  ),

  # 5 squares - Z' shape (chiral)
  "Z5'": (
    ([(0, 0), (1, 0), (1, 1), (1, 2), (2, 2)], R0 | R2),
    ([(0, 1), (0, 2), (1, 1), (2, 0), (2, 1)], R1 | R3),
    ([(0, 2), (1, 0), (1, 1), (1, 2), (2, 0)], M0 | M2),
    ([(0, 0), (0, 1), (1, 1), (2, 1), (2, 2)], M1 | M3),
  ),

}

# aka names, map aka -> polyomino name
akas = {
  "L3": "V3",
  "S4": "Z4'",
  "R4": "L4'",
  "R5": "L5'",
  "Q5": "P5'",
}

# shape selection flags
shape_flags = {
  "ALL": (R0 | R1 | R2 | R3 | M0 | M1 | M2 | M3),
  "ONE_SIDED": (R0 | R1 | R2 | R3),
}

# collect pieces using the names in ps
# return a list of shapes, each shape is represented in all orientations
# if as_map=1, return a dict() mapping names to shapes
def shapes(ps, flags="ALL", as_map=0):
  if isinstance(ps, basestring):
    ps = ps.split()
  else:
    ps = list(ps)
  flags = shape_flags.get(flags, flags)
  s = list()
  for p in ps:
    s.append(list(x for (x, f) in polyominoes[akas.get(p, p)] if f & flags))
  return (dict(zip(ps, s)) if as_map else s)

# look up polyomino names from shape orientations
def names(ss, flags="ALL"):
  flags = shape_flags.get(flags, flags)
  for (k, v) in polyominoes.items():
    v = list(x for (x, f) in v if f & flags)
    if v in ss: yield k

# useful routines for computing orientations of shapes
def orientations(ss, flags="ALL", verbose=0, indent=""):
  flags = shape_flags.get(flags, flags)

  def normalise(ss):
    ss = list(ss)
    (mx, my) = (min(x for (x, y) in ss), min(y for (x, y) in ss))
    return tuple(sorted((x - mx, y - my) for (x, y) in ss))

  def mirror(ss):
    return normalise((-x, y) for (x, y) in ss)

  def rotate(ss):
    return normalise((y, -x) for (x, y) in ss)

  from collections import defaultdict
  from enigma import unpack, join

  # accumulate shapes by orientation
  d = defaultdict(list)

  # consider quarter turns
  ss = normalise(ss)
  ms = mirror(ss)
  for i in (0, 1, 2, 3):
    d[ss].append(i)
    d[ms].append(i + 4)

    # rotate the shape by 90 degrees
    if i == 3: break
    ss = rotate(ss)
    ms = rotate(ms)

  # collect the results in order
  rs = ((k, v) for (k, v) in d.items() if flags & sum(1 << x for x in v))
  rs = sorted(rs, key=unpack(lambda k, v: v[0]))

  if verbose:
    for (k, v) in rs:
      v = list("RM"[i] + "0123"[j] for (i, j) in sorted(divmod(x, 4) for x in v))
      printf("{indent}({k}, {v}),", k=list(k), v=join(v, sep=" | "))

  return tuple(k for (k, v) in rs)

# extend the available shapes using the <name> -> <shape> data in d
def extend(d):
  for (k, v) in d.items():
    polyominoes[k] = orientations(v)
    

# [see enigma321.py for simple fit() algorithm]

# fit using algorithm X
# (see: [ https://www.cs.mcgill.ca/~aassaf9/python/algorithm_x.html ])

# NOTE: input parameters X, soln are modified; Y is unmodified
def algorithmX(X, Y, soln):
  if not X:
    yield soln
  else:
    c = min(X.keys(), key=lambda k: len(X[k]))
    # copy X[c], as X is modified (could use sorted(X[c]) for stability)
    for r in list(X[c]):
      soln.append(r)

      # cols = select(X, Y, r)
      cols = list()
      for j in Y[r]:
        for i in X[j]:
          for k in Y[i]:
            if k != j:
              X[k].remove(i)
        cols.append(X.pop(j))

      for s in algorithmX(X, Y, soln): yield s

      # deselect(X, Y, r, cols)
      for j in reversed(Y[r]):
        X[j] = cols.pop()
        for i in X[j]:
          for k in Y[i]:
            if k != j:
              X[k].add(i)

      soln.pop()

# generate placements for piece <p> in an <x> x <y> grid, avoiding <holes>
# <p> is a sequence of possible orientations for the piece
# return the linear indices of the occupied squares
def placements(p, x, y, holes):
  for q in p:
    # try to place the piece at <x0>, <y0>
    for y0 in range(y):
      for x0 in range(x):
        ss = list()
        for (dx, dy) in q:
          (i, j) = (x0 + dx, y0 + dy)
          if not(i < x and j < y): break
          if (i, j) in holes: break
          ss.append(i + x * j)
        else:
          yield ss

def fit(ps, x, y, holes=set(), fn=None):

  # check the dimensions of the pieces
  assert not(sum(len(p[0]) for p in ps) + len(holes) > x * y), "Impossible!"

  # how to format the results
  # each grid is calculated as a linear list
  # default is to chunk it into rows
  if fn is None: fn = lambda g: list(chunk(g, x))

  # set up the matrix for algorithm X

  # for each piece, in each orientation, and each position we add a
  # row of the form:
  #
  #   <squares occupied> + <indicator for the piece>
  #
  # where:
  #   
  #   <squares> is a boolean <x> * <y> vector
  #   <indicator> is a boolean <len(ps)> vector
  #
  # (if we add no rows for a piece, then the problem is not soluble)

  Y = list()

  # n is the number of available pieces
  n = len(ps)
  xy = x * y
  for (i, p) in enumerate(ps, start=xy):
    k = 0
    for squares in placements(p, x, y, holes):
      squares.append(i)
      Y.append(squares)
      k += 1
    if k == 0: return

  # if there are holes make a final piece to fill them
  if holes:
    squares = list(i + x * j for (i, j) in holes)
    squares.append(xy + n)
    Y.append(squares)

  # set up X as a dict of sets
  X = dict((k, set()) for k in range(n + xy + bool(holes)))
  for (i, y) in enumerate(Y):
    for k in y:
      X[k].add(i)

  # find exact covers using algorithm X
  for rs in algorithmX(X, Y, list()):
    # produce the grid
    g = [None] * xy
    for r in rs:
      # label the grid
      # 0 = hole
      # 1+ = piece number
      k = (Y[r][-1] + 1 - xy) % (n + 1)
      for i in Y[r][:-1]:
        g[i] = k
    yield fn(g)

# pack rectangles into a grid
def rectpack(rs, x, y, holes=set(), fn=None):

  # turn the rectangles into shapes (in both orientations)
  ps = list()
  for (a, b) in rs:
    p = list()
    p.append(list((i, j) for i in range(a) for j in range(b)))
    if a != b: p.append(list((j, i) for i in range(a) for j in range(b)))
    ps.append(p)

  # attempt to fit the shapes into a square
  for g in fit(ps, x, y, holes, fn):
    yield g
  

# output a grid
_labels = "-123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
def output_grid(g, reverse=1, label=_labels, sep=" ", end=""):
  if reverse: g = reversed(g)
  for r in g:
    print("[", join((label[x] for x in r), sep=sep), "]")
  if end is not None: print(end)


if enigma._namecheck(__name__):

  from enigma import arg, args, printf

  r = arg("D", 0)

  if r == "F":
    vs = args(["F5"], 1)
    printf("# [F] compute shape data: {vs}\n", vs=join(vs, sep=" "))
    ps = shapes(vs)
    for (v, p) in zip(vs, shapes(vs)):
      printf('"{v}": (')
      orientations(p[0], verbose=1, indent="  ")
      printf("),\n")

  if r == "E":
    printf("[E] pentominoes into 8x8 rectangle, with 2x2 central hole\n")
    # fit the pentominoes into a 8x8 grid with a 2x2 hole in the centre
    # each solution appears 8x (mirrored and rotated)
    # 260 solutions [65 different]
    # (129168 solutions if O4 is added and the hole removed)
    ps = shapes("F5 I5 L5 N5 P5 T5 U5 V5 W5 X5 Y5 Z5")
    ps[0] = ps[0][:1] # remove duplicates by fixing the orientation of F5
    n = 0
    for g in fit(ps, 8, 8, [(3, 3), (3, 4), (4, 3), (4, 4)]):
      output_grid(g)
      n += 1
    print("[{n} solutions]".format(n=n))

  if r == "D":
    x = arg(20, 1, int)
    y = arg(3, 2, int)
    printf("[D] pentominoes into {x}x{y} rectangle\n")
    # fit the pentominoes into a rectangle
    # each solution appears 4x (mirrored and rotated)
    # 3x20 = 8 solutions [2 different]
    # 4x15 = 1472 solutions [368 different]
    # 5x12 = 4040 solutions [1010 different]
    # 6x10 = 9356 solutions [2339 different]
    ps = shapes("F5 I5 L5 N5 P5 T5 U5 V5 W5 X5 Y5 Z5")
    ps[7] = ps[7][:1] # remove duplicates by fixing the orientation of V5
    n = 0
    for g in fit(ps, x, y):
      output_grid(g)
      n += 1
    print("[{n} solutions]".format(n=n))

  if r == "C":
    # 5x V3s in a 4x4 grid with 1 hole
    printf("[C] 5x V3 in a 4x4 grid with 1 hole\n")
    (O1, V3) = shapes("O1 V3")
    for g in fit([V3] * 5, 4, 4, holes=[(0, 2)]):
      output_grid(g)

  if r == "B":
    # one sided shapes in a 5x5 grid
    ns = args("I2 I3 O4 I4 S4 L4 R4".split(), 1)
    printf("[B] one sided shapes in a 5x5 grid: {ns}\n", ns=join(ns, sep=" "))
    ps = shapes(ns, "ONE_SIDED")
    for g in fit(ps, 5, 5):
      output_grid(g)

  if r == "A":
    # two sided shapes in a 13x2 grid
    ns = args("O1 I2 I3 V3 O4 U5 Z4 L4".split(), 1)
    printf("[A] two sided shapes in a 13x2 grid: {ns}\n", ns=join(ns, sep=" "))
    ps = shapes(ns)
    for g in fit(ps, 13, 2):
      output_grid(g)
