#! python3
# -*- mode: Python; py-indent-offset: 2; -*-

# routines for handling polyominoes

from __future__ import print_function

from enigma import (
  enigma, algorithmX, seq_all_same_r, basestring, chunk, unpack, join,
  arg, args, printf
)

__author__ = "Jim Randell <jim.randell@gmail.com>"
__version__ = "2024-07-01"

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

# shape data
polyominoes = dict()

# polyomino templates
template = dict()

# monominos
template[1] = [("O1", 0, "@")]

# dominos
template[2] = [("I2", 0, "@@")]

# trominoes
template[3] = [
  # (<name>, <chiral?>, <layout>)
  ("I3", 0, "@@@"),
  ("V3", 0, "@@|@ "),
]

# tetrominoes
template[4] = [
  # (<name>, <chiral?>, <layout>)
  ("O4", 0, "@@|@@"),
  ("I4", 0, "@@@@"),
  ("T4", 0, "@@@| @ "),
  ("Z4", 1, " @@|@@ "),
  ("L4", 1, "@@@|  @"),
]

# pentominoes
template[5] = [
  # (<name>, <chiral?>, <layout>)
  ("I5", 0, "@@@@@"),
  ("U5", 0, "@@@|@ @"),
  ("T5", 0, "@@@| @ | @ "),
  ("V5", 0, "@@@|@  |@  "),
  ("X5", 0, " @ |@@@| @ "),
  ("W5", 0, " @@|@@ |@  "),
  ("L5", 1, "@@@@|   @"),
  ("Y5", 1, "@@@@|  @ "),
  ("P5", 1, "@@@|@@ "),
  ("N5", 1, " @@@|@@  "),
  ("F5", 1, " @ |@@@|@  "),
  ("Z5", 1, "@  |@@@|  @"),
]

# hexominoes
template[6] = [
  # (<name>, <chiral?>, <layout>)
  ("A06", 0, "@@@|@@ |@  "), # "Kadon's A"
  ("C06", 0, "@@@@|@  @"),
  ("D06", 0, "@@@@| @@ "),
  ("E06", 0, "@ @|@@@| @ "),
  ("F06", 1, "   @|@@@@|  @ "), # "hi F"
  ("F16", 1, "   @|@@@@| @  "), # "low F"
  ("F26", 1, "  @ |@@@ |  @@"), # "hi 4"
  ("F36", 1, " @  |@@  | @@@"), # "low 4"
  ("G06", 1, "@@@ |@ @@"),
  ("H06", 1, "@ @|@@@|@  "),
  ("I06", 0, "@@@@@@"),
  ("J06", 1, "@@@|@ @|  @"),
  ("K06", 0, "@@ |@@@| @ "),
  ("L06", 1, "@@@@@|    @"),
  ("M06", 1, "@   |@@@ |  @@"),
  ("N06", 1, " @@@|@@@ "), # "short N"
  ("N16", 1, "   @@|@@@@ "), # "long N"
  ("O06", 0, "@@@|@@@"),
  ("P06", 1, "@@@@|@@  "),
  ("Q06", 1, "@@ |@@@|  @"),
  ("R06", 1, "@@@|@@ | @ "),
  ("S06", 1, "@@@  |  @@@"), # "long S"
  ("T06", 0, "@   |@@@@|@   "), # "long T"
  ("T16", 1, "@@@@|  @ |  @ "), # "short T"
  ("U06", 1, "@@@@|@ @ "),
  ("V06", 1, "@@@@|@   |@   "),
  ("W06", 1, " @@@|@@  |@   "), # "Wa"
  ("W16", 1, "@@  | @@ |  @@"), # "Wb"
  ("W26", 1, "@@  | @@@|  @ "), # "Wc"
  ("X06", 0, " @  |@@@@| @  "),
  ("X16", 1, "  @ |@@@@| @  "), # "italic X"
  ("Y06", 1, "@@@@@|   @ "), # "hi Y"
  ("Y16", 0, "@@@@@|  @  "), # "lo Y"
  ("Z06", 1, "@   |@@@@|   @"), # "long Z"
  ("Z16", 1, " @@@| @  |@@  "), # "short Z"
]

populated = set()

# populate the polyominoes data from templates
def populate(ks):
  for k in ks:
    if k in populated: continue
    for (n, f, s) in template[k]:
      # determine layout width
      ss = s.split('|')
      w = seq_all_same_r(len(x) for x in ss)
      assert w.same
      w = w.value
      # construct cooordinates
      ps = list()
      for (y, row) in enumerate(ss):
        for (x, c) in enumerate(row):
          if c != ' ': ps.append((x, y))
      assert len(ps) == k
      # populate the orientations
      polyominoes[n] = orientations(ps)
      # and the mirror image for chiral shapes
      if f:
        ps = list((w - x, y) for (x, y) in ps)
        polyominoes[n + "'"] = orientations(ps)
    populated.add(k)

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
  "R0": R0, "R1": R1, "R2": R2, "R3": R3,
  "M0": M0, "M1": M1, "M2": M2, "M3": M3,
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

  return tuple((ps, sum(1 << f for f in fs)) for (ps, fs) in rs)

# extend the available shapes using the <name> -> <shape> data in d
def extend(d):
  assert 0, "DEPRECATED!"
  for (k, v) in d.items():
    polyominoes[k] = orientations(v)

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
          if not (i < x and j < y): break
          if (i, j) in holes: break
          ss.append(i + x * j)
        else:
          yield ss

def fit(ps, x, y, holes=set(), fn=None):

  # check the dimensions of the pieces
  assert not (sum(len(p[0]) for p in ps) + len(holes) > x * y), "Impossible!"

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
    printf("[ {r} ]", r=join((label[x] for x in r), sep=sep))
  if end is not None: printf("{end}")

# compute the polyominoes
populate([1, 2, 3, 4, 5]) # 6 (hexominoes) are optional

if enigma._namecheck(__name__):

  r = arg("D", 0)

  if r == "G":
    vs = args(["*"], 1)

    # a grid of characters
    class Grid(object):
      def __init__(self):
        self.X = 0
        self.Y = 0
        self.lines = list()
      def place(self, xys, dx=0, dy=0, c='O'):
        for (x, y) in xys:
          x += dx
          y += dy
          while y >= self.Y:
            self.lines.append([' '] * self.X)
            self.Y += 1
          while x >= len(self.lines[y]):
            self.lines[y].append(' ')
            if x > self.X: self.X = x
          self.lines[y][x] = c
      def display(self, pre="", post="", start=None, end=None):
        if start is not None: printf("{start}")
        for r in g.lines[::-1]:
          printf("{pre}{r}{post}", r=str.join('', r))
        if end is not None: printf("{end}")

    populate([1, 2, 3, 4, 5, 6])
    if vs == ['*']: vs = list(polyominoes.keys())
    for (v, ps) in zip(vs, shapes(vs)): # flags="ONE_SIDED"
      printf("\"{v}\" [{n} orientations]", n=len(ps))
      g = Grid()
      for xys in ps:
        g.place(xys, dx=(g.X + 6 if g.X else 0))
      g.display(pre="  ", start="", end="\n")


  if r == "F":
    vs = args(["F5"], 1)
    printf("# [F] compute shape data: {vs}\n", vs=join(vs, sep=" "))
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
    printf("[{n} solutions]")

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
    printf("[{n} solutions]")

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
