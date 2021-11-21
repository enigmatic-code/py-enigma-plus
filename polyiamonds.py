#!/usr/bin/env python3 -t
# -*- mode: Python; py-indent-offset: 2; -*-

from __future__ import print_function

# 

from enigma import basestring, exact_cover, irange, unpack, peek, join, printf

# polyiamonds:
#
# the grid is layed out as follows:
#
#           /.
#          /.  \   /   \   /   \   /   \  /
#         /.    \ /     \ /     \ /     \/...
#        ------------------------------------
#       / \(0,3)/ \(1,3)/ \(2,3)/ \(3,3)/
#      /   \   /   \   /   \   /   \   /
#     /(0,2)\ /(1,2)\ /(2,2)\ /(3,2)\ /...
#    -------------------------------------
#   / \(0,1)/ \(1,1)/ \(2,1)/ \(3,1)/
#  /   \   /   \   /   \   /   \   /
# /(0,0)\ /(1,0)\ /(2,0)\ /(3,0)\ /...
# ------------------------------------

# R0 | R1 | R2 | R3 | R4 | R5 = rotated 0 | 60 | 120 | 180 | 240 | 300 degrees
(R0, R1, R2, R3, R4, R5) = (1, 2, 4, 8, 16, 32)
# corresponding reflections
(M0, M1, M2, M3, M4, M5) = (64, 128, 256, 512, 1024, 2048)

# precalculated polyiamonds (can also be generated from prototypes)
polyiamonds = {}

# prototypes of convex polyiamonds that fit in a 24-hex grid
protos = {
  # 1 cell (moniamond)
  "T1": [(0, 0)], # achiral
  # 2 cells (diamond)
  "D2": [(0, 0), (0, 1)], # achiral
  # 3 cells (triamond)
  "I3": [(0, 0), (0, 1), (1, 0)], # achiral
  # 4 cells (tetriamonds)
  "T4": [(0, 0), (0, 1), (0, 2), (1, 0)], # achiral
  "I4": [(0, 0), (0, 1), (1, 0), (1, 1)],
  "I4'": [(0, 1), (1, 0), (1, 1), (2, 0)],
  # 5 cells (pentiamonds)
  "I5": [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0)], # achiral
  # 6 cells (hexiamonds)
  "O6": [(0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2)], # achiral
  "I6": [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)],
  "I6'": [(0, 1), (1, 0), (1, 1), (2, 0), (2, 1), (3, 0)],
  # non-convex hexiamonds
  "C6": [(0, 0), (0, 1), (0, 2), (0, 3), (1, 2), (1, 3)], # achiral
  "E6": [(0, 1), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1)], # achiral
  "F6": [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)],
  "F6'": [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 2)],
  "G6": [(0, 1), (0, 2), (0, 3), (0, 4), (1, 0), (1, 1)],
  "G6'": [(0, 0), (0, 1), (0, 2), (0, 3), (1, 1), (1, 2)],
  "H6": [(0, 1), (1, 0), (1, 1), (1, 2), (1, 3), (2, 0)],
  "H6'": [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3)],
  "J6": [(0, 1), (1, 0), (1, 1), (2, 0), (2, 1), (2, 2)],
  "J6'": [(0, 1), (0, 2), (1, 0), (1, 1), (2, 0), (2, 1)],
  "P6": [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 0)],
  "P6'": [(0, 0), (0, 1), (1, 0), (1, 1), (1, 2), (2, 0)],
  "S6": [(0, 1), (1, 0), (1, 1), (1, 2), (1, 3), (2, 2)],
  "S6'": [(0, 2), (0, 3), (1, 1), (1, 2), (2, 0), (2, 1)],
  "V6": [(0, 0), (0, 1), (0, 2), (1, 0), (0, 3), (1, 1)], # achiral
  "X6": [(1, 0), (1, 1), (2, 0), (0, 3), (1, 2), (1, 3)], # achiral
  # 7 cells (heptiamond)
  "D7": [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2)], # achiral
  "I7": [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1), (3, 0)], # achiral
  # 8 cells (octiamond)
  "t8": [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2), (2, 0)], # achiral
  "d8": [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2), (1, 3)], # achiral
}

# aka names
akas = {}

# shape selection flags
shape_flags = {
  "ALL": (R0 | R1 | R2 | R3 | R4 | R5 | M0 | M1 | M2 | M3 | M4 | M5),
  "ONE_SIDED": (R0 | R1 | R2 | R3 | R4 | R5),
}

# make a shape from a prototype (and place it into polyiamonds)
def make_shape(name, proto=None):
  name = akas.get(name, name)
  if proto is None: proto = protos[name]
  v = orientations(proto, verbose=0)
  polyiamonds[name] = v
  return v

# collect pieces using the names in ps
def shapes(ps, flags="ALL", as_map=1):
  if isinstance(ps, basestring):
    ps = ps.split()
  else:
    ps = list(ps)
  flags = shape_flags.get(flags, flags)
  s = list()
  for p in ps:
    p = akas.get(p, p)
    v = polyiamonds.get(p)
    if v is None: v = make_shape(p)
    s.append(list(x for (x, f) in v if f & flags))
  return (dict(zip(ps, s)) if as_map else s)

# compute canonical orientations of a shape from a prototype
def orientations(proto, flags="ALL", verbose=0, indent=""):
  flags = shape_flags.get(flags, flags)

  # normalise position of cells
  def normalise(cs):
    cs = list(cs)
    mx = min(x for (x, y) in cs)
    my = min(y for (x, y) in cs)
    my -= my % 2
    return tuple(sorted((x - mx, y - my) for (x, y) in cs))

  # mirror in horizontal axis
  def mirror(cs):
    return normalise((x + (y + 1) // 2, -(y + 1)) for (x, y) in cs)

  # rotate anticlockwise
  def rotate(cs):
    return normalise((-1 - y // 2, 2 * x + y + 1) for (x, y) in cs)

  # add bit v to dict d at key k
  def add(d, k, v):
    v = 1 << v
    if k in d:
      d[k] += v
    else:
      d[k] = v

  # return bits set in v
  def bits(v):
    n = 0
    while v:
      if v & 1: yield n
      v >>= 1
      n += 1
      
  # accumulate shapes by orientation
  d = dict()
  cs = normalise(proto)
  ms = mirror(cs)
  for i in irange(0, 5):
    add(d, cs, i)
    add(d, ms, i + 6)
    if i == 5: break
    cs = rotate(cs)
    ms = rotate(ms)

  # collect the results in order
  rs = ((k, v) for (k, v) in d.items() if flags & v)
  rs = sorted(rs, key=unpack(lambda k, v: peek(bits(v))))

  if verbose:
    for (k, v) in rs:
      v = list("RM"[i] + "012345"[j] for (i, j) in (divmod(x, 6) for x in bits(v)))
      printf("{indent}({k}, {v}),", k=list(k), v=join(v, sep=" | "))
  
  return rs

# generate placements for piece <p> in grid <grid>
def placements(p, grid):
  xmax = max(x for (x, y) in grid)
  ymax = max(y for (x, y) in grid)
  for q in p:
    dx = dy = 0
    q_ = q
    while True:
      while True:
        if grid.issuperset(q_):
          yield q_
        # increase x
        dx += 1
        q_ = list((x + dx, y + dy) for (x, y) in q)
        if all(x > xmax for (x, y) in q_): break
      # increase y
      dx = 0
      dy += 2
      q_ = list((x + dx, y + dy) for (x, y) in q)
      if all(y > ymax for (x, y) in q_): break

# fit pieces <ps> into grid <grid>
# <start> is the starting label for the pieces
# <accept> is used to determine acceptable placements
def fit(ps, grid, start=1, accept=None):
  # check the dimensions of the pieces
  assert sum(len(p[0]) for p in ps) == len(grid)

  # create the sets for exact_cover
  sss = list()
  for p in ps:
    ss = list(filter(accept, placements(p, grid)))
    if not ss: return
    sss.append(ss)

  # solve the exact cover
  for rs in exact_cover(sss, grid):
    # return a map of grid cells to piece number
    g = dict()
    for (i, cs) in enumerate(rs, start=start):
      g.update((c, i) for c in cs)
    yield g

# output a grid
_labels = "-123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
def output_grid(g, label=_labels, end=""):
  ymin = min(y for (x, y) in g.keys())
  ymax = max(y for (x, y) in g.keys())
  for y in irange(ymax, ymin, step=-1):
    xmin = min(x for (x, y_) in g.keys() if y_ == y)
    xmax = max(x for (x, y_) in g.keys() if y_ == y)
    print(" " * (2 * xmin + (y + 1) // 2) + "[" + join((_labels[g.get((x, y), 0)] for x in irange(xmin, xmax)), sep=" ") + "]")
  if end is not None: print(end)

if __name__ == "__main__":

  from enigma import arg

  # can use this to provide orientation data to place in polyiamonds dict
  name = arg("T1", 0)
  orientations(protos[name], verbose=1)
