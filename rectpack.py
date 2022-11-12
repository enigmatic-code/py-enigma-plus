#! python3
# -*- mode: Python; py-indent-offset: 2; -*-

from __future__ import print_function

# routines for packing rectangles
#
# rectangles are represented by (w, h) dimensions
#
# a packing is returned as a list of (x, y, w, h), indicating the
# (x, y) co-ordinate of the top left-hand corner, and the dimensions
# of the rectangle (which may be in a different orientation to the
# corresponding input rectangle)

from enigma import irange, unpack, uniq, join, printf

# sort rectangles <rs> by area (default = largest to smallest)
def by_area(rs, reverse=0):
  return sorted(rs, key=unpack(lambda w, h: w * h), reverse=not(reverse))

# by_area -> smallest to largest, largest to smallest
by_area_stol = lambda rs: by_area(rs, reverse=1)
by_area_ltos = by_area

# determine if rectangle <r> overlaps with a rectangle in <ps>
# return the index in <rs> of an overlapping rectangle, or -1
def overlap(r, ps):
  (i1, j1, p1, q1) = r
  for (k, (i2, j2, p2, q2)) in enumerate(ps):
    if i1 < i2 + p2 and i2 < i1 + p1 and j1 < j2 + q2 and j2 < j1 + q1:
      return k
  return -1

# fit the rectangles <rs> into a <n> x <m> grid (loose packing)
# n, m - the dimensions of the grid
# rs - dimensions of the rectangles [(w, h), ...]
# ps - positions of the rectangles [(x, y, w, h), ...]
def pack_loose(n, m, rs, ps=[]):
  # are we done?
  if not rs:
    yield ps
  else:
    # try to fit the next rectangle into the grid
    r = rs[0]
    for (p, q) in {r, r[::-1]}:
      # consider possible locations for the rectangle
      i = j = 0
      while True:
        if i + p > n:
          i = 0
          j += 1
        if j + q > m: break
        # does this position overlap with any placed rectangles?
        r = (i, j, p, q)
        k = overlap(r, ps)
        if k == -1:
          # try to place the remaining rectangles
          for z in pack_loose(n, m, rs[1:], ps + [r]): yield z
          i += 1
        else:
          (x, y, w, h) = ps[k]
          i = x + w

# find the first empty square, starting at (i, j)
def empty(n, m, ps, i=0, j=0):
  while True:
    if i >= n:
      j += i // n
      i %= n
    if j >= m: break
    k = overlap((i, j, 1, 1), ps)
    if k == -1:
      return (i, j)
    else:
      (x, y, w, h) = ps[k]
      i = x + w

# fit the rectangles <rs> into an <n> x <m> grid (tight packing)
# n, m = the dimensions of the grid
# rs = dimensions of the rectangles [(w, h), ...]
# ps = positions of the rectangles [(x, y, w, h), ...]
# i, j = position to start looking for empty squares
def pack_tight(n, m, rs, ps=[], i=0, j=0):
  # are we done?
  if not rs:
    yield ps
  else:
    # find an empty square
    (i, j) = empty(n, m, ps, i, j)
    # fit one of the remaining rectangles there
    for (k, r) in enumerate(rs):
      for (p, q) in {r, r[::-1]}:
        if not(i + p > n or j + q > m):
          r = (i, j, p, q)
          if overlap(r, ps) == -1:
            # and try to place the remaining rectangles
            for z in pack_tight(n, m, rs[:k] + rs[k + 1:], ps + [r], i + p, j): yield z

def pack(n, m, rs, packer=pack_tight, order=by_area):
  # do some quick checks to look for impossible scenarios
  # total area
  if sum(w * h for (w, h) in rs) > n * m: return ()
  # check all rectangles fit in the grid
  (dmin, dmax) = sorted([n, m])
  if any(max(r) > dmax or min(r) > dmin for r in rs): return ()
  # stack rectangles with min dimension > 1/2 min dimension of grid [suggested by frits]
  if sum(x for x in map(min, rs) if 2 * x > dmin) > dmax: return ()
  # order the rectangles (if required)
  if callable(order): rs = order(rs)
  # determine the packer
  if not callable(packer): packer = locals().get(packer)
  # do the packing
  return packer(n, m, rs)

# reflect a solution about vertical / horizontal axis
soln = lambda s: tuple(sorted(s))
reflect_v = lambda X, Y, s: soln((X - x - w, y, w, h) for (x, y, w, h) in s)
reflect_h = lambda X, Y, s: soln((x, Y - y - h, w, h) for (x, y, w, h) in s)
# rotate a solution 90 degrees (will change shape of grid unless X = Y)
rotate = lambda X, Y, s: soln((Y - y - h, x, h, w) for (x, y, w, h) in s)

# determine canonical form of a packing
def canonical(n, m, s):
  s0 = soln(s)
  s1 = reflect_v(n, m, s0)
  s2 = reflect_h(n, m, s0)
  s3 = reflect_h(n, m, s1)
  if n != m:
    # a non-square rectangle has a symmetry group of order 4
    return min(s0, s1, s2, s3)
  else:
    # a square has a symmetry group of order 8
    r0 = rotate(n, m, s0)
    r1 = rotate(n, m, s1)
    r2 = rotate(n, m, s2)
    r3 = rotate(n, m, s3)
    return min(s0, s1, s2, s3, r0, r1, r2, r3)

# generate symmetrically different packings
# this keeps track of all packings found, so may use a lot of memory
def pack_uniq(n, m, rs, packer=pack_tight, order=by_area, verbose=0):
  return uniq((canonical(n, m, s) for s in pack(n, m, rs, packer, order)), verbose=verbose)

# return an <n> x <m> grid containg rectangles <ps> = [(x, y, w, h) ...]
def make_grid(n, m, ps):
  # make an empty grid
  g = list([0] * n for _ in irange(1, m))
  # fill out the rectangles
  for (k, (x, y, p, q)) in enumerate(ps, start=1):
    for j in irange(y, y + q - 1):
      g[j][x : x + p] = [k] * p
  return g

# output a packed <n> x <m> grid (either from <ps> or <g>)
def output_grid(n, m, ps=None, g=None):
  if g is None: g = make_grid(n, m, ps)
  # output the packing
  k = max(max(r) for r in g)
  z = len(str(k))
  for r in g:
    printf("[ {r} ]", r=join((str(x or 0).zfill(z) for x in r), sep=' '))
