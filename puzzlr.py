#! python3
# -*- mode: Python; python-indent-offset: 2; coding: utf-8 -*-

from __future__ import print_function

# 

from enigma import (
  defaultdict,
  irange, flatten, chunk, subsets, unzip, chain, update, ifirst,
  divisors_pairs, cproduct, singleton, exact_cover, seq_all_same_r,
  join, fmts, printf,
  make_namespace, call, basestring, fail, warn, args, _namecheck
)

###############################################################################

# solver for "Plus Side"/"Number Cross" puzzles (e.g. "Puzzlr 1", "Puzzlr 3")

def __plus_side():

  # solve the puzzle in grid <grid> with row and col sums <rows>, <cols>
  def solve(grid, rows, cols):
    n = len(grid[0])
    # collect the numbers in the grid
    ns = flatten(grid)
    # make a corresponding grid of indices
    g = list(chunk(irange(len(ns)), n))
    # for each row/col construct possible subsets
    d = dict()
    for (ks, t) in chain(zip(g, rows), zip(unzip(g), cols)):
      if t is None: continue
      ss = list()
      for js in subsets(ks):
        if sum(ns[j] for j in js) == t:
          ss.append(js)
      d[ks] = ss
    # find solutions
    for rs in _solve(d, [None] * len(ns)):
      yield list(chunk(rs, n))

  def _solve(d, ss):
    # are we done?
    if not d:
      yield ss
    else:
      # consider a line with the fewest options
      ks = min(d.keys(), key=(lambda ks: len(d[ks])))
      for vs in d[ks]:
        # allocate the indices in ks
        ss_ = update(ss, ks, (k in vs for k in ks))
        # make a new map without <ks> and without any conflicting subsets
        d_ = dict()
        for (k, xs) in d.items():
          if k == ks: continue
          xs_ = list()
          for x in xs:
            if any((j in x) != ss_[j] for j in k if ss_[j] is not None): continue
            xs_.append(x)
          if not xs_: break  # no viable subsets remaining
          d_[k] = xs_
        else:
          # solve remaining viable subsets
          #yield from _solve(d_, ss_)  #[Python 3]
          for z in _solve(d_, ss_): yield z  #[Python 2]

  # solve the puzzle and output solutions
  def output(grid, rows, cols, sol, w=None):
    if w is None: w = max(len(str(x)) for x in chain(flatten(grid), rows, cols) if x is not None)
    (sf, null) = (join([w, "d"]), '-' * w)
    # output solution
    fmt = lambda x: (fmts(sf, x) if x is not None else null)
    for (row, vs, t) in zip(grid, sol, rows):
      xs = list((r if v else None) for (r, v) in zip(row, vs))
      printf("[ {xs} ] = {t}", xs=join(xs, fn=fmt, sep=" "), t=fmt(t))
    printf("= {cols}", cols=join(cols, fn=fmt, sep=" "))
    printf()
    #print(ss)

  def run(grid, rows, cols, w=None, first=None):
    # solve the puzzle
    sols = solve(grid, rows, cols)
    if first: sols = ifirst(sols, count=first)
    for sol in sols:
      # output solution
      output(grid, rows, cols, sol, w)

  # argv = ("<row sums>", "<col sums>", "<grid>")
  def run_command_line(argv):
    fail(len(argv) != 3, "invalid argv")
    # translate any strings into lists of numbers (or None for "-")
    fn = lambda x: (None if x == '-' else int(x))
    argv = list((tuple(map(fn, arg.split())) if isinstance(arg, basestring) else arg) for arg in argv)
    (rows, cols, grid) = argv
    fail(len(rows) * len(cols) != len(grid), "invalid grid spec")
    grid = list(chunk(grid, len(cols)))
    run(grid, rows, cols)

  # return exported names
  return locals()

plus_side = make_namespace('plus_side', __plus_side())


###############################################################################

# Solver for "Block Universe"/"Shikaku" puzzles (e.g. "Puzzlr 2", "Puzzlr 4")

def __block_universe():

  # return indices for a (w, h) rectangle @ (x, y)
  def rect(W, H, w, h, x, y):
    rs = set()
    i = x + y * W
    for _ in irange(h):
      rs.update(irange(i, i + w - 1))
      i += W
    return rs

  # solve (and output) the grid
  def solve(grid):
    # record the grid size
    (W, H) = (len(grid[0]), len(grid))
    # and then turn it into a linear list
    g = flatten(grid)
    # map positions to rectangle size
    pos = dict((i, v) for (i, v) in enumerate(g) if v > 0)
    ps = set(pos.keys())
    # consider possible rectangle area
    ss = defaultdict(list)
    for n in set(g):
      if n == 0: continue
      # consider possible rectangles
      for (w, h) in divisors_pairs(n, every=1):
        # possible (x, y) locations of a (w, h) rectangle
        for (x, y) in cproduct([irange(0, W - w), irange(0, H - h)]):
          rs = rect(W, H, w, h, x, y)
          # does it hit a single numbered cell?
          p = singleton(ps.intersection(rs))
          # is it the right size?
          if p is None or n != pos[p]: continue
          ss[p].append(rs)

    # form the grid into rectangles
    for rs in exact_cover(ss.values()):
      # map positions to rectangles
      m = dict()
      for (i, ps) in enumerate(rs):
        for p in ps:
          m[p] = i

      # output the rectangles in a grid
      fmt = fmts("02d")
      for y in irange(H):
        ns = list(m[i] for i in irange(y * W, (y + 1) * W - 1))
        printf("[ {xs} ]", xs=join(ns, fn=fmt, sep=" "))
      printf()

  run = solve

  # argv = ("<row>", "<row>", ...)
  def run_command_line(argv):
    # translate any strings into lists of numbers
    argv = list((tuple(map(int, arg.split())) if isinstance(arg, basestring) else arg) for arg in argv)
    rows = argv
    # check rows are all the same length
    r = seq_all_same_r(map(len, rows))
    fail(r.empty or not r.same, "invalid grid spec")
    # warn if the size of the rectangles isn't the same as the area of the grid
    ns = flatten(rows)
    warn(sum(ns) != len(ns), "area mismatch")
    run(rows)

  # return exported names
  return locals()

block_universe = make_namespace('block_universe', __block_universe())

###############################################################################

# allow puzzles to be solved from the command line
if _namecheck(__name__):
  argv = args([], 0)
  if argv:
    # extract the command
    cmd = argv.pop(0)
    ns = locals().get(cmd)
    fail(ns is None, "unrecognised command")
    fn = getattr(ns, 'run_command_line', None)
    fail(fn is None, "invalid command")
    # call the function
    fn(argv)
  #else: help()
