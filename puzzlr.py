#! python3
# -*- mode: Python; python-indent-offset: 2; coding: utf-8 -*-

from __future__ import print_function

# solve "New Scientist" Puzzlr puzzles

from enigma import (
  defaultdict,
  irange, flatten, chunk, subsets, unzip, chain, update, ifirst,
  divisors_pairs, cproduct, singleton, exact_cover, seq_all_same_r,
  wrap, fcompose, nl, join, fmts, printf,
  make_namespace, filter2, basestring, base2int, fail, warn, args, lazy_import, _namecheck
)

__author__ = "Jim Randell <jim.randell@gmail.com>"
__version__ = "2026-09-05"

###############################################################################

# provide defaults in the module

def resolve(name, frame):
  if name in frame.f_locals: return frame.f_locals.get(name)
  if name in frame.f_globals: return frame.f_globals.get(name)
  # (we could track back through previous frames)
  raise ValueError("can\'t resolve name " + repr(name))

class Defaults(object):

  def __init__(self, defs, fns):
    self.defs = defs
    self.fns = fns
    # also store a scope for resolving function names
    getframe = lazy_import("sys._getframe")
    self.frame = getframe(2)
    # has data been loaded?
    self.loaded = 0

  # load defaults from environment, command line
  def load(self):
    if self.loaded: return
    # (we could potentially load defaults from a preferences file here)
    # load any defaults defined in $PUZZLR_DEFAULTS
    self.load_env("PUZZLR_DEFAULTS")
    self.load_argv()
    self.loaded = 1

  # load defaults from an environment variable
  def load_env(self, var):
    os = lazy_import('os')
    s = os.getenv(var)
    if not s: return
    kw = dict(map(str.strip, str.split(x, '=')) for x in str.split(s, ';'))
    self.set(kw)

  # load defaults from sys.argv (which is modified)
  def load_argv(self):
    sys = lazy_import('sys')
    argv = sys.argv
    xs = list()
    for (i, x) in enumerate(argv):
      if x.startswith("--"):
        self.set(dict([x[2:].split('=')]))
        xs.insert(0, i)
    for i in xs: del argv[i]

  def set(self, defs, fns=None):
    self.defs.update(defs)
    if fns: self.fns.update(fns)

  def get(self, k, v=None):
    if not self.loaded: self.load()
    # if no value is specified fetch the stored value
    if v is None: v = self.defs.get(k)
    # do we need to parse a string?
    if isinstance(v, basestring):
      # is there a parse function?
      fn = self.fns.get(k)
      if isinstance(fn, basestring):
        fn = resolve(fn, self.frame)
      if callable(fn):
        v = fn(v, self.frame)
    #printf("defaults.get {k!r} -> {v!r}")
    return v

# parse strings to appropriate values
def list_of_str(v, frame=None): return str.split(v, ',')
def list_of_fn(v, frame=None): return list(resolve(str.strip(fn), frame) for fn in str.split(v, ','))

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

  # output solution
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

  ascii = output_ascii = output

  # plot a puzzle (using plot.py, if available)
  def output_plot(grid, rows, cols, sol=None, w=None):
    (W, H) = (len(grid[0]), len(grid))

    Plot = lazy_import('plot.Plot')

    p = Plot(width=500, height=500, xscale=64.0, yscale=-64.0, xoffset=0.40625, yoffset=-6.7345)

    # plot the cells and numbers
    font = defaults.get('plot.font')
    for y in irange(H):
      for x in irange(W):
        n = grid[y][x]
        p.line((x, y, x + 1, y, x + 1, y + 1, x, y + 1, x, y), width=0, tag=3)
        p.label((x + 0.5, y + 0.5), str(n), font=font, tag=2)
        p.label((x + 0.5, H + 0.15), u"\u2193", anchor="n", font=font, tag=2) # DOWN ARROW
        p.label((x + 0.5, H + 1.15), str(cols[x]), anchor="s", font=font, tag=2)
      p.label((W + 0.15, y + 0.5), u"\u2192", anchor="w", font=font, tag=2) # LEFT ARROW
      p.label((W + 1.2, y + 0.5), str(rows[y]), anchor="e", font=font, tag=2)

    # bounding box
    p.line((0, 0, W, 0, W, H, 0, H, 0, 0), width=4, tag=4)

    if sol:
      for y in irange(H):
        for x in irange(W):
          if sol[y][x]:
            p.circle((x + 0.5, y + 0.5), 0.4, fill=None, outline="red", width=4, tag=1)

    p.display()

  plot = output_plot

  # set defaults
  defaults = Defaults(
    { 'output': [output], 'plot.font': ("Helvetica", "22", "bold") },
    { 'output': list_of_fn, 'plot.font': list_of_str },
  )

  def run(grid, rows, cols, w=None, first=None, output=None):
    output_fns = defaults.get('output', output)
    # solve the puzzle
    sols = solve(grid, rows, cols)
    if first: sols = ifirst(sols, count=first)
    for sol in sols:
      # output solution
      for fn in output_fns:
        fn(grid, rows, cols, sol, w)

  # argv = ("<row sums>", "<col sums>", "<grid>")
  def run_command_line(argv):
    fail(len(argv) != 3, "invalid argv")
    # translate any strings into lists of numbers (or None for "-")
    fn = lambda x: (None if x == '-' else base2int(x))
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

      # return solved rows
      rows = list()
      for y in irange(H):
        rows.append(list(m[i] for i in irange(y * W, (y + 1) * W - 1)))
      yield rows

  # output the rectangles in a grid
  def output(grid, sol):
    fmt = fmts("02d")
    for ns in sol:
      printf("[ {xs} ]", xs=join(ns, fn=fmt, sep=" "))
    printf()

  output_rects = output

  # plot a puzzle (using plot.py, if available)
  def output_plot(grid, sol=None):
    (W, H) = (len(grid[0]), len(grid))

    Plot = lazy_import('plot.Plot')

    p = Plot(width=680, height=800, xscale=64.0, yscale=-64.0, xoffset=0.40625, yoffset=-11.5)

    # plot the cells and numbers
    font = defaults.get('plot.font')
    for y in irange(0, H):
      p.line((0, y, W, y), width=0, dash=(2, 4), tag=3)
    for x in irange(0, W):
      p.line((x, 0, x, H), width=0, dash=(2, 4), tag=3)
    for (y, x) in cproduct([irange(H), irange(W)]):
      n = grid[y][x]
      if n != 0:
        p.label((x + 0.5, y + 0.5), str(n), font=font, tag=2)

    # bounding box
    p.line((0, 0, W, 0, W, H, 0, H, 0, 0), width=4, tag=4)

    if sol:
      # plot solution
      rect = dict()
      for y in range(H):
        for x in range(W):
          n = sol[y][x]
          if n not in rect:
            # top-left
            rect[n] = [(x, y), None]
          else:
            # other, last will be bottom-right
            rect[n][1] = (x + 1, y + 1)
      for (n, ((x0, y0), (x1, y1))) in rect.items():
        p.line((x0, y0, x1, y0, x1, y1, x0, y1, x0, y0), width=4, tag=1)

    p.display()

  plot = output_plot

  # pretty print (using ASCII art) [contributed by Ruud van der Ham]
  @wrap(fcompose(join, print))  # join all the bits of output together
  def output_ascii(grid, rows):
    # turn the solution into a dict (to allow out of range indexing)
    (Y, X) = (len(grid), len(grid[0]))
    sol = dict(((y, x), rows[y][x]) for y in irange(Y) for x in irange(X))
    # generate bits of output
    fmt = fmts('2d')
    for y in irange(0, Y):
      for x in irange(0, X):
        v = (sol.get((y - 1, x - 1)) == sol.get((y, x - 1)) and sol.get((y - 1, x)) == sol.get((y, x)))
        h = (sol.get((y - 1, x - 1)) == sol.get((y - 1, x)) and sol.get((y, x - 1)) == sol.get((y, x)))
        yield ((' ' if h else '|') if v else ('-' if h else '+'))
        yield ('  ' if sol.get((y - 1, x)) == sol.get((y, x)) else '--')
      yield nl
      if y != Y:
        for x in irange(0, X):
          yield (' ' if sol.get((y, x - 1)) == sol.get((y, x)) else '|')
          if x != X:
            n = grid[y][x]
            yield (' .' if n == 0 else fmt(n))
      yield nl

  ascii = output_ascii

  defaults = Defaults(
    { 'output': [output_rects], 'plot.font': ("Helvetica", "22", "bold") },
    { 'output': list_of_fn, 'plot.font': list_of_str },
  )

  def run(grid, output=None):
    output_fns = defaults.get('output', output)
    # solve the puzzle
    for sol in solve(grid):
      # output solution
      for fn in output_fns:
        fn(grid, sol)

  # argv = ("<row>", "<row>", ...)
  def run_command_line(argv):
    # process any options "--<arg>=<val>"
    (optv, argv) = filter2((lambda x: x.startswith("--")), argv)
    for opt in optv:
      kw = dict([opt[2:].split('=')])
      defaults.set(kw)
    # translate any strings into lists of numbers
    argv = list((tuple(map(base2int, arg.split())) if isinstance(arg, basestring) else arg) for arg in argv)
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
