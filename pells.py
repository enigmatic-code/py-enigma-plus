#! python3
# -*- mode: Python; python-indent-offset: 2; coding: utf-8 -*-

# solving Pell type equations via continued fractions
#
# See:
# [ https://en.wikipedia.org/wiki/Pell%27s_equation ]
# [ "Solving the generalized Pell equation x^2 − D.y^2 = N", John P. Robertson, 2004 ]

from __future__ import print_function

from enigma import (
  enigma, irange, inf, is_square, sqrtf, sqrtc, gcd, div, divf, divc, multiply,
  prime_factor, divisors_pairs, sq, merge, as_int, cache, printf,
)

__author__ = "Jim Randell <jim.randell@gmail.com>"
__version__ = "2025-07-02"

pells = enigma.module(__name__)
verbose = ('v' in enigma._PY_ENIGMA)

######################################################################

# simple continued fractions:

# a continued fraction [a; b, (c, d)...] -> (i=a; nr=[b], rr=[c, d])

# continued fraction of sqrt(n)
@cache
def cf_sqrt(n):
  m = sqrtf(n)
  if m * m == n: return (m, [], [])
  (a, p, q, cs, end) = (m, 0, 1, [], m * 2)
  while a != end:
    p = a * q - p
    q = (n - p * p) // q
    a = (p + m) // q
    cs.append(a)
  return (m, [], cs)

# generate terms in the continued fraction
def cf_terms(cf, first=1):
  (i, nr, rr) = cf
  if first: yield i
  if nr:
    #yield from nr  #[Python 3]
    for x in nr: yield x  #[Python 2]
  if rr:
    while 1:
      #yield from rr  #[Python 3]
      for x in rr: yield x  #[Python 2]

# generate convergents for continued fraction cf
def cf_convergents(cf):
  (i, nr, rr) = cf
  (p0, q0, p, q) = (1, 0, i, 1)
  yield (p, q)
  for a in cf_terms(cf, first=0):
    (p0, q0, p, q) = (p, q, a * p + p0, a * q + q0)
    yield (p, q)

######################################################################

# solve Pell's equations using continued fractions:
# (only non-negative (X, Y) solutions are generated)

# find the fundamental solution for: X^2 - D.Y^2 = 1
def pells1_fundamental(D):
  cf = cf_sqrt(D)
  for (p, q) in cf_convergents(cf):
    if p * p - D * q * q == 1:
      return (p, q)

# find all (X, Y) solutions for: X^2 - D.Y^2 = 1
def pells1(D, trivial=1):
  if trivial: yield (1, 0)  # trivial solution
  # find the fundamental solution
  (x, y) = (x1, y1) = pells1_fundamental(D)
  (A, B, C) = (x1, D * y1, y1)
  if verbose: printf("[pells] [x^2 - {D}y^2 = 1] (x', y') = ({A}x + {B}y, {C}x + {A}y); (x1, y1) = ({x1}, {y1})")
  # generate all solutions
  while 1:
    yield (x, y)
    (x, y) = (A * x + B * y, A * y + C * x)

# find fundamental solution for: X^2 - D.Y^2 = -1 [D > 0, non-square]
def pells1n_fundamental(D):
  cf = cf_sqrt(D)
  # must have an odd period
  if len(cf[2]) % 2 == 0: return None
  # find the fundamental solution
  for (p, q) in cf_convergents(cf):
    if p * p - D * q * q == -1:
      return (p, q)

# find all (X, Y) solutions for: X^2 - D.Y^2 = -1
def pells1n(D):
  # find the fundamental solution
  s = pells1n_fundamental(D)
  if s is None: return
  (x, y) = (x1, y1) = s
  # generate all solutions
  (A, B, C) = (x1 * x1 + D * y1 * y1, 2 * D * x1 * y1, 2 * x1 * y1)
  if verbose: printf("[pells] [x^2 - {D}y^2 = -1] (x', y') = ({A}x + {B}y, {C}x + {A}y); (x1, y1) = ({x1}, {y1})")
  while 1:
    yield (x, y)
    (x, y) = (A * x + B * y, A * y + C * x)

# generate a solution family based on D, (x, y), (u, v)
def pells_sol(D, xy, uv):
  ((x, y), (u, v)) = (xy, uv)
  (A, B, C) = (u, D * v, v)
  if verbose: printf("[pells] [x^2 - {D}y^2 = {N}] (x', y') = ({A}x + {B}y, {C}x + {A}y); (x1, y1) = ({x}, {y})", N=sq(x) - D * sq(y))
  while 1:
    yield (x, y)
    (x, y) = (A * x + B * y, C * x + A * y)

# find all (X, Y) solutions for: X^2 - D.Y^2 = N [D > 0, non-square; N != 0]
def pellsN(D, N):
  # find the fundamental solution to the resolvant: X^2 - D.Y^2 = 1
  (u, v) = pells1_fundamental(D)

  # the brute force approach is probably OK for now,
  # but there are more sophisticated approaches [see: LMM algorithm]
  if N > 0:
    (a, b) = (0, sqrtf(divc(N * (u - 1), 2 * D)))
  else:
    (a, b) = (sqrtc(divc(-N, D)), sqrtf(divc(-N * (u + 1), 2 * D)))
  if b - a > 10000: printf("pells: WARNING: attempting brute force y = [{a} .. {b}]")

  # find solution families
  sols = list()
  for y in irange(a, b):
    x = is_square(N + D * y * y)
    if x is not None:
      sols.append(pells_sol(D, (x, y), (u, v)))
      # find minimal positive equivalent solution for (-x, y)
      (X, Y) = (x * u - y * v * D, x * v - y * u)
      if X < 0: (X, Y) = (-X, -Y)
      if (X, Y) != (x, y):
        sols.append(pells_sol(D, (X, Y), (u, v)))

  return merge(sols, uniq=1)

######################################################################

# find solutions to the quadratic Diophantine equation: a.X^2 + b.Y^2 = c

# a.X^2 + b.Y^2 = 0 [a > 0]
def _diop_quad_c0(a, b):
  # X=0 Y=0 is the first solution
  yield (0, 0)
  # if a and b have different signs there may be further solutions
  if b > 0: return
  (a, b) = (is_square(a), is_square(-b))
  if a is None or b is None: return
  for t in irange(1, inf):
    yield (t * b, t * a)

# a.X^2 + b.Y^2 = c [a > 0, b > 0, c > 0]
def _diop_quad_bp_old(a, b, c):
  # a brute force search
  X = 0
  while True:
    r = c - a * X * X
    if r < 0: break
    (Y2, z) = divmod(r, b)
    if z == 0:
      Y = sqrtf(Y2)
      if Y * Y == Y2: yield (X, Y)
    X += 1

# a.X^2 + b.Y^2 = c [a > 0, b > 0, c > 0]
def _diop_quad_bp_new(a, b, c):
  # start by considering possible Y values
  Y0 = sqrtf(divf(c, b))
  Y1 = sqrtc(divc(a * c, b * (a + b)))
  #printf("Y = [{Y0}, {Y1}]")
  X = None
  for Y in irange(Y0, Y1, step=-1):
    r = c - b * Y * Y
    (X2, z) = divmod(r, a)
    if z == 0:
      X = sqrtf(X2)
      if X * X == X2: yield (X, Y)
  # and then consider possible X values
  if X is None:
    X = sqrtf(divf(b * c, a * (a + b)))
  else:
    X += 1
  while True:
    r = c - a * X * X
    if r < 0: break
    (Y2, z) = divmod(r, b)
    if z == 0:
      Y = sqrtf(Y2)
      if Y * Y == Y2: yield (X, Y)
    X += 1

_diop_quad_bp = _diop_quad_bp_new

# X^2 - (dY)^2 = N [d > 0, N != 0]
def _diop_quad_d2(d, N):
  if N > 0:
    pqs = list(divisors_pairs(N))
    for (p, q) in pqs[::-1]:
      (X, Y) = (div(p + q, 2), div(q - p, 2 * d))
      if X is not None and Y is not None:
        yield (X, Y)
  else:
    pqs = list(divisors_pairs(-N))
    for (p, q) in pqs[::-1]:
      (X, Y) = (div(q - p, 2), div(p + q, 2 * d))
      if X is not None and Y is not None:
        yield (X, Y)

# X^2 - D.Y^2 = N [D > 0, non-square; N != 0]
def _diop_pells(D, N):
  if N == 1: return pells1(D)
  if N == -1: return pells1n(D)
  return pellsN(D, N)

# X^2 - D.Y^2 = N [D > 0; N != 0]
def _diop_quad_a1(D, N):
  d = is_square(D)
  if d is not None: return _diop_quad_d2(d, N)
  return _diop_pells(D, N)

def _diop_empty():
  if 0: yield None

# results (X, Y) for increasing X
def diop_quad(a, b, c, maxC=10000, validate=0):
  if validate: (a, b, c) = map(as_int, (a, b, c))
  if a == 0 or b == 0: raise ValueError("diop_quad: invalid equation")
  if a < 0: (a, b, c) = (-a, -b, -c)
  g = gcd(a, b)
  if c % g != 0: return _diop_empty()
  g = gcd(g, c)
  if g > 1: (a, b, c) = (a // g, b // g, c // g)
  if c == 0: return _diop_quad_c0(a, b)
  if b > 0: return _diop_quad_bp(a, b, c)
  if a == 1: return _diop_quad_a1(-b, c)

  # find a multiplier m, such that m.a is a square
  m = multiply(p for (p, e) in prime_factor(a) if e % 2 == 1)
  r = is_square(m * a)
  def _fn():
    for (i, (X, Y)) in enumerate(_diop_quad_a1(m * -b, m * c), start=1):
      #printf("[diop_quad: considering: X={X} Y={Y}]")
      X = div(X, r)
      if X is not None: yield (X, Y)
      if i == maxC: printf("diop_quad: WARNING: terminating search (after {maxC} candidates)"); break
  return _fn()

######################################################################

if enigma._namecheck(__name__):
  from enigma import (arg, number as num)

  (a, b, c, N, v) = (arg(-13, 0, num), arg(1, 1, num), arg(4, 2, num), arg(20, 3, num), arg(0, 4, num))
  printf("solving: {a}.X^2 + {b}.Y^2 = {c}")
  if v: verbose = 1
  for (i, (X, Y)) in enumerate(diop_quad(a, b, c, validate=1), start=1):
    r = a * X * X + b * Y * Y
    printf("[{i}] X={X} Y={Y} -> {r}")
    if i == N: printf("[first {N} solutions]"); break
  else:
    printf("[all solutions]")
