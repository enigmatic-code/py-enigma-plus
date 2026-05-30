#! python3
# -*- mode: Python; python-indent-offset: 2; coding: utf-8 -*-

# solving Pell type equations via continued fractions
#
# See:
# [ https://en.wikipedia.org/wiki/Pell%27s_equation ]
# [ "Solving the generalized Pell equation x^2 − D.y^2 = N", John P. Robertson, 2004 ]

"""
The primary interface to this package is:

  diop_quad(a, b, c)

which finds non-negative integer solutions to the equation:

  a.X^2 + b.Y^2 = c

(X, Y) solutions are generated in increasing X order.

the behaviour of the algorithm can be modified with the following variables:

  pells_threshold = threshold at which it switches to using the LMM
    algorithm to solve the general Pell's equation. Below this a brute
    force algorithm is used. (default = 100000)

  cornacchia_threshold = threshold at which it switches to Cornacchia's
    algorithm when a > 0, b > 0, c > 0. Below this value a brute force
    algorithm is used. (default = 500000)

  diop_quad_threshold = threshold at which candidate solutions are
    considered exhausted. (default = 10000)

Auxiliary routines which may be useful are:

  sqrtmod(a, m) = calculate modular square roots of <a> mod <m>

  sqrtmodp(a, p, k) = calculate modular square roots of <a> mod <p>^<k>,
    where <p> is prime.
"""

from __future__ import print_function

from enigma import (
  enigma, irange, inf, is_square, sqrtf, sqrtc, gcd, div, divf, divc, multiply, invmod,
  divisors_pairs, sq, rev, uniq, merge, multiset, cproduct, crt, identity, catch, as_int, cache, printf,
)

__author__ = "Jim Randell <jim.randell@gmail.com>"
__version__ = "2026-05-26"

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

# threshold above which we switch to LMM (or issue a warning)
pells_threshold = 100000

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

  fn = None
  if b - a > pells_threshold:
    rD = sqrtf(D)
    if N * N < D and rD * rD < D:
      if verbose: printf("[pells] switching to simplified LMM (rather then brute force y = [{a} .. {b}])")
      fn = pells_LMMs(D, N, rD)
    else:
      if verbose: printf("[pells] switching to LMM (rather than brute force y = [{a} .. {b}])")
      fn = pells_LMM(D, N, rD, u, v)
  if not fn:
    if verbose or b - a > pells_threshold: printf("[pells] attempting brute force y = [{a} .. {b}]")
    fn = pells_BF(D, N, a, b, u, v)

  # find solution families
  sols = (pells_sol(D, xy, (u, v)) for xy in uniq(fn))
  # merge solutions from all families
  return merge(sols, uniq=1)

# brute force
def pells_BF(D, N, a, b, u, v):
  for y in irange(a, b):
    x = is_square(N + D * y * y)
    if x is not None:
      yield (x, y)
      # find minimal positive equivalent solution for (-x, y)
      (X, Y) = (u * x - D * v * y, v * x - u * y)
      if X < 0: (X, Y) = (-X, -Y)
      if (X, Y) != (x, y):
        yield (X, Y)

# simplified LMM (for N^2 < D)
def pells_LMMs(D, N, rD=None):
  if rD is None: rD = sqrtf(D)
  fs = dict()
  for (f, _) in divisors_pairs(abs(N)):
    (m, r) = divmod(N, f * f)
    if r == 0:
      fs[m] = f
  (P, Q, G0, G1, B0, B1) = (0, 1, 0, 1, 1, 0)
  while 1:
    for _ in (0, 1):
      a = (P + rD) // Q
      (G0, G1, B0, B1) = (G1, a * G1 + G0, B1, a * B1 + B0)
      P = a * Q - P
      Q = (D - P * P) // Q
      f = fs.get(G1 * G1 - D * B1 * B1)
      if f is not None:
        yield (f * G1, f * B1)
    if Q == 1: break

######################################################################

# Tonelli-Shanks algorithm for modular square roots
# (used in Cornacchia's algorithm, and general LMM algorithm)

# we could use:
# [[ sqrtmod = lambda a, m, fs=None: enigma.poly_roots_mod.sqrtmod(a, m) ]]
# (which saves a lot of coding)

# but the following is more efficient for large numbers (and is based on sympy.ntheory.sqrt_mod_iter)

# we need an efficient factorisation implementation
prime_factor = enigma.prime_factor
# although for large numbers we might use:
#prime_factor = enigma.partial(enigma.prime_factor_h, ps=enigma.primes, end=1000000, mr=1)

# check for numbers that have modular square roots (i.e. is a quadratic residue)
# there are several functions that can do this: legendre(), jacobi(), kronecker()

# we use the Legendre symbol (a|p); a > 0; p is an odd prime
# = 0 if a = 0 (mod p)
# = 1 if a is a quadratic residue (mod p)
# = -1 otherwise
def legendre(a, p):
  r = pow(a, (p - 1) // 2, p)  # calculate using Euler's criteria
  return (-1 if r == p - 1 else r)

# check if a number is a quadratic residue
# [see also: gmpy2.legendre(), gmpy2.jacobi(), gmpy2.kronecker()]
qr_check = legendre

# Tonelli-Shanks algorithm for modular square roots
# returns x, such that x^2 = a (mod p)
# -x (mod p) is also a root
# p must be prime, with p % 8 = 1
def tonelli_shanks(a, p):
  assert p % 8 == 1

  # s = gmpy2.bitscan1(p - 1)
  # t = p >> s
  (s, t) = (0, p - 1)
  while t & 1 == 0:
    t >>= 1
    s += 1

  # find a non-quadratic residue
  if p % 12 == 5:
    d = 3
  elif p % 5 in {2, 3}:
    d = 5
  else:
    d = 6
    while d < p:
      if qr_check(d, p) == -1: break
      d += 1

  (A, D, m) = (pow(a, t, p), pow(d, t, p), 0)
  for i in irange(s):
    u = (A * pow(D, m, p)) % p
    u = pow(u, 1 << (s - 1 - i), p)
    if u % p == p - 1:
      m += (1 << i)
  return (pow(a, (t + 1) // 2, p) * pow(D, m // 2, p)) % p

# now find all solutions to x^2 = a (mod p^k)

# case 1: a and p are coprime (there are at most 4 roots)
# returns a list of the roots
def _sqrtmodp1(a, p, k):
  assert a % p != 0

  pk = p**k
  a %= pk
  if p == 2:
    if a % 8 != 1: return
    if k < 4: return list(irange(1, pk, step=2))
    # hensel lift
    r = 1
    for n in irange(3, k - 1):
      if ((r*r - a) >> n) & 1:
        r += 1 << (n - 1)
    h = r + (1 << (k - 1))
    return [r, pk - r, h % pk, -h % pk]

  if qr_check(a, p) != 1: return []
  if p % 4 == 3:
    r = pow(a, (p + 1) // 4, p)
  elif p % 8 == 5:
    r = pow(a, (p + 3) // 8, p)
    if pow(r, 2, p) != a % p:
      r = (r * pow(2, (p - 1) // 4, p)) % p
  else:
    r = tonelli_shanks(a, p)

  if k > 1:
    # hensel lift
    px = p
    for _ in irange(k.bit_length() - 1):
      px *= px
      r = (r - (r*r - a) * invmod(2*r, px)) % px
    if k & (k - 1):
      # k is not a power of 2
      r = (r - (r*r - a) * invmod(2*r, px)) % pk

  return [r, pk - r]

# case 2: a % p = 0, there may be many solutions
# returns an iterable or None
def _sqrtmodp2(a, p, k):
  assert a % p == 0
  pk = p**k
  a %= pk

  if a == 0: return irange(0, pk - 1, step=p**((k + 1) // 2))

  m = 0
  while a % p == 0:
    (a, m) = (a // p, m + 1)
  if m % 2 == 1: return None
  rs = _sqrtmodp1(a, p, k - m)
  if not rs: return None
  m //= 2
  return (x for r in rs for x in irange(r * p**m, pk - 1, step=p**(k - m)))

# combine the cases
# returns an iterable of roots
def sqrtmodp(a, p, k=1):
  """
  find square roots of <a> mod <p>^<k> where <p> is a prime number.

  an iterator of roots is returned.
  """
  return (_sqrtmodp1(a, p, k) if a % p != 0 else _sqrtmodp2(a, p, k))

# find square roots of <a> mod <m> (i.e x such that pow(x, 2, m) = a (mod m)
# fs is (optionally) the prime factorisation of n
# returns an iterable of roots
def sqrtmod(a, m, fs=None):
  """
   find square roots of <a> mod <m>.

  i.e. values <x> such that (<x> * <x>) is congurent to <a> (mod <m>).

  an iterator of roots is returned.

  if you know the prime factorisation of <m> it can be passed as <fs>
  (which should be a multiset).

  if you know <m> is prime (or a power of a prime) you can use: sqrtmodp().
  """
  if m == 1: return [0]
  # solve for each prime power in the factorisation
  if not fs: fs = multiset.from_pairs(prime_factor(m))
  (roots, ms) = (list(), list())
  for (p, k) in fs.items():
    rs = sqrtmodp(a, p, k)
    if not rs: return []
    roots.append(rs)
    ms.append(p**k)
  if roots:
    # combine the results using CRT
    if len(roots) == 1: return roots[0]
    return (crt(vs, ms).x for vs in cproduct(roots))

######################################################################

# Cornacchia's algorithm:
# find all (X, Y) solutions for: X^2 + D.Y^2 = N [0 < D < N, gcd(D, N) = 1]

# generate primitive (= co-prime) solutions to: X^2 + D.Y^2 = N
# fs is (optionally) the prime factorisation of N
def cornacchia_primitive(D, N, fs=None):
  # deal with some simple cases
  if N == 1:
    yield (1, 0)
    if D == 1: yield (0, 1)
    return
  if D == N: yield (0, 1)
  r = is_square(N)
  if r is not None: yield (r, 0)
  if D >= N: return

  # now the general case ...
  assert 0 < D < N and gcd(D, N) == 1
  # look for square roots of -D (mod N)
  for z in sqrtmod(-D, N, fs):
    # we don't need roots > N/2
    if 2 * z > N: continue
    # euclidean descent
    (r, X) = (N, z)
    while X * X >= N:
      (r, X) = (X, r % X)
    # validate candidate solution
    Y = is_square(div(N - X * X, D))
    if Y is not None:
      yield (X, Y)
      if D == 1: yield (Y, X)

# find all solutions for X^2 + D.Y^2 = N
def cornacchia(D, N):
  # find prime factors of N
  fs = multiset.from_pairs(prime_factor(N))
  # collect primitive solutions
  ss = set(cornacchia_primitive(D, N, fs))
  # determine non-primitive solutions
  sqs = multiset.from_pairs((p, e // 2) for (p, e) in fs.items())
  for vs in sqs.subsets(min_size=1):
    g = multiply((v if k == 1 else v**k) for (v, k) in vs.items())
    n = N // sq(g)
    ss.update((g * X, g * Y) for (X, Y) in cornacchia_primitive(D, n, fs.difference(vs.multiply(2))))
  # return solutions in order
  return sorted(ss)

###############################################################################

# general LMM algorithm

# general PQa algorithm
def pells_PQa(P, Q, D, rD=None):
  if rD is None: rD = sqrtf(D)
  (A1, A2, B1, B2, G1, G2) = (1, 0, 0, 1, Q, -P)
  while 1:
    a = (P + rD) // Q
    (A1, A2, B1, B2, G1, G2) = (a * A1 + A2, A1, a * B1 + B2, B1, a * G1 + G2, G1)
    yield (P, Q, a, A1, B1, G1)
    P = a * Q - P
    Q = (D - P*P) // Q

# general LLM (for N^2 > D)
def pells_LMM(D, N, rD, u, v):
  (s, t) = pells1n_fundamental(D)
  for (f, _) in divisors_pairs(abs(N)):
    (M, r) = divmod(N, f * f)
    if r != 0: continue
    m = abs(M)
    for z in sqrtmod(D, m):
      if 2 * z > m: z = z - m
      # run the PQa algorithm until abs(Q) = 1, or (P, Q) repeats
      (pB, pG, seen) = (0, m, set())
      for (i, (P, Q, _, _, B, G)) in enumerate(pells_PQa(z, m, D, rD)):
        if abs(Q) == 1:
          if pG * pG - D * pB * pB == M:
            (x, y) = (f * pG, f * pB)
          else:
            assert pG * pG - D * pB * pB == -M
            (x, y) = (f * (pG * s + pB * D * t), f * (pB * s + pG * t))
          # find minimum non-negative solution for this family
          # (-x, -y) is in the same family as (x, y)
          if (N > 0 and x < 0) or (N < 0 and y < 0): (x, y) = (-x, -y)
          # advance until both x any y are non-negative
          while x < 0 or y < 0:
            (x, y) = (u * x + D * v * y, v * x + u * y)
          yield (x, y)
          break
        # have we hit the repeating section?
        if (P, Q) in seen: break
        seen.add((P, Q))
        (pB, pG) = (B, G)

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

# threshold at which we switch to Cornacchia's algorithm (or issue a warning)
cornacchia_threshold = 50000

# a.X^2 + b.Y^2 = c [a > 0, b > 0, c > 0]
def _diop_quad_bp(a, b, c):
  if c < 0: return

  # start by considering possible Y values
  Y0 = sqrtf(divf(c, b))
  Y1 = sqrtc(divc(a * c, b * (a + b)))
  # are there too many values to brute force?
  if Y0 - Y1 >= cornacchia_threshold:
    if a == 1 and b < c and gcd(b, c) == 1:
      if verbose: printf("[pells] switching to Cornacchia (instead of brute forcing {n} values)", n=Y0 - Y1)
      for s in cornacchia(b, c): yield s
      return
    # otherwise fall back to brute force

  # brute force search
  if verbose or Y0 - Y1 >= cornacchia_threshold: printf("[pells] brute forcing {n} values", n=Y0 - Y1)
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

# X^2 - (dY)^2 = N [d > 0, N != 0]
def _diop_quad_d2(d, N):
  if N > 0:
    pqs = rev(divisors_pairs(N))
    for (p, q) in pqs:
      (X, Y) = (div(p + q, 2), div(q - p, 2 * d))
      if X is not None and Y is not None:
        yield (X, Y)
  else:
    pqs = rev(divisors_pairs(-N))
    for (p, q) in pqs:
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

# threshold at which candidate solutions are considered exhausted
diop_quad_threshold = 10000

# results (X, Y) for increasing X
# [consider adding: fnX, fnY to validate (X, Y) solutions
def diop_quad(a, b, c, maxC=None, validate=0):
  """
  find non-negative integer solution to the equation: a.X^2 + b.Y^2 = c.

  an iterator of (X, Y) solutions are returned in increasing X order.
  """
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

  # find a multiplier m, such that (m * a) is a square
  m = multiply(p for (p, e) in prime_factor(a) if e % 2 == 1)
  r = is_square(m * a)
  if maxC is None: maxC = diop_quad_threshold
  def fn():
    for (i, (X, Y)) in enumerate(_diop_quad_a1(m * -b, m * c), start=1):
      #printf("[diop_quad: considering: X={X} Y={Y}]")
      X = div(X, r)
      if X is not None: yield (X, Y)
      if i == maxC:
        printf("diop_quad: WARNING: terminating search (after {maxC} candidates)")
        return
  return fn()

# filter the results of a call to diop_quad using function <f>
# return the first <n> results (default = all)
def diop_quadf(a, b, c, n=None, f=identity, maxC=None, validate=0):
  j = 0
  for (i, (X, Y)) in enumerate(diop_quad(a, b, c, maxC=maxC, validate=validate)):
    r = catch(f, X, Y)
    if r is None: continue
    yield r
    j += 1
    if j == n: break
    if i == maxC:
      printf("diop_quadf: WARNING: terminating search (after {maxC} candidates)")
      return

######################################################################

if enigma._namecheck(__name__):
  from enigma import (timer, arg, number as num)

  # c can be an expression
  (a, b, c, N, v) = (arg(-13, 0, num), arg(1, 1, num), arg(4, 2, eval), arg(20, 3, num), arg(0, 4, num))
  printf("solving: {a}.X^2 + {b}.Y^2 = {c}")
  if v: verbose = 1
  if verbose: timer.start()
  for (i, (X, Y)) in enumerate(diop_quad(a, b, c, validate=1), start=1):
    r = a * X * X + b * Y * Y
    printf("[{i}] X={X} Y={Y} -> {r}{v}", v=(" [FAIL!]" if r != c else ""))
    if i == N: printf("[first {N} solutions]"); break
  else:
    printf("[all solutions]")
  if verbose: timer.stop()
