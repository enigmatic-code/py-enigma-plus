#! python3
# -*- mode: Python; python-indent-offset: 2; coding: utf-8 -*-

# routines for dealing with (simple undirected) graphs

from __future__ import print_function

from enigma import (
  enigma, namedtuple, defaultdict, static, group, cproduct, subsets, is_disjoint, update, fail
)

__author__ = "Jim Randell <jim.randell@gmail.com>"
__version__ = "2025-10-29"

graph = enigma.module(__name__)

######################################################################

# edges -> adjacency matrix
def edges2adj(es, vs=()):
  adj = dict()
  for (x, y) in es:
    if x not in adj: adj[x] = set()
    adj[x].add(y)
    if y not in adj: adj[y] = set()
    adj[y].add(x)
  for v in vs:
    if v not in adj: adj[v] = set()
  return adj

# adjacency matrix -> edges
def adj2edges(adj):
  for x in sorted(adj.keys()):
    for y in adj[x]:
      if not (x > y):
        yield (x, y)

######################################################################

# isomoprhisms

# find an isomorphism for graph <adj> to a graph in <adjs>
# return (<index>, <map>) where <index> is an index into <adjs>
# and <map> maps nodes in <adj> to <adj1>
@static(rtype=None, fail=None)
def find_isomorphism(adj, adjs):
  if find_isomorphism.rtype is None:
    # set up return type
    find_isomorphism.rtype = namedtuple('ISO', 'adj map')
    find_isomorphism.fail = find_isomorphism.rtype(None, None)
  # group nodes by degree
  deg = lambda adj: group(adj.keys(), by=(lambda k: len(adj[k])))
  g0 = deg(adj)
  ks = sorted(g0.keys())
  # we only need to consider graphs with the same degree signature
  adjs_ = list()
  for (i, adj1) in enumerate(adjs):
    g1 = deg(adj1)
    if sorted(g1.keys()) == ks:
      adjs_.append((i, adj1, g1))
  if adjs_:
    xss = list(g0[k] for k in ks)
    # permute each degree separately
    for (i, adj1, g1) in adjs_:
      for yss in cproduct(subsets(g1[k], size=len, select='P') for k in ks):
        # construct the map from adj -> adj1
        m = dict()
        for (xs, ys) in zip(xss, yss):
          m.update(zip(xs, ys))
        # remap adj -> adj_
        adj_ = dict((m[k], set(m[v] for v in vs)) for (k, vs) in adj.items())
        if adj_ == adj1:
          return find_isomorphism.rtype(i, m)
  return find_isomorphism.fail

# check two graphs (adjacency matrix) are isomorphic
# return a map of nodes in adj1 to nodes in adj2, or None
def is_isomorphic(adj1, adj2): return find_isomorphism(adj1, [adj2]).map

######################################################################

# bipartite graphs

# subgraph of (x, y) edges that connect X and Y
def bipartite_subgraph(edges, X, Y):
  (X, Y) = (set(X), set(Y))
  return set((x, y) for (x, y) in edges if x in X and y in Y)

# matchings

# (k, vs) -> len(vs)
_key_fn = (lambda k_v: len(k_v[1]))

# remove s -> t from adj
_adj_remove = lambda adj, s, t: dict((k, vs.difference({t})) for (k, vs) in adj.items() if k != s)

def _matching(adj_xy, adj_yx, d):
  # are we done?
  if adj_xy and adj_yx:
    # find minimal x->y and y->x sets
    ((x, ys), (y, xs)) = (min(adj_xy.items(), key=_key_fn), min(adj_yx.items(), key=_key_fn))
    if not (xs and ys): return
    # process the smallest choice
    if len(xs) < len(ys):
      ys = {y}
    else:
      xs = {x}
    for (x, y) in cproduct([xs, ys]):
      adj_xy_ = _adj_remove(adj_xy, x, y)
      adj_yx_ = _adj_remove(adj_yx, y, x)
      #yield from _matching(adj_xy_, adj_yx_, update(d, [(x, y)]))  # [Python 3]
      for z in _matching(adj_xy_, adj_yx_, update(d, [(x, y)])): yield z  # [Python 2]
  elif not (adj_xy or adj_yx):
    yield d

# find (perfect) matchings in the bipartite graph specified by (x, y) edges
def find_bipartite_matching(edges, X=None, Y=None):
  # construct x -> y, y -> x adjacency matrices
  adj_xy = (dict((k, set()) for k in X) if X else defaultdict(set))
  adj_yx = (dict((k, set()) for k in Y) if Y else defaultdict(set))
  for (x, y) in edges:
    adj_xy[x].add(y)
    adj_yx[y].add(x)
  fail(not is_disjoint([adj_xy.keys(), adj_yx.keys()]), "find_bipartite_matching: graph is not bipartite")
  return _matching(adj_xy, adj_yx, dict())

######################################################################
