#! python3
# -*- mode: Python; python-indent-offset: 2; coding: utf-8 -*-

from __future__ import print_function

from enigma import (namedtuple, static, group, cproduct, subsets)

# routines for dealing with (simple undirected) graphs

__author__ = "Jim Randell <jim.randell@gmail.com>"
__version__ = "2024-12-31"

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
