#! python3
# -*- mode: Python; python-indent-offset: 2; coding: utf-8 -*-

from __future__ import print_function

from enigma import (group, cproduct, subsets)

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

# check two graphs (adjacency matrix) are isomorphic
# return a map of nodes in adj1 to nodes in adj2, or None
def is_isomorphic(adj1, adj2):
  # group nodes by degree
  (g1, g2) = (group(adj.keys(), by=(lambda k: len(adj[k]))) for adj in (adj1, adj2))
  # check there are the same number of nodes by degree
  ks = sorted(g1.keys())
  if ks != sorted(g2.keys()): return
  if any(len(g1[k]) != len(g2[k]) for k in ks): return

  # permute each degree separately
  xss = list(g1[k] for k in ks)
  for yss in cproduct(subsets(g2[k], size=len, select='P') for k in ks):
    # construct the map between graphs (adj1 -> adj2)
    m = dict()
    for (xs, ys) in zip(xss, yss):
      m.update(zip(xs, ys))
    # remap adj1 -> adj3 using the map
    adj3 = dict((m[k], set(m[v] for v in vs)) for (k, vs) in adj1.items())
    # is it the same as adj2?
    if adj3 == adj2: return m
