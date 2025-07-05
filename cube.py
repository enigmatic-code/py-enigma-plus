#! python3
# -*- mode: Python; py-indent-offset: 2; -*-

# the 24 rotations of the cube

from __future__ import print_function

__author__ = "Jim Randell <jim.randell@gmail.com>"
__version__ = "2025-07-02"

#from enigma import (enigma)
#cube = enigma.module(__name__)

###############################################################################

# labels for the faces
face_label = "UDLRFB"
(U, D, L, R, F, B) = (0, 1, 2, 3, 4, 5)

# the 24 rotational positions of the cube
# we record the positions of the faces
# and their orientations (in clockwise quarter turns)
_rotations = (
  # the first 6 match up with a rotation of the face U, D, L, R, F, B
  # faces               orientations
  # U  D  L  R  F  B    U  D  L  R  F  B      U F
  ((U, D, F, B, R, L), (1, 3, 0, 0, 0, 0)), # U R = rotate U
  ((U, D, B, F, L, R), (3, 1, 0, 0, 0, 0)), # U L = rotate D = rotate UUU
  ((B, F, L, R, U, D), (2, 0, 1, 3, 0, 2)), # B U = rotate L
  ((F, B, L, R, D, U), (0, 2, 3, 1, 0, 2)), # F D = rotate R = rotate LLL
  ((L, R, D, U, F, B), (1, 1, 1, 1, 1, 3)), # L F = rotate F
  ((R, L, U, D, F, B), (3, 3, 3, 3, 3, 1)), # R F = rotate B = rotate FFF
  # the remaining rotations can be derived from those above
  ((U, D, L, R, F, B), (0, 0, 0, 0, 0, 0)), # U F = rotate []
  ((U, D, R, L, B, F), (2, 2, 0, 0, 0, 0)), # U B = rotate UU
  ((D, U, R, L, F, B), (2, 2, 2, 2, 2, 2)), # D F = rotate FF
  ((D, U, L, R, B, F), (0, 0, 2, 2, 2, 2)), # D B = rotate LL
  ((D, U, F, B, L, R), (3, 1, 2, 2, 2, 2)), # D L = rotate ULL
  ((D, U, B, F, R, L), (1, 3, 2, 2, 2, 2)), # D R = rotate LLU
  ((L, R, F, B, U, D), (2, 0, 1, 3, 1, 1)), # L U = rotate FU
  ((L, R, B, F, D, U), (0, 2, 3, 1, 1, 1)), # L D = rotate FD
  ((L, R, U, D, B, F), (3, 3, 1, 1, 3, 1)), # L B = rotate FUU
  ((R, L, B, F, U, D), (2, 0, 1, 3, 3, 3)), # R U = rotate BD
  ((R, L, F, B, D, U), (0, 2, 3, 1, 3, 3)), # R D = rotate BU
  ((R, L, D, U, B, F), (1, 1, 3, 3, 1, 3)), # R B = rotate BUU
  ((F, B, R, L, U, D), (2, 0, 1, 3, 2, 0)), # F U = rotate RUU
  ((F, B, U, D, L, R), (3, 3, 2, 0, 3, 1)), # F L = rotate RD
  ((F, B, D, U, R, L), (1, 1, 0, 2, 1, 3)), # F R = rotate RU
  ((B, F, R, L, D, U), (0, 2, 3, 1, 2, 0)), # B D = rotate LUU
  ((B, F, D, U, L, R), (1, 1, 2, 0, 1, 3)), # B L = rotate LD
  ((B, F, U, D, R, L), (3, 3, 0, 2, 3, 1)), # B R = rotate LU
)

# map names to rotations U, D, L, R, F, B
_names = dict(zip(face_label, (U, D, L, R, F, B)))


# a class representing the rotations of the cube
class Cube(object):

  def __init__(self, faces=(U, D, L, R, F, B), orientations=(0, 0, 0, 0, 0, 0)):
    self.faces = tuple(faces)
    self.orientations = tuple(orientations)

  # a new cube derived from the old one by applying the specified transformation
  def transform(self, faces, orientations):
    return Cube(
      faces=tuple(self.faces[i] for i in faces),
      orientations=tuple((self.orientations[i] + r) % 4 for (i, r) in zip(faces, orientations)),
    )

  # generate all rotations of the cube
  def rotations(self):
    for (faces, orientations) in _rotations:
      yield self.transform(faces, orientations)

  # apply specific rotations
  def rotate(self, ts):
    cube = self
    for t in ts:
      cube = cube.transform(*(_rotations[_names.get(t, t)]))
    return cube

  # make a copy of this cube, with the specified face/orientation updates
  def update(self, faces=None, orientations=None):
    (fs, os) = (self.faces, self.orientations)
    if faces:
      fs = list(fs)
      for (k, v) in faces: fs[k] = v
    if orientations:
      os = list(os)
      for (k, v) in orientations: os[k] = v
    return Cube(faces=fs, orientations=os)
