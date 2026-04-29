"""
lithicscarpattern_py
====================
Python backend for the lithicscarpattern R package.

The R side imports this package and calls functions exposed in `api`.
End users should not import this directly; use the R package interface
(spharm_from_directions, spharm_from_meshes).
"""

from . import api  # noqa: F401
