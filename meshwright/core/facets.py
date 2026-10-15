"""Coplanar clustering: group a mesh's triangles into flat facets.

A facet is a set of triangles lying on one plane (within an angle and an
offset tolerance). Each facet gets a plane refitted through all its
vertices, so a facet made of many slightly noisy triangles still gets one
clean plane.

How it works: triangles are taken largest first. Each one either joins an
existing facet whose seed plane it matches, or seeds a new facet. Candidate
facets are found through a hash of the seed normal, so the search stays fast
on large meshes. Membership is always tested against the facet's SEED plane
(not a drifting running average), and against all three corners of the
triangle (not just its centre), so a facet can never creep around a gently
curved surface.

The hash is sign agnostic: a triangle whose winding is reversed still finds
its facet. The facet's normal direction is then taken from the area-weighted
windings of its triangles.
"""

from __future__ import absolute_import, division

import math

from . import vec
from .topology import bounding_size, face_area, face_normal

ANGLE_DEG = 1.0        # default: triangles within 1 degree of the seed normal
OFFSET_RATIO = 0.003   # default: corners within this x mesh size of the seed plane


class Facet(object):
    """One flat region of a mesh.

    faces    indices of the triangles in it
    area     total area
    normal   unit normal, pointing the way most of its area faces
    offset   plane offset: points p on the plane satisfy dot(normal, p) = offset
    centre   area-weighted centre of its triangles
    rms      root-mean-square distance of its vertices from the refitted plane
    """

    def __init__(self, faces, area, normal, offset, centre, rms):
        self.faces = faces
        self.area = area
        self.normal = normal
        self.offset = offset
        self.centre = centre
        self.rms = rms

    def distance(self, point):
        """Signed distance of a point from the facet's plane."""
        return vec.dot(self.normal, point) - self.offset

    def __repr__(self):
        return "Facet(%d faces, area=%.4g, normal=(%.3f, %.3f, %.3f))" % (
            len(self.faces), self.area, self.normal[0], self.normal[1], self.normal[2])


def cluster_planes(vertices, faces, angle_deg=ANGLE_DEG, offset_tol=None):
    """Group triangles into flat facets. Returns facets sorted by area, largest first.

    angle_deg   how far a triangle's normal may turn from the facet's seed normal
    offset_tol  how far a triangle's corners may sit from the seed plane
                (default: OFFSET_RATIO x the mesh's bounding size)
    Zero-area triangles are skipped.
    """
    size = bounding_size(vertices)
    if offset_tol is None:
        offset_tol = OFFSET_RATIO * size
    cos_tol = math.cos(math.radians(angle_deg))
    cell = 2.0 * math.sin(math.radians(angle_deg) / 2.0)  # chord length of the angle

    order = []
    for i, f in enumerate(faces):
        n = face_normal(vertices, f)
        if n is not None:
            order.append((face_area(vertices, f), i, n))
    order.sort(key=lambda t: -t[0])

    seeds = []    # per facet: (seed normal, seed offset)
    members = []  # per facet: list of face indices
    grid = {}     # cell of the sign-normalised seed normal -> facet ids seeded there
    cache = {}    # cell -> nearby facet ids, valid until the grid changes
    for area, i, n in order:
        corners = [vertices[k] for k in faces[i]]
        joined = -1
        m, near_flip = _normalised(n, cell)
        nearby = _nearby(grid, cache, _cell_key(m, cell))
        if near_flip:  # close to where the sign rule flips: look on the other side too
            nearby = nearby + _nearby(grid, cache, _cell_key(vec.scale(m, -1.0), cell))
        for candidate in nearby:
            sn, so = seeds[candidate]
            c = vec.dot(n, sn)
            if abs(c) < cos_tol:
                continue
            if all(abs(vec.dot(sn, p) - so) <= offset_tol for p in corners):
                joined = candidate
                break
        if joined < 0:
            seed_offset = vec.dot(n, vec.mean(corners))
            joined = len(seeds)
            seeds.append((n, seed_offset))
            members.append([])
            grid.setdefault(_cell_key(m, cell), []).append(joined)
            cache.clear()
        members[joined].append(i)

    facets = [_fit(vertices, faces, group) for group in members]
    facets.sort(key=lambda fc: -fc.area)
    return facets


def _cell_key(n, cell):
    return (int(math.floor(n[0] / cell)), int(math.floor(n[1] / cell)), int(math.floor(n[2] / cell)))


def _normalised(n, cell):
    """n or -n, whichever has its largest component positive, so a normal and
    its reverse hash to the same cell. Also says when n is close enough to the
    flip-over (two components of nearly equal size) that both sides need a look."""
    a = [abs(n[0]), abs(n[1]), abs(n[2])]
    k = a.index(max(a))
    m = n if n[k] >= 0.0 else vec.scale(n, -1.0)
    second = sorted(a)[1]
    return m, a[k] - second <= 2.0 * cell


def _nearby(grid, cache, key):
    found = cache.get(key)
    if found is None:
        kx, ky, kz = key
        found = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    found.extend(grid.get((kx + dx, ky + dy, kz + dz), ()))
        cache[key] = found
    return found


def _fit(vertices, faces, group):
    """Refit one facet: area-weighted direction, least-squares plane through its corners."""
    area = 0.0
    wound = (0.0, 0.0, 0.0)
    cx = cy = cz = 0.0
    corner_ids = set()
    for i in group:
        f = faces[i]
        a = face_area(vertices, f)
        area += a
        n = face_normal(vertices, f)
        wound = vec.add(wound, vec.scale(n, a))
        c = vec.mean([vertices[k] for k in f])
        cx += c[0] * a
        cy += c[1] * a
        cz += c[2] * a
        corner_ids.update(f)
    centre = (cx / area, cy / area, cz / area)
    points = [vertices[k] for k in sorted(corner_ids)]
    if len(points) >= 3:
        _, normal, rms, _ = vec.plane_fit(points)
    else:
        normal, rms = vec.unit(wound), 0.0
    if vec.dot(normal, wound) < 0.0:
        normal = vec.scale(normal, -1.0)
    return Facet(sorted(group), area, normal, vec.dot(normal, centre), centre, rms)
