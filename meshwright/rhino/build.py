"""Polyhedron faces to a closed Rhino Brep.

First try: one planar Brep per face, joined. That gives the cleanest result
(one Brep face per polyhedron face, exact edges). When Rhino will not make or
join the planar faces, the same faces go through a welded mesh instead, which
cannot leave gaps, and its coplanar triangles are merged back into flat faces.
"""

from __future__ import absolute_import, division

import Rhino.Geometry as rg

from ..core.report import Report


def polyhedron_to_brep(poly, tolerance):
    """(brep, report) for a meshwright.core.polyhedron.Polyhedron.

    The Brep is closed and points outward, or None when Rhino could not make
    a closed solid from the faces (the report says why).
    """
    report = Report("brep")
    brep = _from_planar_faces(poly, tolerance)
    if brep is None:
        brep = _from_mesh(poly, tolerance, report)
        if brep is None:
            return None, report
        report.warn("the planar faces did not join; built through a welded mesh instead")
    if brep.SolidOrientation == rg.BrepSolidOrientation.Inward:
        brep.Flip()
    if not brep.IsValid:
        report.messy("Rhino reports the Brep as not valid")
    report.note("%d Brep faces" % brep.Faces.Count)
    return brep, report


def _from_planar_faces(poly, tolerance):
    """One planar Brep per face, joined; None unless that gives one closed solid."""
    pieces = []
    for polygon in poly.polygons():
        points = [rg.Point3d(p[0], p[1], p[2]) for p in polygon]
        points.append(points[0])
        made = rg.Brep.CreatePlanarBreps(rg.PolylineCurve(points), tolerance)
        if not made or len(made) != 1:
            return None
        pieces.append(made[0])
    joined = rg.Brep.JoinBreps(pieces, tolerance)
    if not joined or len(joined) != 1 or not joined[0].IsSolid:
        return None
    return joined[0]


def _from_mesh(poly, tolerance, report):
    """The faces as one welded mesh, turned into a Brep with coplanar faces merged."""
    vertices, triangles = poly.triangles()
    mesh = rg.Mesh()
    mesh.Vertices.UseDoublePrecisionVertices = True
    for x, y, z in vertices:
        mesh.Vertices.Add(x, y, z)
    for a, b, c in triangles:
        mesh.Faces.AddFace(a, b, c)
    mesh.Normals.ComputeNormals()
    mesh.Compact()
    brep = rg.Brep.CreateFromMesh(mesh, True)
    if brep is None or not brep.IsSolid:
        report.fail("Rhino could not make a closed solid from the %d faces" % len(poly.faces))
        return None
    # The return value only says whether anything merged, so count the faces instead.
    brep.MergeCoplanarFaces(tolerance)
    if brep.Faces.Count > len(poly.faces):
        report.messy("coplanar triangles did not all merge: %d Brep faces for %d flat faces" % (
            brep.Faces.Count, len(poly.faces)))
    return brep
