"""Rhino mesh to plain lists and back, plus document tolerance.

The only module that knows how RhinoCommon stores meshes. Everything it
returns is plain Python, ready for meshwright.core.
"""

from __future__ import absolute_import, division

import Rhino.Geometry as rg


def mesh_to_lists(mesh):
    """(vertices, triangles) from a Rhino mesh.

    Vertices come from the mesh's topology vertices, so coincident vertices
    (split for texture or normal seams) arrive welded. Topology vertices are
    single precision, so their positions are read from the mesh's own
    double-precision vertices when it has them (Rhino 8 meshes usually do).
    Quads are split along their shorter diagonal, which keeps the two
    triangles as close to the quad's real shape as possible.
    """
    top = mesh.TopologyVertices
    exact = getattr(mesh.Vertices, "UseDoublePrecisionVertices", False)
    vertices = []
    for i in range(top.Count):
        p = mesh.Vertices.Point3dAt(top.MeshVertexIndices(i)[0]) if exact else top[i]
        vertices.append((float(p.X), float(p.Y), float(p.Z)))
    faces = []
    for fi in range(mesh.Faces.Count):
        f = mesh.Faces[fi]
        a = top.TopologyVertexIndex(f.A)
        b = top.TopologyVertexIndex(f.B)
        c = top.TopologyVertexIndex(f.C)
        if f.IsQuad:
            d = top.TopologyVertexIndex(f.D)
            if _dist_sq(vertices[a], vertices[c]) <= _dist_sq(vertices[b], vertices[d]):
                faces.append((a, b, c))
                faces.append((a, c, d))
            else:
                faces.append((a, b, d))
                faces.append((b, c, d))
        else:
            faces.append((a, b, c))
    return vertices, faces


def lists_to_mesh(vertices, faces):
    """A Rhino mesh (double-precision vertices) from vertex tuples and triangle index triples."""
    mesh = rg.Mesh()
    mesh.Vertices.UseDoublePrecisionVertices = True
    for x, y, z in vertices:
        mesh.Vertices.Add(x, y, z)
    for a, b, c in faces:
        mesh.Faces.AddFace(a, b, c)
    mesh.Normals.ComputeNormals()
    mesh.Compact()
    return mesh


def to_point(p):
    return rg.Point3d(p[0], p[1], p[2])


def to_polyline(points, closed=False):
    """A Rhino polyline through (x, y, z) points, closed back to the first when asked."""
    line = rg.Polyline([to_point(p) for p in points])
    if closed and points:
        line.Add(to_point(points[0]))
    return line


def to_vector(v):
    return rg.Vector3d(v[0], v[1], v[2])


def from_vector(v):
    """(x, y, z) from a Vector3d, a Plane (its Z axis) or a 3-sequence; None if none of those."""
    if v is None:
        return None
    # Not hasattr(v, "ZAxis"): under pythonnet a Vector3d also answers that,
    # with the static Vector3d.ZAxis (world Z).
    if isinstance(v, rg.Plane):
        v = v.ZAxis
    if hasattr(v, "X"):
        return (float(v.X), float(v.Y), float(v.Z))
    try:
        return (float(v[0]), float(v[1]), float(v[2]))
    except (TypeError, IndexError, ValueError):
        return None


def to_plane(origin, frame):
    """A Rhino plane at origin with axes frame = (x, y, z)."""
    x, y, _ = frame
    return rg.Plane(rg.Point3d(origin[0], origin[1], origin[2]), to_vector(x), to_vector(y))


def document_tolerance(default=0.001):
    """The active document's absolute tolerance, or default outside Rhino."""
    try:
        import Rhino
        doc = Rhino.RhinoDoc.ActiveDoc
        if doc is not None:
            return doc.ModelAbsoluteTolerance
    except Exception:
        pass
    return default


def _dist_sq(a, b):
    dx, dy, dz = a[0] - b[0], a[1] - b[1], a[2] - b[2]
    return dx * dx + dy * dy + dz * dz
