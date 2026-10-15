"""Rhino mesh to plain lists and back, plus document tolerance.

The only module that knows how RhinoCommon stores meshes. Everything it
returns is plain Python, ready for meshwright.core.
"""

from __future__ import absolute_import, division

import Rhino.Geometry as rg


def mesh_to_lists(mesh):
    """(vertices, triangles) from a Rhino mesh.

    Vertices come from the mesh's topology vertices, so coincident vertices
    (split for texture or normal seams) arrive welded. Quads are split along
    their shorter diagonal, which keeps the two triangles as close to the
    quad's real shape as possible.
    """
    top = mesh.TopologyVertices
    vertices = []
    for i in range(top.Count):
        p = top[i]
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
    """A Rhino mesh from vertex tuples and triangle index triples."""
    mesh = rg.Mesh()
    for x, y, z in vertices:
        mesh.Vertices.Add(x, y, z)
    for a, b, c in faces:
        mesh.Faces.AddFace(a, b, c)
    mesh.Normals.ComputeNormals()
    mesh.Compact()
    return mesh


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
