import math

import numpy as np
import pytest
from scipy.spatial import ConvexHull, HalfspaceIntersection

from solidwright.core import polyhedron as ph
from solidwright.core import vec
from solidwright.core.facets import cluster_planes

import shapes


def unit_box():
    return ph.box((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))


def from_trimesh(mesh):
    v, f = mesh.as_lists()
    return ph.Polyhedron(v, f)


def euler(poly):
    edges = set()
    for f in poly.faces:
        for k in range(len(f)):
            a, b = f[k], f[(k + 1) % len(f)]
            edges.add((min(a, b), max(a, b)))
    return len(poly.vertices) - len(edges) + len(poly.faces)


def assert_solid(poly):
    assert ph.closed_problems(poly) == []
    assert euler(poly) == 2
    assert poly.volume() > 0.0


# ---------------------------------------------------------------- polyhedron basics

def test_box_is_closed_and_points_outward():
    b = unit_box()
    assert_solid(b)
    assert abs(b.volume() - 1.0) < 1e-12 and abs(b.area() - 6.0) < 1e-12
    for k in range(6):
        n, d, flat = b.plane(k)
        centre = vec.mean(b.polygon(k))
        assert vec.dot(n, vec.sub(centre, (0.5, 0.5, 0.5))) > 0.0 and flat == 0.0


def test_triangles_keep_the_volume():
    b = unit_box()
    v, tris = b.triangles()
    assert len(tris) == 12
    total = sum(vec.dot(v[a], vec.cross(v[b_], v[c])) for a, b_, c in tris) / 6.0
    assert abs(total - 1.0) < 1e-12


# ---------------------------------------------------------------- clipping

def test_clip_in_half():
    cut, problem = ph.clip(unit_box(), (1.0, 0.0, 0.0), 0.5, label=7)
    assert problem is None
    assert_solid(cut)
    assert abs(cut.volume() - 0.5) < 1e-12
    k = cut.labels.index(7)
    n, d, _ = cut.plane(k)
    assert abs(n[0] - 1.0) < 1e-12 and abs(d - 0.5) < 1e-12
    assert len(cut.vertices) == 8 and len(cut.faces) == 6


def test_clip_that_misses_returns_the_same_solid():
    b = unit_box()
    same, problem = ph.clip(b, (1.0, 0.0, 0.0), 2.0)
    assert same is b and problem is None
    # Touching at a face is not a cut either.
    same, _ = ph.clip(b, (1.0, 0.0, 0.0), 1.0, eps=1e-12)
    assert same is b


def test_clip_that_removes_everything():
    gone, problem = ph.clip(unit_box(), (1.0, 0.0, 0.0), -0.5)
    assert gone is None and "whole solid" in problem


def test_clip_through_corners_uses_them():
    # The plane x = y runs through two vertical edges of the box: no new corners there.
    n = vec.unit((1.0, -1.0, 0.0))
    cut, _ = ph.clip(unit_box(), n, 0.0, eps=1e-12)
    assert_solid(cut)
    assert len(cut.vertices) == 6 and len(cut.faces) == 5
    assert abs(cut.volume() - 0.5) < 1e-12


def test_clip_keeps_corners_on_straight_runs_of_the_rim():
    # Every face of a triangulated box has a diagonal, and a cut at x = 0.5 crosses
    # four of them halfway: those corners sit on straight runs of the new face's rim.
    cut, _ = ph.clip(from_trimesh(shapes.box()), (1.0, 0.0, 0.0), 0.5, label="cap")
    assert_solid(cut)
    rim = cut.faces[cut.labels.index("cap")]
    assert len(rim) == 8
    assert abs(cut.volume() - 0.5) < 1e-12


def raised_edge_box(rise):
    """Unit box with an extra corner halfway along its top front edge, raised by rise."""
    b = unit_box()
    v = b.vertices + [(0.5, 0.0, 1.0 + rise)]
    faces = [list(f) for f in b.faces]
    faces[1] = [4, 8, 5, 6, 7]   # top
    faces[2] = [0, 1, 5, 8, 4]   # front
    return ph.Polyhedron(v, faces)


def test_corner_a_hair_outside_is_left_alone():
    # Only the raised corner is outside, by 5e-10 (within ISLAND_EPS x eps), and
    # its neighbours lie on the plane: there is nothing to cut off.
    tilt = vec.unit((0.0, -1e-3, 1.0))
    b = raised_edge_box(5e-10)
    same, problem = ph.clip(b, tilt, vec.dot(tilt, (0.0, 0.0, 1.0)), eps=1e-12)
    assert same is b and problem is None
    # Further out, such a cut would leave a face with no area: refused, not guessed.
    gone, problem = ph.clip(raised_edge_box(1e-8), tilt, vec.dot(tilt, (0.0, 0.0, 1.0)), eps=1e-12)
    assert gone is None and "no face to close it" in problem
    # Far out, it is cut off, which leaves the plain unit box.
    cut, problem = ph.clip(raised_edge_box(0.1), (0.0, 0.0, 1.0), 1.0, eps=1e-12)
    assert problem is None and len(cut.vertices) == 8 and abs(cut.volume() - 1.0) < 1e-12
    assert_solid(cut)


def test_stray_outside_corner_counts_as_on_the_plane():
    # The back top edge is clearly outside, the raised front corner only by 1e-10
    # (within ISLAND_EPS x eps) and not connected to it: it stays where it is.
    b = raised_edge_box(1e-10)
    n = vec.unit((0.0, 0.1, 1.0))
    cut, problem = ph.clip(b, n, vec.dot(n, (0.0, 0.0, 1.0)), eps=1e-12)
    assert problem is None
    assert_solid(cut)
    assert (0.5, 0.0, 1.0 + 1e-10) in cut.vertices
    # Further out than that, the two separate cuts are refused rather than guessed.
    gone, problem = ph.clip(raised_edge_box(1e-6), n, vec.dot(n, (0.0, 0.0, 1.0)), eps=1e-12)
    assert gone is None and "more than one place" in problem


def test_shared_edges_get_one_corner():
    cut, _ = ph.clip(unit_box(), vec.unit((1.0, 1.0, 1.0)), 0.9)
    used = [i for f in cut.faces for i in f]
    coords = set(cut.vertices)
    assert len(coords) == len(cut.vertices) == len(set(used))


@pytest.mark.parametrize("seed", range(6))
def test_many_planes_match_scipy(seed):
    g = shapes.rng(seed)
    normals = [tuple(x / np.linalg.norm(x)) for x in g.normal(size=(40, 3))]
    offsets = [float(o) for o in g.uniform(0.8, 1.2, size=40)]
    planes = list(zip(normals, offsets))
    poly, problem = ph.intersect_planes(planes, (-3.0, -3.0, -3.0), (3.0, 3.0, 3.0), 1e-12)
    assert problem is None
    assert_solid(poly)
    assert all(lab is not None for lab in poly.labels)
    for p in poly.vertices:
        assert all(vec.dot(n, p) <= d + 1e-9 for n, d in planes)
    halfspaces = np.array([list(n) + [-d] for n, d in planes])
    truth = ConvexHull(HalfspaceIntersection(halfspaces, np.zeros(3)).intersections).volume
    assert abs(poly.volume() - truth) < 1e-9 * truth


def test_unclosed_directions_leave_box_faces():
    planes = [((1.0, 0.0, 0.0), 1.0), ((-1.0, 0.0, 0.0), 1.0)]
    poly, _ = ph.intersect_planes(planes, (-5.0, -5.0, -5.0), (5.0, 5.0, 5.0))
    assert poly.labels.count(None) == 4


# ---------------------------------------------------------------- merge and checks

def test_merge_joins_a_tiny_cut_back_into_one_corner():
    # Shaving a corner by 1e-6 leaves three corners 1e-6 apart: one corner again after the merge.
    n = vec.unit((1.0, 1.0, 1.0))
    shaved, _ = ph.clip(unit_box(), n, vec.dot(n, (1.0, 1.0, 1.0)) - 1e-6)
    assert len(shaved.vertices) == 10 and len(shaved.faces) == 7
    merged = ph.merge_corners(shaved, 1e-4)
    assert_solid(merged)
    assert len(merged.vertices) == 8 and len(merged.faces) == 6
    problems, planes = ph.validate(merged, 1e-6)
    assert problems == [] and len(planes) == 6


def test_merge_is_transitive():
    # Corners 0.6 apart in a row: each is within 1 of its neighbour only, yet all join.
    poly = ph.Polyhedron([(0.0, 0, 0), (0.6, 0, 0), (1.2, 0, 0), (1.8, 0, 0), (0, 5, 0)],
                         [[0, 1, 2, 3, 4]])
    merged = ph.merge_corners(poly, 1.0)
    assert merged.faces == []  # the face collapsed to two corners
    assert ph.merge_corners(poly, 0.5).faces == [[0, 1, 2, 3, 4]]


def test_closed_problems_are_named():
    b = unit_box()
    assert ph.closed_problems(b) == []
    opened = ph.Polyhedron(b.vertices, b.faces[:-1])
    assert "4 open edges" in "; ".join(ph.closed_problems(opened))
    flipped = ph.Polyhedron(b.vertices, b.faces[:-1] + [list(reversed(b.faces[-1]))])
    text = "; ".join(ph.closed_problems(flipped))
    assert "run the same way twice" in text
    pinched = ph.Polyhedron(b.vertices, b.faces[:-1] + [b.faces[-1] + [b.faces[-1][0]]])
    assert "pinched" in "; ".join(ph.closed_problems(pinched))


def test_validate_flags_bent_faces_and_lost_clearance():
    b = unit_box()
    bent = ph.Polyhedron([(x, y, z + 0.01 if (x, y, z) == (1.0, 1.0, 1.0) else z)
                          for x, y, z in b.vertices], b.faces)
    problems, _ = ph.validate(bent, 1e-3)
    assert any("off flat" in p for p in problems)
    # Points 0.1 inside the faces, with a clearance of 0.2 asked for.
    pts = [(x, y, z) for x in (0.1, 0.9) for y in (0.1, 0.9) for z in (0.1, 0.9)]
    normals = [b.plane(k)[0] for k in range(6)]
    ext = ph.Extent(pts, normals)
    problems, _ = ph.validate(b, 1e-6, ext, clearance=0.2)
    assert any("6 faces closer to the part" in p for p in problems)
    assert ph.validate(b, 1e-6, ext, clearance=0.1)[0] == []
    # With a vertical axis only the four sides need the clearance.
    problems, _ = ph.validate(b, 1e-6, ext, clearance=0.2, axis=(0.0, 0.0, 1.0))
    assert any("4 faces closer" in p for p in problems)


def test_support_planes_touch_then_push_out():
    pts = [(x, y, z) for x in (0.0, 1.0) for y in (0.0, 1.0) for z in (0.0, 1.0)]
    slope = vec.unit((1.0, 0.0, 1.0))
    planes = ph.support_planes(pts, [(1.0, 0.0, 0.0), (0.0, 0.0, 1.0), slope], 0.1)
    assert [round(d, 12) for _, d in planes] == [1.1, 1.1, round(math.sqrt(2.0) + 0.1, 12)]
    planes = ph.support_planes(pts, [(1.0, 0.0, 0.0), (0.0, 0.0, 1.0), slope], 0.1, axis=(0.0, 0.0, 1.0))
    assert abs(planes[0][1] - 1.1) < 1e-12
    assert abs(planes[1][1] - 1.0) < 1e-12
    assert abs(planes[2][1] - (math.sqrt(2.0) + 0.1 * math.sqrt(0.5))) < 1e-12


def test_extent_reach_matches_brute_force():
    g = shapes.rng(3)
    pts = [tuple(p) for p in g.normal(size=(500, 3))]
    normals = [tuple(x / np.linalg.norm(x)) for x in g.normal(size=(20, 3))]
    ext = ph.Extent(pts, normals)
    for k, n in enumerate(normals):
        for turn in (1e-6, 1e-3, 0.3):
            m = vec.unit(vec.add(n, tuple(turn * x for x in g.normal(size=3))))
            assert ext.reach(m, k) == max(vec.dot(m, p) for p in pts)


# ---------------------------------------------------------------- offset_solid

def test_zero_clearance_rebuilds_a_convex_part_exactly():
    m = shapes.faceted_part()
    v, f = m.as_lists()
    res = ph.offset_solid(v, f)
    assert res.status == 0
    assert_solid(res.solid)
    assert abs(res.solid.volume() - m.volume()) < 1e-12
    assert len(res.solid.faces) == len(cluster_planes(v, f)) == len(res.planes)
    assert res.merged is None


def test_box_with_clearance_all_round():
    v, f = shapes.box().as_lists()
    res = ph.offset_solid(v, f, clearance=0.1)
    assert res.status == 0
    assert len(res.solid.vertices) == 8
    assert abs(res.solid.volume() - 1.2 ** 3) < 1e-12


def test_box_with_clearance_square_to_an_axis():
    v, f = shapes.box().as_lists()
    res = ph.offset_solid(v, f, clearance=0.1, axis=(0.0, 0.0, 2.0))
    assert res.status == 0
    assert abs(res.solid.volume() - 1.2 * 1.2 * 1.0) < 1e-12


def test_sideways_clearance_keeps_the_flat_end_and_clears_every_point():
    m = shapes.faceted_part()
    v, f = m.as_lists()
    clr = 0.05
    res = ph.offset_solid(v, f, clearance=clr, axis=(0.0, 0.0, 1.0), tolerance=0.001)
    assert res.status == 0
    assert_solid(res.solid)
    assert abs(max(p[2] for p in res.solid.vertices) - max(p[2] for p in v)) < 1e-9
    for n, d in res.planes:
        need = clr * ph.side_share(n, (0.0, 0.0, 1.0))
        gap = d - max(vec.dot(n, p) for p in v)
        assert gap >= ph.KEEP_CLEARANCE * need - 1e-12
    for k in range(len(res.solid.faces)):
        assert res.solid.plane(k)[2] <= 0.5 * 0.001
    assert min(res.solid.edge_lengths()) >= 0.001


def test_split_apex_is_merged():
    # The four sides of an off-centre pyramid slope differently, so a sideways
    # offset moves them by different amounts and the apex splits apart.
    v, f = shapes.pyramid().as_lists()
    res = ph.offset_solid(v, f, clearance=0.05, axis=(0.0, 0.0, 1.0), tolerance=0.001)
    assert res.status == 0
    assert len(res.solid.vertices) == 5 and len(res.solid.faces) == 5
    assert res.merged == pytest.approx(0.0125)
    assert "split corners merged" in str(res.report)
    planes = ph.support_planes(v, [fc.normal for fc in res.facets], 0.05, (0.0, 0.0, 1.0))
    exact, _ = ph.intersect_planes(planes, (-5.0, -5.0, -5.0), (5.0, 5.0, 5.0))
    assert len(exact.vertices) > 5  # the split the merge removed


def test_merge_that_breaks_the_solid_is_not_used():
    # A tolerance of 0.2 on a part 1 across: merging within it would collapse the
    # solid, so the exact unmerged solid is kept and its short edges reported.
    m = shapes.faceted_part()
    v, f = m.as_lists()
    res = ph.offset_solid(v, f, tolerance=0.2)
    assert res.status == 1
    text = str(res.report)
    assert "short edges kept: merging their corners within 0.2 would break the solid" in text
    assert "below the tolerance" in text
    assert_solid(res.solid)
    assert abs(res.solid.volume() - m.volume()) < 1e-12 and res.merged is None


def test_rotation_does_not_change_the_solid():
    m = shapes.faceted_part()
    v, f = m.as_lists()
    base = ph.offset_solid(v, f, clearance=0.03, axis=(0.0, 0.0, 1.0)).solid.volume()
    g = shapes.rng(5)
    placed, rot = shapes.randomly_placed(m, g)
    pv, pf = placed.as_lists()
    res = ph.offset_solid(pv, pf, clearance=0.03, axis=tuple(rot.dot((0.0, 0.0, 1.0))))
    assert res.status == 0
    assert abs(res.solid.volume() - base) < 1e-9 * base


def test_cut_along_a_chain_of_near_corners():
    # Found by a stress run: tilted and float32, the lower apex of this part
    # becomes a chain of corners within 2e-9 of a later plane, some just
    # outside and some on it, which used to split the cut into two rims.
    g = shapes.rng(91)
    shapes.noisy(shapes.subdivided(shapes.faceted_part(), 1), g, 2e-5)  # advances g as the run did
    placed, rot = shapes.randomly_placed(shapes.faceted_part(), g, spread=10.0)
    v, f = shapes.as_float32(placed).as_lists()
    res = ph.offset_solid(v, f, clearance=0.05, axis=tuple(rot.dot((0.0, 0.0, 1.0))), tolerance=0.001)
    assert res.status == 0, str(res.report)
    assert_solid(res.solid)


@pytest.mark.parametrize("seed", [0, 3, 7, 9, 14])
def test_placed_float32_parts(seed):
    # Rhino stores mesh vertices in single precision. Placed far from the origin
    # and tilted, that noise splits corners into clusters about 1e-7 across,
    # which a cut can pass right through; these seeds once broke the rim chaining
    # and the clearance check.
    g = shapes.rng(seed)
    for part in (shapes.faceted_part(), shapes.faceted_part(8), shapes.pyramid()):
        placed, rot = shapes.randomly_placed(part, g, spread=10.0)
        v, f = shapes.as_float32(placed).as_lists()
        axis = tuple(rot.dot((0.0, 0.0, 1.0)))
        for clr, a in ((0.05, axis), (0.05, None), (0.3, axis)):
            res = ph.offset_solid(v, f, clearance=clr, axis=a, tolerance=0.001)
            assert res.status == 0, str(res.report)
            assert_solid(res.solid)


def test_float32_input():
    v, f = shapes.as_float32(shapes.faceted_part()).as_lists()
    res = ph.offset_solid(v, f, clearance=0.02, axis=(0.0, 0.0, 1.0))
    assert res.report.ok
    assert_solid(res.solid)


def test_a_part_that_is_not_convex_gets_its_hollows_bridged():
    m = shapes.l_block()
    v, f = m.as_lists()
    res = ph.offset_solid(v, f)
    assert res.status == 2
    assert "not convex: 2 facet planes pass through it" in str(res.report)
    # Only the part's own facet directions are used: the L's bounding block, not its hull.
    assert abs(res.solid.volume() - 4.0) < 1e-12
    assert ConvexHull(m.vertices).volume < 4.0


def test_curved_part_is_refused():
    v, f = shapes.faceted_part(sides=200).as_lists()
    res = ph.offset_solid(v, f)
    assert res.status == 3 and res.solid is None
    assert "looks curved rather than faceted" in str(res.report)


def test_missing_flat_end_is_reported():
    m = shapes.faceted_part()
    v, f = m.as_lists()
    facets = [fc for fc in cluster_planes(v, f) if fc.normal[2] < 0.99]
    res = ph.offset_solid(v, f, facets=facets)
    assert res.status == 2 and "a facet may be missing" in str(res.report)


def test_open_direction_fails():
    m = shapes.box()
    v, f = m.as_lists()
    f = [t for t in f if not all(v[i][2] == 1.0 for i in t)]
    res = ph.offset_solid(v, f)
    assert res.status == 3 and res.solid is None
    assert "open in some direction" in str(res.report)


def test_noisy_facets_still_clear_every_point():
    g = shapes.rng(7)
    m = shapes.noisy(shapes.subdivided(shapes.faceted_part(), 2), g, 1e-5)
    v, f = m.as_lists()
    res = ph.offset_solid(v, f, clearance=0.05, axis=(0.0, 0.0, 1.0), tolerance=0.001)
    assert res.report.ok
    for n, d in res.planes:
        need = 0.05 * ph.side_share(n, (0.0, 0.0, 1.0))
        assert d - max(vec.dot(n, p) for p in v) >= ph.KEEP_CLEARANCE * need - 0.0005


@pytest.mark.parametrize("kwargs, text", [
    ({"clearance": -0.1}, "zero or more"),
    ({"tolerance": 0.0}, "more than zero"),
    ({"axis": (0.0, 0.0, 0.0)}, "zero vector"),
])
def test_bad_inputs(kwargs, text):
    v, f = shapes.box().as_lists()
    res = ph.offset_solid(v, f, **kwargs)
    assert res.status == 3 and text in str(res.report)


def test_no_faces_and_no_facets():
    assert "no faces" in str(ph.offset_solid([], []).report)
    v, f = shapes.box().as_lists()
    res = ph.offset_solid(v, f, facets=[])
    assert res.status == 3 and "no flat facets" in str(res.report)


# ---------------------------------------------------------------- faces with holes

def square_frame():
    """A 2 x 2 x 1 block with a 1 x 1 square hole through it, as a Polyhedron with holes."""
    v = []
    for z in (0.0, 1.0):
        v += [(0, 0, z), (2, 0, z), (2, 2, z), (0, 2, z)]                  # outer: 0..3 and 8..11
        v += [(0.5, 0.5, z), (1.5, 0.5, z), (1.5, 1.5, z), (0.5, 1.5, z)]  # inner: 4..7 and 12..15
    faces = [[8, 9, 10, 11], [3, 2, 1, 0]]
    holes = [[[15, 14, 13, 12]], [[4, 5, 6, 7]]]  # top and bottom, holes clockwise seen from outside
    for k in range(4):  # outer walls
        a, b = k, (k + 1) % 4
        faces.append([a, b, b + 8, a + 8])
        holes.append([])
    for k in range(4):  # walls of the hole, facing into it
        a, b = 4 + k, 4 + (k + 1) % 4
        faces.append([b, a, a + 8, b + 8])
        holes.append([])
    return ph.Polyhedron(v, faces, holes=holes)


def test_faces_with_holes_measure_right():
    p = square_frame()
    assert ph.closed_problems(p) == []
    assert abs(p.volume() - 3.0) < 1e-12
    assert abs(p.area() - (2 * 3.0 + 4 * 2.0 + 4 * 1.0)) < 1e-12
    n, d, flat = p.plane(0)
    assert n == (0.0, 0.0, 1.0) and d == 1.0 and flat == 0.0
    assert len(p.edge_lengths()) == 24


def test_triangles_of_faces_with_holes_close_the_solid():
    from solidwright.core.topology import EdgeCensus, signed_volume
    v, tris = square_frame().triangles()
    assert EdgeCensus(tris).closed
    assert abs(signed_volume(v, tris) - 3.0) < 1e-12


def test_merge_drops_a_collapsed_hole():
    # Shrink the hole to 2e-5 across: merging within 1e-4 removes it, walls and all.
    p = square_frame()
    v = [(1 + (x - 1) * 2e-5, 1 + (y - 1) * 2e-5, z) if 0.4 < x < 1.6 and 0.4 < y < 1.6 else (x, y, z)
         for x, y, z in p.vertices]
    merged = ph.merge_corners(ph.Polyhedron(v, p.faces, holes=p.holes), 1e-4)
    assert merged.holes[0] == [] and merged.holes[1] == []
    assert len(merged.faces) == 6
    assert ph.closed_problems(merged) == []


def test_clip_refuses_faces_with_holes():
    with pytest.raises(ValueError):
        ph.clip(square_frame(), (1.0, 0.0, 0.0), 1.0)


# ---------------------------------------------------------------- coplanar faces

def test_merge_coplanar_joins_a_triangulated_box_into_six_faces():
    tri_box = from_trimesh(shapes.box((2.0, 3.0, 1.0)))
    assert len(tri_box.faces) == 12
    merged = ph.merge_coplanar(tri_box, 1e-12)
    assert len(merged.faces) == 6 and all(len(f) == 4 for f in merged.faces)
    assert ph.closed_problems(merged) == []
    assert abs(merged.volume() - 6.0) < 1e-12


def test_merge_coplanar_makes_a_face_with_a_hole():
    # The top of a frame given as four trapezoids around the hole.
    p = square_frame()
    v, faces, holes = p.vertices, p.faces[1:], p.holes[1:]
    outer, inner = [8, 9, 10, 11], [12, 13, 14, 15]
    for k in range(4):
        a, b = outer[k], outer[(k + 1) % 4]
        c, d = inner[(k + 1) % 4], inner[k]
        faces = faces + [[a, b, c, d]]
        holes = holes + [[]]
    split = ph.Polyhedron(v, faces, holes=holes)
    assert ph.closed_problems(split) == []
    merged = ph.merge_coplanar(split, 1e-12)
    assert len(merged.faces) == 10
    top = [k for k in range(len(merged.faces)) if merged.plane(k)[0][2] > 0.999]
    assert len(top) == 1 and len(merged.holes[top[0]]) == 1
    assert abs(merged.volume() - 3.0) < 1e-12


def test_merge_coplanar_leaves_faces_at_an_angle_alone():
    b = unit_box()
    assert ph.merge_coplanar(b, 1e-12) is b
