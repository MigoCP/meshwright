"""Small 3D vector toolkit on plain tuples.

Pure Python on purpose: it runs unchanged in GHPython (IronPython 2.7), in
Rhino 8 CPython and in plain CPython for testing. Vectors and points are
3-tuples of floats. Functions never mutate their inputs.

Also here: a symmetric 3x3 eigen solver and the fits built on it (plane fit,
direction of least variance), because several tools need them.
"""

from __future__ import absolute_import, division

import math

TINY = 1e-12


# ---------------------------------------------------------------- basics

def add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def scale(a, s):
    return (a[0] * s, a[1] * s, a[2] * s)


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def length_sq(a):
    return a[0] * a[0] + a[1] * a[1] + a[2] * a[2]


def length(a):
    return math.sqrt(length_sq(a))


def distance(a, b):
    return length(sub(a, b))


def unit(a):
    """a scaled to length 1, or None when a is (nearly) zero."""
    n = length(a)
    if n <= TINY:
        return None
    return (a[0] / n, a[1] / n, a[2] / n)


def lerp(a, b, t):
    return (a[0] + (b[0] - a[0]) * t,
            a[1] + (b[1] - a[1]) * t,
            a[2] + (b[2] - a[2]) * t)


def mean(points):
    """Average of a non-empty sequence of points."""
    n = 0
    sx = sy = sz = 0.0
    for p in points:
        sx += p[0]
        sy += p[1]
        sz += p[2]
        n += 1
    if n == 0:
        raise ValueError("mean of no points")
    return (sx / n, sy / n, sz / n)


def angle_between(a, b):
    """Angle in radians between two non-zero vectors, in [0, pi]."""
    ua, ub = unit(a), unit(b)
    if ua is None or ub is None:
        raise ValueError("angle with a zero vector")
    c = max(-1.0, min(1.0, dot(ua, ub)))
    return math.acos(c)


def rotate(v, axis, angle):
    """Rotate v about the (non-zero) axis through the origin by angle radians."""
    k = unit(axis)
    if k is None:
        raise ValueError("rotation about a zero axis")
    c, s = math.cos(angle), math.sin(angle)
    kv = cross(k, v)
    kd = dot(k, v) * (1.0 - c)
    return (v[0] * c + kv[0] * s + k[0] * kd,
            v[1] * c + kv[1] * s + k[1] * kd,
            v[2] * c + kv[2] * s + k[2] * kd)


def area_vector(polygon):
    """Newell's method: normal of a closed polygon scaled by its area.

    Works for any simple planar polygon (convex or not) and is robust for
    slightly non-planar ones. Vertex order sets the direction.
    """
    nx = ny = nz = 0.0
    count = len(polygon)
    for i in range(count):
        x0, y0, z0 = polygon[i]
        x1, y1, z1 = polygon[(i + 1) % count]
        nx += (y0 - y1) * (z0 + z1)
        ny += (z0 - z1) * (x0 + x1)
        nz += (x0 - x1) * (y0 + y1)
    return (nx * 0.5, ny * 0.5, nz * 0.5)


# ---------------------------------------------------------------- frames

def perpendicular_pair(n):
    """Two unit vectors (u, v) so that (u, v, n) is a right-handed frame."""
    z = unit(n)
    if z is None:
        raise ValueError("frame around a zero vector")
    # Start from the world axis least aligned with z: always well conditioned.
    ax = [abs(z[0]), abs(z[1]), abs(z[2])]
    i = ax.index(min(ax))
    seed = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)][i]
    u = unit(cross(seed, z))
    v = cross(z, u)
    return u, v


def frame_from_axis(axis, x_hint=(1.0, 0.0, 0.0)):
    """Right-handed unit frame (x, y, z) with z along axis.

    x is the hint projected square to the axis, so frames stay predictable
    (for a vertical axis, x follows world X). When the hint is parallel to the
    axis, any perpendicular x is used.
    """
    z = unit(axis)
    if z is None:
        raise ValueError("frame from a zero axis")
    x = unit(sub(x_hint, scale(z, dot(x_hint, z))))
    if x is None:
        x, _ = perpendicular_pair(z)
    y = cross(z, x)
    return x, y, z


def angular_order(points, normal=None):
    """Indices of points sorted counter-clockwise around their centre.

    Counter-clockwise as seen looking down the normal (from its tip). When no
    normal is given, the best-fit plane normal is used. Useful for ordering
    unordered outline points.
    """
    if len(points) < 3:
        return list(range(len(points)))
    c = mean(points)
    if normal is None:
        normal = plane_fit(points)[1]
    u, v = perpendicular_pair(normal)
    keyed = []
    for i, p in enumerate(points):
        d = sub(p, c)
        keyed.append((math.atan2(dot(d, v), dot(d, u)), i))
    keyed.sort()
    return [i for _, i in keyed]


# ---------------------------------------------------------------- 3x3 symmetric eigen

def eig_sym3(m):
    """Eigen decomposition of a symmetric 3x3 matrix (cyclic Jacobi).

    Returns (values, vectors): values ascending, vectors[i] the unit
    eigenvector of values[i]. Only the upper triangle of m is read.
    """
    a = [[float(m[0][0]), float(m[0][1]), float(m[0][2])],
         [float(m[0][1]), float(m[1][1]), float(m[1][2])],
         [float(m[0][2]), float(m[1][2]), float(m[2][2])]]
    v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    for _ in range(60):
        off = a[0][1] ** 2 + a[0][2] ** 2 + a[1][2] ** 2
        diag = a[0][0] ** 2 + a[1][1] ** 2 + a[2][2] ** 2
        if off <= 1e-30 * diag or off == 0.0:
            break
        for p, q in ((0, 1), (0, 2), (1, 2)):
            apq = a[p][q]
            if apq == 0.0:
                continue
            theta = (a[q][q] - a[p][p]) / (2.0 * apq)
            if abs(theta) > 1e150:
                t = 0.5 / theta
            else:
                t = 1.0 / (abs(theta) + math.sqrt(theta * theta + 1.0))
                if theta < 0.0:
                    t = -t
            c = 1.0 / math.sqrt(t * t + 1.0)
            s = t * c
            for k in range(3):  # columns p, q of A
                akp, akq = a[k][p], a[k][q]
                a[k][p] = c * akp - s * akq
                a[k][q] = s * akp + c * akq
            for k in range(3):  # rows p, q of A
                apk, aqk = a[p][k], a[q][k]
                a[p][k] = c * apk - s * aqk
                a[q][k] = s * apk + c * aqk
            for k in range(3):  # accumulate the rotation
                vkp, vkq = v[k][p], v[k][q]
                v[k][p] = c * vkp - s * vkq
                v[k][q] = s * vkp + c * vkq
    pairs = []
    for i in range(3):
        vec = unit((v[0][i], v[1][i], v[2][i]))
        pairs.append((a[i][i], vec))
    pairs.sort(key=lambda pr: pr[0])
    return [pr[0] for pr in pairs], [pr[1] for pr in pairs]


# ---------------------------------------------------------------- fits

def covariance(points, weights=None):
    """Weighted centre and 3x3 covariance (divided by the total weight)."""
    if weights is None:
        weights = [1.0] * len(points)
    if len(points) == 0 or len(points) != len(weights):
        raise ValueError("covariance needs points with one weight each")
    total = float(sum(weights))
    if total <= 0.0:
        raise ValueError("covariance needs a positive total weight")
    cx = cy = cz = 0.0
    for p, w in zip(points, weights):
        cx += p[0] * w
        cy += p[1] * w
        cz += p[2] * w
    c = (cx / total, cy / total, cz / total)
    xx = xy = xz = yy = yz = zz = 0.0
    for p, w in zip(points, weights):
        dx, dy, dz = p[0] - c[0], p[1] - c[1], p[2] - c[2]
        xx += w * dx * dx
        xy += w * dx * dy
        xz += w * dx * dz
        yy += w * dy * dy
        yz += w * dy * dz
        zz += w * dz * dz
    m = [[xx / total, xy / total, xz / total],
         [xy / total, yy / total, yz / total],
         [xz / total, yz / total, zz / total]]
    return c, m


def plane_fit(points, weights=None):
    """Least-squares plane through points.

    Returns (origin, normal, rms, flatness):
      origin    the weighted centre
      normal    unit normal (direction of least spread; its sign is arbitrary)
      rms       root-mean-square distance of the points from the plane
      flatness  smallest / middle spread, 0 for perfectly flat points. Near 1
                means the points do not define a plane (a ball or a line).
    """
    if len(points) < 3:
        raise ValueError("plane fit needs at least 3 points")
    origin, m = covariance(points, weights)
    values, vectors = eig_sym3(m)
    smallest = max(values[0], 0.0)
    middle = values[1]
    flatness = smallest / middle if middle > TINY else 1.0
    return origin, vectors[0], math.sqrt(smallest), flatness


def least_variance_direction(directions, weights=None):
    """Unit direction most square to a set of direction vectors.

    The smallest eigenvector of sum(w * d d^T). Typical use: the face normals
    of an elongated rod-like part are all roughly square to its long axis, so
    this returns that axis. Returns (direction, ratio) where ratio is the
    smallest / middle eigenvalue: near 0 the answer is sharp, near 1 there is
    no preferred direction (a ball-like part) and the caller should not trust it.
    Sign of the direction is arbitrary.
    """
    if not directions:
        raise ValueError("least variance of no directions")
    if weights is None:
        weights = [1.0] * len(directions)
    xx = xy = xz = yy = yz = zz = 0.0
    for d, w in zip(directions, weights):
        xx += w * d[0] * d[0]
        xy += w * d[0] * d[1]
        xz += w * d[0] * d[2]
        yy += w * d[1] * d[1]
        yz += w * d[1] * d[2]
        zz += w * d[2] * d[2]
    values, vectors = eig_sym3([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])
    ratio = values[0] / values[1] if values[1] > TINY else 1.0
    return vectors[0], max(ratio, 0.0)
