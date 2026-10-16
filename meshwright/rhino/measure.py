"""Sizes of Rhino geometry for meshwright.core.cluster1d.

A closed mesh or a closed Brep is measured by its volume, anything else by
its area. Volume works for parts in any orientation (tilted parts on a
curve sort correctly); when some parts are open, every part is measured by
area instead so all sizes stay comparable.
"""

from __future__ import absolute_import, division

import Rhino.Geometry as rg

from ..core.report import Report


def is_closed(geometry):
    if isinstance(geometry, rg.Mesh):
        return bool(geometry.IsClosed)
    if isinstance(geometry, rg.Brep):
        return bool(geometry.IsSolid)
    return False


def volume(geometry):
    try:
        props = rg.VolumeMassProperties.Compute(geometry)
    except Exception:
        return None
    return None if props is None else abs(props.Volume)


def area(geometry):
    try:
        props = rg.AreaMassProperties.Compute(geometry)
    except Exception:
        return None
    return None if props is None else abs(props.Area)


def sizes(items):
    """(measures, power, report) for a list of meshes or Breps (None allowed).

    power turns a measure into a linear size: 1/3 for volumes, 1/2 for areas.
    """
    report = Report("measure")
    present = [g for g in items if g is not None]
    open_ones = [i for i, g in enumerate(items) if g is not None and not is_closed(g)]
    if present and not open_ones:
        return [None if g is None else volume(g) for g in items], 1.0 / 3.0, report
    if present and len(open_ones) < len(present):
        report.warn("input%s %s open: every part measured by area instead of volume" % (
            "" if len(open_ones) == 1 else "s", ", ".join(str(i) for i in open_ones) + (
                " is" if len(open_ones) == 1 else " are")))
    return [None if g is None else area(g) for g in items], 0.5, report


def extent(geometry, axis):
    """(low, high) of a piece of geometry along an axis (a plain 3-tuple), from
    its bounding box in a plane square to the axis; None when it has none."""
    if geometry is None:
        return None
    plane = rg.Plane(rg.Point3d(0.0, 0.0, 0.0), rg.Vector3d(axis[0], axis[1], axis[2]))
    try:
        box = geometry.GetBoundingBox(plane)
    except Exception:
        return None
    if box is None or (hasattr(box, "IsValid") and not box.IsValid):
        return None
    return box.Min.Z, box.Max.Z
