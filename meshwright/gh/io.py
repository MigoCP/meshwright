"""Input reading and status messages for Grasshopper script components.

In a GH script component every input is a global variable. If an input was
removed or renamed on the component, reading it directly raises NameError,
and an unconnected optional input arrives as None. read_input() handles both
by falling back to a default.
"""

from __future__ import absolute_import, division


def read_input(namespace, name, default=None):
    """namespace[name], or default when it is missing or None.

    Call it with the component's globals(): read_input(globals(), "Tol", 0.0).
    """
    value = namespace.get(name, None)
    if value is None:
        return default
    return value


def show_report(ghenv, report):
    """Show a Report as the component's balloon (warning or error colour).

    Good results show nothing, warnings and messy results show an orange
    balloon, failures a red one. Does nothing outside Grasshopper.
    """
    if ghenv is None or report.status == 0:
        return
    try:
        from Grasshopper.Kernel import GH_RuntimeMessageLevel as Level
    except Exception:
        return
    level = Level.Error if report.failed else Level.Warning
    ghenv.Component.AddRuntimeMessage(level, report.summary())
