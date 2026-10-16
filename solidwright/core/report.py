"""Graded results.

Every solidwright tool that can go wrong returns a Report next to its result,
instead of raising or quietly handing back bad geometry. A Report holds a
status and the reasons behind it:

    0  GOOD     nothing to say
    1  WARNING  usable, but something deserves a look
    2  MESSY    a result was produced from questionable input or by a fallback
    3  FAILED   no trustworthy result

The status is the worst severity among the reasons, so adding reasons can only
make a report worse, never better.
"""

from __future__ import absolute_import, division

GOOD = 0
WARNING = 1
MESSY = 2
FAILED = 3

STATUS_NAMES = {
    GOOD: "good",
    WARNING: "warning",
    MESSY: "messy",
    FAILED: "failed",
}


class Report(object):
    """Status 0..3 plus the readable reasons that produced it."""

    def __init__(self, title=""):
        self.title = title
        self.reasons = []  # list of (severity, text), in the order they were added

    @property
    def status(self):
        worst = GOOD
        for severity, _ in self.reasons:
            if severity > worst:
                worst = severity
        return worst

    @property
    def ok(self):
        """True when the result can be used as is (good or warning)."""
        return self.status <= WARNING

    @property
    def failed(self):
        return self.status >= FAILED

    def add(self, severity, text):
        """Record one reason. Severity GOOD is allowed for plain notes."""
        severity = int(severity)
        if severity < GOOD or severity > FAILED:
            raise ValueError("severity must be 0..3, got %r" % (severity,))
        self.reasons.append((severity, text))
        return self

    def note(self, text):
        return self.add(GOOD, text)

    def warn(self, text):
        return self.add(WARNING, text)

    def messy(self, text):
        return self.add(MESSY, text)

    def fail(self, text):
        return self.add(FAILED, text)

    def merge(self, other, prefix=""):
        """Copy another report's reasons into this one, optionally prefixed."""
        for severity, text in other.reasons:
            self.reasons.append((severity, prefix + text if prefix else text))
        return self

    def worst_reasons(self):
        """The reasons that set the current status."""
        s = self.status
        return [text for severity, text in self.reasons if severity == s and s > GOOD]

    def summary(self):
        """One line: status name plus the reasons that set it."""
        s = self.status
        head = "%d %s" % (s, STATUS_NAMES[s])
        if self.title:
            head = "%s: %s" % (self.title, head)
        worst = self.worst_reasons()
        if worst:
            head += " (" + "; ".join(worst) + ")"
        return head

    def __str__(self):
        lines = [self.summary()]
        for severity, text in self.reasons:
            lines.append("  [%d] %s" % (severity, text))
        return "\n".join(lines)

    def __repr__(self):
        return "Report(status=%d, reasons=%d)" % (self.status, len(self.reasons))
