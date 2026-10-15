"""Repository rules, checked on every test run.

1. The library parses as code IronPython 2.7 can run: no syntax that only
   Python 3 has (f-strings, nonlocal, annotations, keyword-only arguments,
   math.inf, max/min default=, starred assignment, walrus, matrix multiply),
   and true division switched on in every module that divides.
2. core/ stays pure: no Rhino, .NET, Grasshopper, numpy or scipy imports.
3. Word list: words from a list kept outside the repository must not
   appear in it. The list comes from the SOLIDWRIGHT_WORDS environment
   variable (comma separated) or the file ~/.solidwright-words (one word per
   line); without either, that check is skipped.
4. Style: no em dash characters anywhere.
"""

import ast
import io
import os
import re

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PACKAGE = os.path.join(ROOT, "meshwright")
THIS_FILE = os.path.abspath(__file__)

SKIP_DIRS = {".git", ".venv", "venv", "__pycache__", ".pytest_cache", "build", "dist"}
TEXT_EXT = {".py", ".md", ".txt", ".toml", ".cfg", ".yml", ".yaml", ".ini", ".json", ""}


def _load_words():
    text = os.environ.get("SOLIDWRIGHT_WORDS", "").replace(",", "\n")
    path = os.path.join(os.path.expanduser("~"), ".solidwright-words")
    if not text.strip() and os.path.isfile(path):
        with io.open(path, encoding="utf-8") as f:
            text = f.read()
    return [w.strip() for w in text.splitlines() if w.strip()]


WORDS = _load_words()
WORDS_RE = re.compile(r"\b(" + "|".join(re.escape(w) for w in WORDS) + r")s?\b", re.IGNORECASE) if WORDS else None
EM_DASH = chr(0x2014)


def repo_files(exts=None):
    for base, dirs, files in os.walk(ROOT):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS and not d.endswith(".egg-info")]
        for name in files:
            path = os.path.join(base, name)
            if exts is None or os.path.splitext(name)[1] in exts:
                yield path


COMPONENTS = os.path.join(ROOT, "components")


def package_modules():
    """Library modules and Grasshopper component scripts: everything that runs inside Rhino."""
    return sorted(p for p in repo_files({".py"}) if p.startswith(PACKAGE) or p.startswith(COMPONENTS))


def read(path):
    with io.open(path, encoding="utf-8") as f:
        return f.read()


def rel(path):
    return os.path.relpath(path, ROOT)


# ---------------------------------------------------------------- 1. IronPython 2.7

class _Py3OnlyFinder(ast.NodeVisitor):
    def __init__(self):
        self.problems = []

    def _flag(self, node, what):
        self.problems.append("line %d: %s" % (getattr(node, "lineno", 0), what))

    def visit_JoinedStr(self, node):
        self._flag(node, "f-string")

    def visit_Nonlocal(self, node):
        self._flag(node, "nonlocal")

    def visit_AnnAssign(self, node):
        self._flag(node, "variable annotation")

    def visit_NamedExpr(self, node):
        self._flag(node, "walrus operator")

    def visit_YieldFrom(self, node):
        self._flag(node, "yield from")

    def visit_AsyncFunctionDef(self, node):
        self._flag(node, "async def")

    def visit_MatMult(self, node):
        self._flag(node, "@ operator")

    def visit_FunctionDef(self, node):
        args = node.args
        if args.kwonlyargs:
            self._flag(node, "keyword-only arguments")
        if node.returns is not None:
            self._flag(node, "return annotation")
        every = args.args + args.kwonlyargs + getattr(args, "posonlyargs", [])
        if args.vararg is not None:
            every.append(args.vararg)
        if args.kwarg is not None:
            every.append(args.kwarg)
        for a in every:
            if a.annotation is not None:
                self._flag(node, "argument annotation")
                break
        if getattr(args, "posonlyargs", []):
            self._flag(node, "positional-only arguments")
        self.generic_visit(node)

    def visit_Lambda(self, node):
        if node.args.kwonlyargs:
            self._flag(node, "keyword-only arguments")
        self.generic_visit(node)

    def visit_Attribute(self, node):
        if isinstance(node.value, ast.Name) and node.value.id == "math" and node.attr in ("inf", "nan", "tau", "isclose"):
            self._flag(node, "math.%s" % node.attr)
        self.generic_visit(node)

    def visit_Call(self, node):
        if isinstance(node.func, ast.Name) and node.func.id in ("max", "min"):
            if any(k.arg == "default" for k in node.keywords):
                self._flag(node, "%s(..., default=)" % node.func.id)
        self.generic_visit(node)

    def visit_Starred(self, node):
        if isinstance(node.ctx, ast.Store):
            self._flag(node, "starred assignment")
        self.generic_visit(node)


@pytest.mark.parametrize("path", package_modules(), ids=rel)
def test_module_is_ironpython_compatible(path):
    source = read(path)
    tree = ast.parse(source)
    finder = _Py3OnlyFinder()
    finder.visit(tree)
    assert not finder.problems, "%s uses Python 3 only syntax: %s" % (rel(path), finder.problems)
    divides = any(isinstance(n, (ast.BinOp, ast.AugAssign)) and isinstance(n.op, ast.Div)
                  for n in ast.walk(tree))
    if divides:
        futures = set()
        for node in tree.body:
            if isinstance(node, ast.ImportFrom) and node.module == "__future__":
                futures.update(a.name for a in node.names)
        assert "division" in futures, "%s uses / without from __future__ import division" % rel(path)


# ---------------------------------------------------------------- 2. pure core

FORBIDDEN_IN_CORE = {"Rhino", "System", "clr", "Grasshopper", "ghpythonlib", "rhinoscriptsyntax",
                     "scriptcontext", "numpy", "scipy"}


@pytest.mark.parametrize("path", [p for p in package_modules() if os.sep + "core" + os.sep in p], ids=rel)
def test_core_imports_nothing_external(path):
    tree = ast.parse(read(path))
    for node in ast.walk(tree):
        names = []
        if isinstance(node, ast.Import):
            names = [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            names = [node.module]
        for name in names:
            assert name.split(".")[0] not in FORBIDDEN_IN_CORE, "%s imports %s" % (rel(path), name)


# ---------------------------------------------------------------- 3 and 4. words and dashes

def _text_files():
    return [p for p in repo_files(TEXT_EXT) if os.path.abspath(p) != THIS_FILE]


def test_no_listed_words():
    if WORDS_RE is None:
        pytest.skip("no word list: set SOLIDWRIGHT_WORDS or write ~/.solidwright-words")
    hits = []
    for path in _text_files():
        for n, line in enumerate(read(path).splitlines(), 1):
            if WORDS_RE.search(line):
                hits.append("%s:%d" % (rel(path), n))  # never the word: CI logs are public
    assert not hits, "words from the list found on these lines:\n" + "\n".join(hits)


def test_no_em_dashes():
    hits = []
    for path in repo_files(TEXT_EXT):
        for n, line in enumerate(read(path).splitlines(), 1):
            if EM_DASH in line:
                hits.append("%s:%d" % (rel(path), n))
    assert not hits, "em dashes found:\n" + "\n".join(hits)
