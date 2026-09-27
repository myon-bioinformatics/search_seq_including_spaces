"""Import-direction and no-dependency guard for the protein package.

Every module under protein/ (subpackages included) must be registered in
LAYERS by its dotted name relative to the package ("parse", "sub.mod"; a
subpackage's __init__.py is "sub", the package root is "__init__"). A module
may import only the standard library and package modules in a lower layer
(plus the explicit same-layer edges in SAME_LAYER_ALLOWED). Nothing may run
I/O or print at import time.
"""

import ast
import pathlib
import subprocess
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parent.parent
PACKAGE = ROOT / "protein"
ROOT_INIT = "__init__"

# Lower number = lower layer. Modules not written yet are listed so the
# intended structure is enforced as soon as they appear.
LAYERS = {
    "residues": 0,
    "provenance": 0,
    "parse": 1,
    "features": 2,
    "hints": 3,
    "roi": 3,
    "ss": 4,
}
SAME_LAYER_ALLOWED = {("roi", "hints")}
# Modules that may not import any package module (package init, output sink).
NO_PACKAGE_IMPORTS = {ROOT_INIT, "render"}

FORBIDDEN_TOP_LEVEL_CALLS = {"print", "open", "input", "exec", "eval", "breakpoint",
                             "__import__"}
FORBIDDEN_ANY_CALLS = {"__import__", "import_module"}


def _name_from_parts(parts):
    return ".".join(parts) or ROOT_INIT


def _parts_from_name(name):
    return [] if name == ROOT_INIT else name.split(".")


def _path_for(parts):
    """File of a module given as parts relative to protein/, or None."""
    base = PACKAGE.joinpath(*parts)
    if (base / "__init__.py").is_file():
        return base / "__init__.py"
    if parts and base.with_suffix(".py").is_file():
        return base.with_suffix(".py")
    return None


def package_modules():
    names = []
    for path in PACKAGE.rglob("*.py"):
        parts = list(path.relative_to(PACKAGE).with_suffix("").parts)
        if parts[-1] == "__init__":
            parts = parts[:-1]
        names.append(_name_from_parts(parts))
    return sorted(names)


def parse_module(name):
    path = _path_for(_parts_from_name(name))
    return ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def _resolve(base, alias):
    """`from <base> import <alias>`: a submodule if one exists, else base."""
    candidate = base + [alias]
    return _name_from_parts(candidate if _path_for(candidate) else base)


def collect_imports(tree, name):
    """Return (external_top_level_names, package_module_names) for module `name`."""
    parts = _parts_from_name(name)
    is_package = _path_for(parts).name == "__init__.py"
    package_parts = parts if is_package else parts[:-1]
    external, internal = set(), set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                dotted = alias.name.split(".")
                if dotted[0] == "protein":
                    internal.add(_name_from_parts(dotted[1:]))
                else:
                    external.add(dotted[0])
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                up = node.level - 1
                if up > len(package_parts):
                    internal.add("<outside-package>")
                    continue
                base = package_parts[:len(package_parts) - up]
                if node.module:
                    base = base + node.module.split(".")
            else:
                dotted = (node.module or "").split(".")
                if dotted[0] != "protein":
                    external.add(dotted[0])
                    continue
                base = dotted[1:]
            internal.update(_resolve(base, a.name) for a in node.names)
    return external, internal


def call_name(node):
    func = node.func
    if isinstance(func, ast.Name):
        return func.id
    if isinstance(func, ast.Attribute):
        return func.attr
    return None


def top_level_calls(tree):
    """Calls executed at import time (function bodies excluded)."""
    found = []

    def visit(node):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda)):
            return
        if isinstance(node, ast.Call):
            found.append(call_name(node))
        for child in ast.iter_child_nodes(node):
            visit(child)

    for stmt in tree.body:
        visit(stmt)
    return found


class TestLayers(unittest.TestCase):
    def test_every_module_is_registered(self):
        for name in package_modules():
            self.assertTrue(name in LAYERS or name in NO_PACKAGE_IMPORTS,
                            f"protein module {name!r} is not registered in LAYERS")

    def test_only_standard_library_is_imported(self):
        for name in package_modules():
            external, _ = collect_imports(parse_module(name), name)
            third_party = {m for m in external if m not in sys.stdlib_module_names}
            self.assertEqual(third_party, set(), f"protein module {name!r} imports {third_party}")

    def test_package_imports_point_to_lower_layers_only(self):
        for name in package_modules():
            _, internal = collect_imports(parse_module(name), name)
            internal.discard(name)
            if name in NO_PACKAGE_IMPORTS:
                self.assertEqual(internal, set(), f"protein module {name!r} must not import {internal}")
                continue
            self.assertIn(name, LAYERS, f"protein module {name!r} is not registered in LAYERS")
            for target in internal:
                self.assertIn(target, LAYERS, f"protein module {name!r} imports unknown {target!r}")
                ok = (LAYERS[target] < LAYERS[name]
                      or (name, target) in SAME_LAYER_ALLOWED)
                self.assertTrue(ok, f"protein module {name!r} (L{LAYERS[name]}) imports "
                                    f"{target!r} (L{LAYERS[target]})")

    def test_no_dynamic_imports(self):
        for name in package_modules():
            calls = {call_name(n) for n in ast.walk(parse_module(name))
                     if isinstance(n, ast.Call)}
            self.assertFalse(calls & FORBIDDEN_ANY_CALLS, f"protein module {name!r}")

    def test_no_io_at_import_time(self):
        for name in package_modules():
            calls = set(top_level_calls(parse_module(name)))
            self.assertFalse(calls & FORBIDDEN_TOP_LEVEL_CALLS,
                             f"protein module {name!r} calls "
                             f"{calls & FORBIDDEN_TOP_LEVEL_CALLS} at import")


class TestCommandLineEntryPoint(unittest.TestCase):
    """sequence_tool.py is the only place that combines package modules with I/O."""

    def test_imports_only_stdlib_and_protein(self):
        path = ROOT / "sequence_tool.py"
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        imported = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                self.assertEqual(node.level, 0, "sequence_tool.py uses a relative import")
                imported.add(node.module.split(".")[0])
        third_party = {m for m in imported
                       if m != "protein" and m not in sys.stdlib_module_names}
        self.assertEqual(third_party, set())
        self.assertIn("protein", imported)

    def test_package_does_not_import_the_entry_point(self):
        for name in package_modules():
            external, _ = collect_imports(parse_module(name), name)
            self.assertNotIn("sequence_tool", external, f"protein module {name!r}")


class TestRuntimeWithoutSitePackages(unittest.TestCase):
    def test_imports_cleanly_with_site_disabled(self):
        modules = ", ".join(f"protein.{m}" for m in package_modules() if m != ROOT_INIT)
        code = (
            "import sys\n"
            "assert not any(p.endswith(('site-packages', 'dist-packages')) for p in sys.path), sys.path\n"
            f"import {modules}\n"
        )
        result = subprocess.run([sys.executable, "-S", "-E", "-c", code], cwd=ROOT,
                                capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, "")
        self.assertEqual(result.stderr, "")
