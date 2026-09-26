"""Import-direction and no-dependency guard for the protein package.

Every module in protein/ must be registered in LAYERS. A module may import
only the standard library and package modules in a lower layer (plus the
explicit same-layer edges in SAME_LAYER_ALLOWED). Nothing may run I/O or
print at import time.
"""

import ast
import pathlib
import subprocess
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parent.parent
PACKAGE = ROOT / "protein"

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
NO_PACKAGE_IMPORTS = {"__init__", "render"}

FORBIDDEN_TOP_LEVEL_CALLS = {"print", "open", "input", "exec", "eval", "breakpoint",
                             "__import__"}
FORBIDDEN_ANY_CALLS = {"__import__", "import_module"}


def package_modules():
    return sorted(p.stem for p in PACKAGE.glob("*.py"))


def parse_module(name):
    path = PACKAGE / f"{name}.py"
    return ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def collect_imports(tree):
    """Return (external_top_level_names, package_module_names)."""
    external, internal = set(), set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                parts = alias.name.split(".")
                if parts[0] == "protein":
                    internal.add(parts[1] if len(parts) > 1 else "__init__")
                else:
                    external.add(parts[0])
        elif isinstance(node, ast.ImportFrom):
            if node.level > 1:
                internal.add("<outside-package>")
            elif node.level == 1:
                if node.module:
                    internal.add(node.module.split(".")[0])
                else:
                    internal.update(alias.name for alias in node.names)
            else:
                parts = (node.module or "").split(".")
                if parts[0] == "protein":
                    if len(parts) > 1:
                        internal.add(parts[1])
                    else:
                        internal.update(alias.name for alias in node.names)
                else:
                    external.add(parts[0])
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
                            f"protein/{name}.py is not registered in LAYERS")

    def test_only_standard_library_is_imported(self):
        for name in package_modules():
            external, _ = collect_imports(parse_module(name))
            third_party = {m for m in external if m not in sys.stdlib_module_names}
            self.assertEqual(third_party, set(), f"protein/{name}.py imports {third_party}")

    def test_package_imports_point_to_lower_layers_only(self):
        for name in package_modules():
            _, internal = collect_imports(parse_module(name))
            internal.discard(name)
            if name in NO_PACKAGE_IMPORTS:
                self.assertEqual(internal, set(), f"protein/{name}.py must not import {internal}")
                continue
            for target in internal:
                self.assertIn(target, LAYERS, f"protein/{name}.py imports unknown {target!r}")
                ok = (LAYERS[target] < LAYERS[name]
                      or (name, target) in SAME_LAYER_ALLOWED)
                self.assertTrue(ok, f"protein/{name}.py (L{LAYERS[name]}) imports "
                                    f"{target} (L{LAYERS[target]})")

    def test_no_dynamic_imports(self):
        for name in package_modules():
            calls = {call_name(n) for n in ast.walk(parse_module(name))
                     if isinstance(n, ast.Call)}
            self.assertFalse(calls & FORBIDDEN_ANY_CALLS, f"protein/{name}.py")

    def test_no_io_at_import_time(self):
        for name in package_modules():
            calls = set(top_level_calls(parse_module(name)))
            self.assertFalse(calls & FORBIDDEN_TOP_LEVEL_CALLS,
                             f"protein/{name}.py calls {calls & FORBIDDEN_TOP_LEVEL_CALLS} at import")


class TestRuntimeWithoutSitePackages(unittest.TestCase):
    def test_imports_cleanly_with_site_disabled(self):
        modules = ", ".join(f"protein.{m}" for m in package_modules() if m != "__init__")
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


if __name__ == "__main__":
    unittest.main()
