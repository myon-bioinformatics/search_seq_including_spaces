"""Keep the documented source-only container and Python examples runnable."""
import ast
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class ServiceDistributionTests(unittest.TestCase):
    def test_docker_copy_layout_starts_service_without_site_packages(self):
        # Exercise exactly the Dockerfile's simple COPY source selections, without
        # requiring a daemon or pulling an image during the stdlib test suite.
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory)
            for line in (ROOT / 'examples/service/Dockerfile').read_text().splitlines():
                if not line.startswith('COPY '):
                    continue
                _, source, destination = line.split()
                matches = list(ROOT.glob(source))
                self.assertTrue(matches, source)
                for item in matches:
                    output = target / destination
                    if item.is_dir():
                        shutil.copytree(item, output, dirs_exist_ok=True)
                    else:
                        output.mkdir(parents=True, exist_ok=True)
                        shutil.copy2(item, output / item.name)
            run = subprocess.run([sys.executable, '-S', str(target / 'sequence_service.py'), '--help'],
                                 cwd=target, capture_output=True, text=True, timeout=20)
            self.assertEqual(run.returncode, 0, run.stderr)
            self.assertIn('http', run.stdout)

    def test_helper_package_imports_without_root_aliases(self):
        modules = sorted(path.stem for path in (ROOT / 'seqlib').glob('*.py')
                         if path.stem != '__init__')
        self.assertEqual(len(modules), 19)
        self.assertEqual({path.name for path in ROOT.glob('*.py')},
                         {'sequence_tool.py', 'sequence_service.py', 'definition_main.py'})
        code = 'import importlib\n' + '\n'.join(
            'importlib.import_module(' + repr('seqlib.' + name) + ')'
            for name in modules)
        run = subprocess.run([sys.executable, '-S', '-c', code], cwd=ROOT,
                             capture_output=True, text=True, timeout=20)
        self.assertEqual(run.returncode, 0, run.stderr)

    def test_documented_library_imports_resolve(self):
        for name in ('anchor-matcher', 'sequence-readers', 'sequence-matcher'):
            document = (ROOT / 'docs' / (name + '.md')).read_text()
            imports = []
            for block in re.findall(r'```python\n(.*?)```', document, re.S):
                tree = ast.parse(block)
                imports.extend(ast.unparse(node) for node in ast.walk(tree)
                               if isinstance(node, (ast.Import, ast.ImportFrom)))
            self.assertTrue(imports, name)
            run = subprocess.run([sys.executable, '-S', '-c', '\n'.join(imports)],
                                 cwd=ROOT, capture_output=True, text=True, timeout=20)
            self.assertEqual(run.returncode, 0, name + ': ' + run.stderr)


if __name__ == '__main__':
    unittest.main()
