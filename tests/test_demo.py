"""Exercise the documented demo command and its expected result tables."""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


class DemoTests(unittest.TestCase):
    @unittest.skipUnless(shutil.which("bash"), "Demo runner requires Bash")
    def test_demo_command_and_expected_outputs(self):
        samtools = os.environ.get("RIVER_TEST_SAMTOOLS") or shutil.which("samtools")
        if not samtools:
            self.skipTest("Set RIVER_TEST_SAMTOOLS to run the demo integration test")
        repo = Path(__file__).resolve().parents[1]
        env = dict(os.environ, RIVER_PYTHON=sys.executable, RIVER_SAMTOOLS=samtools)
        env["PYTHONPATH"] = str(repo / "src") + os.pathsep + env.get("PYTHONPATH", "")
        with tempfile.TemporaryDirectory() as temporary:
            # A directory containing spaces also checks shell/path handling.
            outdir = Path(temporary) / "synthetic demo"
            command = ["bash", str(repo / "examples/run_demo.sh"), str(outdir)]
            result = subprocess.run(command, env=env, cwd=temporary,
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn("Demo matches expected tables", result.stdout)
            saved = (outdir / "results/summary.tsv").read_bytes()
            repeated = subprocess.run(command, env=env, cwd=temporary,
                                      capture_output=True, text=True)
            self.assertNotEqual(repeated.returncode, 0)
            self.assertIn("already exists", repeated.stderr)
            self.assertEqual((outdir / "results/summary.tsv").read_bytes(), saved)
