import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "2DPPIViewer.py"
EXAMPLE = REPO_ROOT / "tests" / "data" / "example.drw"


class Test2DPPIViewerSmoke(unittest.TestCase):
    def test_script_exists(self):
        self.assertTrue(SCRIPT.exists(), f"Missing script: {SCRIPT}")

    def test_example_exists(self):
        self.assertTrue(EXAMPLE.exists(), f"Missing example: {EXAMPLE}")

    def test_generates_html(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out_html = Path(tmpdir) / "out.html"

            cmd = [
                sys.executable,
                str(SCRIPT),
                str(EXAMPLE),
                "--out",
                str(out_html),
            ]
            res = subprocess.run(cmd, capture_output=True, text=True)

            self.assertEqual(
                res.returncode, 0,
                f"Command failed.\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
            )
            self.assertTrue(out_html.exists(), "Output HTML was not created.")
            self.assertGreater(out_html.stat().st_size, 1000, "Output HTML seems too small.")

            txt = out_html.read_text(encoding="utf-8", errors="ignore").lower()
            self.assertIn("<html", txt)
            self.assertIn("<svg", txt)


if __name__ == "__main__":
    unittest.main()
