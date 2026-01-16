import os
import sys

REPO_ROOT = "/Users/ootagakitakumi/RNAknotDetector"
PYTHON_DIR = os.path.join(REPO_ROOT, "python")

if PYTHON_DIR not in sys.path:
    sys.path.insert(0, PYTHON_DIR)

try:
    cmd.do(f"run {os.path.join(PYTHON_DIR, 'pymol_debug.py')}")
except Exception as exc:
    print(f"Failed to load pymol_debug.py: {exc}")
