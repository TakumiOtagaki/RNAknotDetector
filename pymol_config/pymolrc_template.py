import os
import sys

REPO_ROOT = "/Users/ootagakitakumi/RNAknotDetector"
PYTHON_DIR = os.path.join(REPO_ROOT, "python")
RNA_TOOLS_DIR = os.path.join(REPO_ROOT, "rna-tools")

if PYTHON_DIR not in sys.path:
    sys.path.insert(0, PYTHON_DIR)
if RNA_TOOLS_DIR not in sys.path:
    sys.path.insert(0, RNA_TOOLS_DIR)

from pymol import cmd

print("Loading RNAknotDetector PyMOL helpers...")
# cmd.do(f"run {os.path.join(PYTHON_DIR, 'pymol_debug.py')}")
cmd.do(f"run {os.path.join(PYTHON_DIR, 'pymol_rnaknot_hit.py')}")
try:
    from rna_tools import rna_tools_config as rtc
    wrapper = os.path.join(REPO_ROOT, "apps", "x3dna-dssr-wrapper")
    dssr = os.path.join(REPO_ROOT, "apps", "x3dna-dssr")
    rtc.X3DNA = wrapper if os.path.isfile(wrapper) else dssr
except Exception as exc:
    print(f"[warning] Failed to set DSSR path: {exc}")
