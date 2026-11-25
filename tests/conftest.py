import os
import sys

# Make the license-free tools importable in tests and via their CLIs.
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))

EXAMPLE_DIR = os.path.join(ROOT, "examples", "egfr")
PDB_1M17 = os.path.join(EXAMPLE_DIR, "1M17.pdb")
