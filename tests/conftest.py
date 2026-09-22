"""pytest path setup: put sympy_layer/, cadabra/, and the repo root on
sys.path so tests can import project modules directly, matching the same
sys.path handling verification/run_all_checks.py already uses."""
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SYMPY_LAYER = REPO_ROOT / 'sympy_layer'
CADABRA_DIR = REPO_ROOT / 'cadabra'

for p in (REPO_ROOT, SYMPY_LAYER, CADABRA_DIR):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
