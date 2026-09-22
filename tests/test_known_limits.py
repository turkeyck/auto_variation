"""docs/SPEC.md test case: test_known_limits (H3 catalog integrity)."""
import pytest

from verification.run_all_checks import run_sympy_layer_module
from verification.known_limits import covered_limits


@pytest.mark.parametrize('entry', covered_limits(), ids=[e['name'] for e in covered_limits()])
def test_known_limit_still_present_and_passing(entry):
    results = run_sympy_layer_module(entry['module'])
    matches = [(name, ok) for name, ok, _ in results if entry['check_substring'] in name]
    assert matches, (
        f"known_limits.py catalogs '{entry['name']}' as covered by "
        f"{entry['module']}.run_checks(), but no check name there contains "
        f"'{entry['check_substring']}' anymore -- was it renamed or removed?")
    assert all(ok for _, ok in matches), (
        f"known limit '{entry['name']}' ({entry['module']}) is failing: {matches}")
