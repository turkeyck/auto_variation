"""docs/SPEC.md test case: test_no_rename_dummies_violation."""
from check_no_rename_dummies import find_violations


def test_no_rename_dummies_violation():
    violations = find_violations()
    assert violations == [], (
        'rename_dummies() is banned outside cadabra_utils.py (D2) -- '
        f'found: {violations}')
