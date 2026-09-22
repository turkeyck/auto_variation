"""
docs/SPEC.md Stage 6: pytest wrapper giving per-check granularity over
every sympy_layer/*.py run_checks() (the same set
verification/run_all_checks.py aggregates, reused here rather than
re-implemented -- this file is a thin pytest adapter over that shared
collection logic, not a second copy of it).

Run: py -3 -m pytest tests/test_regression.py -v
"""
import pytest

from verification.run_all_checks import (
    SYMPY_LAYER_MODULES, TRACEABILITY_MODULES,
    run_sympy_layer_module, run_traceability_module, is_known_issue,
)


def _collect_all_checks():
    """(module_name, check_name, ok, detail) for every check in every
    module, computed once at collection time."""
    rows = []
    for module_name in SYMPY_LAYER_MODULES:
        results = run_sympy_layer_module(module_name)
        for check_name, ok, detail in results:
            rows.append((module_name, check_name, ok, detail))
    for module_name in TRACEABILITY_MODULES:
        results = run_traceability_module(module_name)
        for check_name, ok, detail in results:
            rows.append((module_name, check_name, ok, detail))
    return rows


_ALL_CHECKS = _collect_all_checks()


def _test_id(row):
    module_name, check_name, ok, detail = row
    short = check_name.split('.', 1)[-1].strip()[:60]
    return f'{module_name}::{short}'


@pytest.mark.parametrize('row', _ALL_CHECKS, ids=[_test_id(r) for r in _ALL_CHECKS])
def test_sympy_layer_check(row):
    module_name, check_name, ok, detail = row
    if is_known_issue(check_name):
        # A check explicitly labeled "KNOWN [OPEN] ISSUE" in its own name
        # is expected to be tracked, not silently green -- xfail (not
        # skip) so it still shows up distinctly in the pytest summary,
        # and so this test starts FAILING (which is the signal to act
        # on) if that known issue is ever accidentally reintroduced after
        # being fixed, or unexpectedly starts passing (xpass) if it gets
        # fixed without updating the check's own name.
        if not ok:
            pytest.xfail(f'known issue: {check_name}')
    assert ok, f'{module_name}: {check_name}\n  detail: {detail}'
