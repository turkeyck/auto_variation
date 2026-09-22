"""
Module H: one-shot aggregate runner over every sympy_layer/*.py
run_checks() (13 modules) plus the cadabra/ lint check, producing a
human-readable report and a machine-readable JSON summary.

docs/SPEC.md "狀態模型與揭露策略": known_issue-tagged failures must always
be shown first, never buried under passing checks -- this is the direct
fix for the historical failure mode where "6/7 checks pass" reporting let
the 1 open issue go unnoticed.

Scope (Stage 0, docs/SPEC.md): this script runs the sympy_layer/ half
unconditionally (importable from the Windows-side `py -3` interpreter).
The cadabra/ half's run_checks()-equivalents require cadabra2, which only
exists in WSL and cannot be imported from here -- those are reported as
explicitly skipped with the reason, not silently omitted. The one part of
cadabra/ that IS pure Python (check_no_rename_dummies.py's lint) is run
directly, since it needs no cadabra2 kernel.

Usage:
    py -3 verification/run_all_checks.py [--verbose]
"""
import argparse
import importlib
import json
import sys
from pathlib import Path
from datetime import date

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
SYMPY_LAYER = REPO_ROOT / 'sympy_layer'
CADABRA_DIR = REPO_ROOT / 'cadabra'
TRACEABILITY_DIR = REPO_ROOT / 'traceability'

# Every sympy_layer/*.py file that defines run_checks(). tensor_utils.py is
# a pure utility module (Christoffel/Riemann/Ricci helpers) with no
# run_checks() of its own -- it is exercised indirectly via fR_gravity.py
# and proca_minisuperspace.py, so it is intentionally excluded here.
SYMPY_LAYER_MODULES = [
    'proca_minisuperspace',
    'proca_background_1703_09573',
    'proca_background_direct_route',
    'delta_gamma_check',
    'proca_L6_background',
    'proca_tensor',
    'fR_gravity',
    'vector_sector',
    'vector_field_equation',
    'adm_scalar',
    'scalar_sector',
    'scalar_plus_fluid_sector',
    'schutz_sorkin_vector',
    'schutz_sorkin_scalar',
    'fluid_sector',
    'kinetic_matrix',
    'svt_builder',
    'pert_engine',
    'linear_order_check',
]

# traceability/latex_export.py lives outside sympy_layer/ but follows the
# same run_checks() contract; imported separately since it needs the repo
# root (for verification/registry.json), not sympy_layer/, on sys.path.
TRACEABILITY_MODULES = ['latex_export']

KNOWN_ISSUE_MARKERS = ('KNOWN OPEN ISSUE', 'KNOWN ISSUE')


def run_sympy_layer_module(module_name):
    """Import one sympy_layer module and call its run_checks(). Returns
    exactly what the module's own run_checks() returns -- no
    reinterpretation, so a failure here means the underlying module's own
    check failed, not an aggregator bug."""
    if str(SYMPY_LAYER) not in sys.path:
        sys.path.insert(0, str(SYMPY_LAYER))
    mod = importlib.import_module(module_name)
    return mod.run_checks()


def run_traceability_module(module_name):
    """Same contract as run_sympy_layer_module(), for traceability/*.py
    (module I) -- kept as a separate function since it needs
    traceability/ (not sympy_layer/) on sys.path, plus the repo root
    (for its own `from verification.paper_equations...` import)."""
    if str(TRACEABILITY_DIR) not in sys.path:
        sys.path.insert(0, str(TRACEABILITY_DIR))
    if str(REPO_ROOT) not in sys.path:
        sys.path.insert(0, str(REPO_ROOT))
    if str(SYMPY_LAYER) not in sys.path:
        sys.path.insert(0, str(SYMPY_LAYER))
    mod = importlib.import_module(module_name)
    return mod.run_checks()


def run_cadabra_lint():
    """Run cadabra/check_no_rename_dummies.py's find_violations() directly
    -- it is pure-Python regex linting with no cadabra2 kernel dependency,
    so unlike the rest of cadabra/ it can run under this interpreter."""
    if str(CADABRA_DIR) not in sys.path:
        sys.path.insert(0, str(CADABRA_DIR))
    lint_mod = importlib.import_module('check_no_rename_dummies')
    return lint_mod.find_violations()


def is_known_issue(check_name):
    return any(marker in check_name for marker in KNOWN_ISSUE_MARKERS)


def collect_report():
    module_results = {}
    module_errors = {}
    for name in SYMPY_LAYER_MODULES:
        try:
            module_results[name] = run_sympy_layer_module(name)
        except Exception as exc:  # noqa: BLE001 -- keep going, report which module broke
            module_errors[name] = f'{type(exc).__name__}: {exc}'

    for name in TRACEABILITY_MODULES:
        try:
            module_results[name] = run_traceability_module(name)
        except Exception as exc:  # noqa: BLE001
            module_errors[name] = f'{type(exc).__name__}: {exc}'

    lint_violations = None
    lint_error = None
    try:
        lint_violations = run_cadabra_lint()
    except Exception as exc:  # noqa: BLE001
        lint_error = f'{type(exc).__name__}: {exc}'

    return module_results, module_errors, lint_violations, lint_error


def build_json_summary(module_results, module_errors, lint_violations, lint_error):
    known_issues, failures, passes = [], [], []
    for module_name, results in module_results.items():
        for check_name, ok, detail in results:
            entry = {'module': module_name, 'check': check_name}
            if not ok:
                (known_issues if is_known_issue(check_name) else failures).append(entry)
            else:
                passes.append(entry)

    return {
        'generated': date.today().isoformat(),
        'known_issues': known_issues,
        'unexpected_failures': failures,
        'pass_count': len(passes),
        'module_import_errors': module_errors,
        'cadabra_lint': {
            'error': lint_error,
            'violations': lint_violations or [],
        },
        'cadabra_run_checks': {
            'skipped': True,
            'reason': ('requires WSL + cadabra2 2.5.14; not importable from this interpreter '
                       '(see docs/SPEC.md 執行環境). Run individually via: '
                       'wsl.exe -- python3 cadabra/<file>.py'),
        },
    }


def print_report(module_results, module_errors, lint_violations, lint_error, verbose):
    def rule(title):
        print('=' * 78)
        print(title)
        print('=' * 78)

    rule('KNOWN ISSUES (always shown first, never buried under passing checks)')
    any_known = False
    for module_name, results in module_results.items():
        for check_name, ok, detail in results:
            if not ok and is_known_issue(check_name):
                any_known = True
                print(f'[KNOWN ISSUE] {module_name}: {check_name}')
                if verbose:
                    print(f'    detail: {detail}')
    if not any_known:
        print('(none)')
    print()

    rule('UNEXPECTED FAILURES')
    any_fail = False
    for module_name, results in module_results.items():
        for check_name, ok, detail in results:
            if not ok and not is_known_issue(check_name):
                any_fail = True
                print(f'[FAIL] {module_name}: {check_name}')
                print(f'    detail: {detail}')
    if not any_fail:
        print('(none)')
    print()

    rule('MODULE IMPORT ERRORS')
    if module_errors:
        for name, err in module_errors.items():
            print(f'[IMPORT ERROR] {name}: {err}')
    else:
        print('(none)')
    print()

    rule('PER-MODULE SUMMARY')
    total, passed = 0, 0
    for module_name, results in module_results.items():
        n_pass = sum(1 for _, ok, _ in results if ok)
        total += len(results)
        passed += n_pass
        print(f'  {module_name}: {n_pass}/{len(results)} passed')
    print(f'\nTOTAL (sympy_layer): {passed}/{total} passed', end='')
    if module_errors:
        print(f', {len(module_errors)} module(s) failed to import')
    else:
        print()
    print()

    rule('CADABRA LAYER')
    if lint_error:
        print(f'[SKIPPED] check_no_rename_dummies.py lint: {lint_error}')
    elif lint_violations:
        print(f'[LINT FAIL] {len(lint_violations)} rename_dummies() violation(s) found:')
        for name, lineno, line in lint_violations:
            print(f'  {name}:{lineno}: {line}')
    else:
        print('[LINT OK] no standalone rename_dummies() calls found.')
    print('[SKIPPED] cadabra run_checks() -- requires WSL + cadabra2 2.5.14, not run from this '
          'interpreter. Run individually via: wsl.exe -- python3 cadabra/<file>.py')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verbose', action='store_true', help='show detail for known issues too')
    args = parser.parse_args()

    module_results, module_errors, lint_violations, lint_error = collect_report()
    print_report(module_results, module_errors, lint_violations, lint_error, args.verbose)

    summary = build_json_summary(module_results, module_errors, lint_violations, lint_error)
    out_path = HERE / 'last_run_summary.json'
    out_path.write_text(json.dumps(summary, indent=2, default=str), encoding='utf-8')
    print(f'\nJSON summary written to {out_path}')

    if module_errors:
        sys.exit(1)


if __name__ == '__main__':
    main()
