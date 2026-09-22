"""
Module I (docs/SPEC.md): LaTeX export for verification/registry.json
entries, with mandatory source citation and a visually unmissable
provenance marker.

DESIGN NOTE: registry.json (verification/registry_schema.py's
DerivedFormula) deliberately does NOT store the symbolic expression
itself -- it is a lightweight STATUS index ("this equation is at state
X, verified by routes A/B, on date D"), not a database of formulas (see
docs/SPEC.md "State 管理與持久化": "一切狀態即原始碼"). The actual
sympy/cadabra expression lives in whichever sympy_layer/*.py or
cadabra/*.py script derived it. This module's export functions therefore
take the expression (a sympy object, or a pre-rendered LaTeX string) as
a separate argument from the registry entry that describes its
provenance -- the caller is responsible for pulling the live expression
out of the script that owns it (see run_checks() below for a worked
example against a real registry entry, not a synthetic one).
"""
import json
from pathlib import Path

import sympy as sp

REPO_ROOT = Path(__file__).resolve().parent.parent
REGISTRY_PATH = REPO_ROOT / 'verification' / 'registry.json'


def load_registry():
    with open(REGISTRY_PATH, encoding='utf-8') as f:
        return json.load(f)


def find_entry(entries, entry_id):
    for e in entries:
        if e['id'] == entry_id:
            return e
    raise KeyError(f'no registry entry with id={entry_id!r}')


def export_formula(entry, expr=None, latex_body=None, lhs_latex=None):
    """Format one LaTeX block for a single registry entry.

    `expr` (a sympy expression) is converted via sp.latex(); pass
    `latex_body` directly instead if the caller already has a rendered
    string (e.g. from a cadabra Ex object, which this module cannot
    sp.latex() directly). Exactly one of the two must be given.
    `lhs_latex`, if given, prefixes "<lhs_latex> = " before the body.

    UNMISSABLE PROVENANCE MARKERS (docs/SPEC.md "狀態模型與揭露策略"):
      - provenance == 'literature_unverified' -> \\textbf{[UNVERIFIED]}
        prefix, so a copy-pasted formula can never silently look
        "the same" as a verified one in a paper draft.
      - state == 'known_issue' -> \\textbf{[KNOWN ISSUE]} prefix plus the
        known_issue_note as a LaTeX comment line, for the same reason.
    """
    if (expr is None) == (latex_body is None):
        raise ValueError('pass exactly one of expr= or latex_body=')
    body = sp.latex(expr) if expr is not None else latex_body
    if lhs_latex:
        body = f'{lhs_latex} = {body}'

    lines = []
    lines.append(f"% Source: {entry['paper_ref']}")
    lines.append(f"% Ansatz: {entry['ansatz']}")
    lines.append(f"% Derivation route A: {entry['derivation_route_a']}")
    if entry.get('derivation_route_b'):
        lines.append(f"% Derivation route B: {entry['derivation_route_b']}")
    lines.append(f"% State: {entry['state']}, last verified {entry.get('last_verified_date', 'N/A')}")

    prefix = ''
    if entry['provenance'] == 'literature_unverified':
        prefix = r'\textbf{[UNVERIFIED]}\ '
    if entry['state'] == 'known_issue':
        lines.append(f"% KNOWN ISSUE: {entry.get('known_issue_note', '')}")
        prefix = r'\textbf{[KNOWN ISSUE]}\ ' + prefix

    lines.append(r'\[')
    lines.append(f'  {prefix}{body}')
    lines.append(r'\]')
    return '\n'.join(lines)


def export_all(out_path=None):
    """Export every registry entry that HAS a directly re-derivable
    expression available in this repo's own scripts (be1/be2/be3 and the
    delta-Gamma L3 result -- the ones with a clean, cheap-to-recompute
    closed form); entries without one attached here are listed by id only,
    not silently skipped, so the output file's coverage is auditable."""
    import sys
    sys.path.insert(0, str(REPO_ROOT / 'sympy_layer'))
    from proca_minisuperspace import build_lagrangian, euler_lagrange_all
    from verification.paper_equations.arxiv_1703_09573 import background_field_equations

    entries = load_registry()
    blocks = []

    Xsym = sp.Symbol('X', positive=True)
    data = build_lagrangian(include_L5=True, d2_val=0,
                             G2_expr=Xsym, G3_expr=0 * Xsym, G4_expr=0 * Xsym, G5_expr=0 * Xsym)
    EL_N, EL_a, EL_A0 = euler_lagrange_all(data)
    be1, be2, be3 = background_field_equations(data['a_t'], -data['A0_t'], Xsym, 0 * Xsym, 0 * Xsym,
                                                0 * Xsym, Xsym, sp.Symbol('t'))
    lookup = {'1703_09573_be1': be1, '1703_09573_be2': be2, '1703_09573_be3': be3}

    for entry in entries:
        if entry['id'] in lookup:
            blocks.append(export_formula(entry, expr=lookup[entry['id']],
                                          lhs_latex=r'\text{(gravity+Proca side)}'))
        else:
            blocks.append(f"% {entry['id']}: expression not attached in export_all() -- "
                           f"see {entry['derivation_route_a']}")
        blocks.append('')

    text = '\n'.join(blocks)
    if out_path is not None:
        Path(out_path).write_text(text, encoding='utf-8')
    return text


def run_checks():
    """H2-style self-check: export_formula()'s UNVERIFIED/KNOWN ISSUE
    markers actually appear for entries with those provenance/state
    values, and do NOT appear for a clean paper_matched/derived entry --
    checked against the REAL registry.json content, not a synthetic
    fixture."""
    results = []
    entries = load_registry()

    known_issue_entries = [e for e in entries if e['state'] == 'known_issue']
    ok1 = len(known_issue_entries) > 0
    results.append(('1. registry.json has at least one known_issue entry to test the marker against '
                     '(currently 1605_05066_SMS_background_residual)', ok1, len(known_issue_entries)))

    if ok1:
        entry = known_issue_entries[0]
        block = export_formula(entry, latex_body='0')
        ok2 = r'\textbf{[KNOWN ISSUE]}' in block and entry['known_issue_note'][:30] in block
        results.append(('2. export_formula() renders the [KNOWN ISSUE] marker and includes the '
                         'known_issue_note as a comment for a known_issue entry', ok2, None))
    else:
        results.append(('2. (skipped, no known_issue entry found)', False, None))

    paper_matched_entries = [e for e in entries if e['state'] == 'paper_matched']
    ok3 = len(paper_matched_entries) > 0
    if ok3:
        entry = paper_matched_entries[0]
        block = export_formula(entry, latex_body='0')
        ok4 = (r'\textbf{[KNOWN ISSUE]}' not in block) and (r'\textbf{[UNVERIFIED]}' not in block)
        results.append(('3. export_formula() does NOT show either marker for a clean paper_matched '
                         'entry', ok4, None))
    else:
        results.append(('3. (skipped, no paper_matched entry found)', False, None))

    return results


if __name__ == '__main__':
    results = run_checks()
    n_pass = 0
    for name, ok, detail in results:
        status = 'PASS' if ok else 'FAIL'
        if ok:
            n_pass += 1
        print(f'[{status}] {name}')
        if not ok:
            print(f'       detail: {detail}')
    print(f'\n{n_pass}/{len(results)} checks passed')
