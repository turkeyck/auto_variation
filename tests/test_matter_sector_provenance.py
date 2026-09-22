"""
docs/SPEC.md test case: test_schutz_sorkin_scalar_no_kessence_shortcut.

Static check (not a symbolic one -- this guards against a REGRESSION in
which sympy_layer/schutz_sorkin_scalar.py's genuine J^mu-based derivation
is quietly replaced by a k-essence P(X)=X duality shortcut again, exactly
the class of substitution docs/SPEC.md module G explicitly forbids and
that sympy_layer/scalar_plus_fluid_sector.py's own docstring documents
having used). A symbolic check cannot catch this (a k-essence stand-in
can numerically/symbolically resemble parts of the real derivation in
some limits -- that is precisely why the professor's spec insists on
provenance, not just numerical agreement).
"""
from pathlib import Path

SYMPY_LAYER = Path(__file__).resolve().parent.parent / 'sympy_layer'

REQUIRED_MARKERS = ['J0', 'Jz', 'JJg']  # genuine Schutz-Sorkin J^mu construction
FORBIDDEN_MARKERS = ['k-essence', 'P(X)=X', 'P(X) = X']  # duality-shortcut fingerprints


def _code_body_without_module_docstring(src):
    """Strip the leading module docstring (the FIRST '\"\"\"...\"\"\"' block)
    before scanning for forbidden markers -- mentioning "k-essence" in
    PROSE to explain what a file deliberately does NOT do (as this file's
    own docstring does, contrasting itself with
    scalar_plus_fluid_sector.py) is exactly the kind of honest
    provenance note docs/SPEC.md module I wants, not a violation; only a
    literal k-essence-style CODE pattern in the executable body should
    fail this check."""
    parts = src.split('"""')
    if len(parts) >= 3:
        return '"""'.join(parts[2:])  # everything after the closing docstring quote
    return src


def test_schutz_sorkin_scalar_no_kessence_shortcut():
    path = SYMPY_LAYER / 'schutz_sorkin_scalar.py'
    src = path.read_text(encoding='utf-8')
    code_body = _code_body_without_module_docstring(src)

    missing = [m for m in REQUIRED_MARKERS if m not in src]
    assert not missing, (
        f'{path.name} no longer builds J^mu explicitly (missing markers: {missing}) -- '
        'has it regressed to a duality shortcut instead of a genuine Schutz-Sorkin derivation?')

    present_forbidden = [m for m in FORBIDDEN_MARKERS if m in code_body]
    assert not present_forbidden, (
        f'{path.name}\'s CODE BODY (not its docstring) contains k-essence-duality language '
        f'({present_forbidden}) -- docs/SPEC.md module G forbids using a k-essence stand-in for '
        'the genuine Schutz-Sorkin scalar matter sector.')


def test_scalar_plus_fluid_sector_marked_superseded():
    """The k-essence-duality file itself must stay clearly labeled as
    superseded (not silently re-promoted to "the" scalar matter sector
    reference) -- see docs/SPEC.md Stage 4."""
    path = SYMPY_LAYER / 'scalar_plus_fluid_sector.py'
    src = path.read_text(encoding='utf-8')
    assert 'SUPERSEDED' in src, (
        f'{path.name} is missing its SUPERSEDED marker -- it must stay flagged as not '
        'satisfying the Schutz-Sorkin requirement (see schutz_sorkin_scalar.py instead).')
