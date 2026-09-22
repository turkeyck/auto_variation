"""
Shared helpers for working around the cadabra2 2.4.5.4 "Free indices in
different terms in a sum do not match" bug (see fR_variation.py for the
root-cause writeup: the checker cannot see contractions that cross a
\\delta{...} Accent boundary, so it misjudges otherwise-valid sums).

Workaround: never let two additive terms with a delta{...} factor sit
in the same Ex()-parsed string. Build every additive piece as its own
single-term Ex, transform it with a single-term substitute() rule, and
combine already-built trees with Python's `+` (never `Ex(str(tree))` --
that re-parses and re-triggers the bug; use `.copy()` instead).

SECOND, SEPARATE bug found empirically in this same cadabra2 2.4.5.4
install: `rename_dummies()` raises "No index set for index ... known"
on perfectly valid contractions of a plain (undeclared-property)
one-index tensor against an InverseMetric/Symmetric two-index tensor
(e.g. `g^{mu nu} A_mu A_nu`), regardless of what property (if any) is
declared on the one-index tensor. `canonicalise()` alone handles the
exact same expression without complaint (it evidently does its own
internal dummy-renaming). Workaround: skip `rename_dummies()` and call
`canonicalise()` directly wherever this occurs.

THIRD, SEPARATE bug found empirically in this same install: calling
`substitute()` with a replacement RHS that introduces a dummy pair
whose name (e.g. \mu) is already in use as an unrelated dummy pair
elsewhere in the target expression raises "Failed to find dummy
property for $\mu$ while renaming dummies" -- even though the two
pairs are logically unrelated. Workaround: never hand-pick index
letters for a new dummy pair (that requires tracking, by eye, every
name already "spoken for" elsewhere in the expression, which does not
scale to deeply nested derivative terms). Instead draw every new dummy
pair's names from `fresh_indices()` below, which hands out names
guaranteed to be distinct from every name this kernel session has
already produced.

Re-verified against upstream cadabra2 2.5.14 (2025-07-31, ~10 releases
past 2.4.5.4) on 2026-09-12: all three bugs are still present. This is
not a stale-package artifact -- treat the three workarounds below as
permanent operating rules for this project, not stopgaps to remove on
the next upgrade.

RULES (apply in every new cadabra script that touches metric variation
or nested covariant derivatives):
  1. Never let two terms that each carry a `\delta{...}` (Accent) factor
     appear together in one `Ex()`-parsed string or one `substitute()`
     RHS string. Build each additive piece as its own single-term
     Ex/substitute, then combine the built trees with Python `+`. Use
     `sub_copy`/`replace_term` below rather than re-deriving this by hand.
  2. Never call `rename_dummies()` on its own. Always use
     `canonicalise()` (it performs equivalent dummy-renaming internally
     without hitting the bug).
  3. Never hand-pick index letters for a dummy pair introduced by a
     `substitute()` RHS. Always draw the names from `fresh_indices()`.
"""
from cadabra2 import Ex, substitute

_fresh_index_counter = [0]


def fresh_indices(n, letter=r'\mu', position='fixed'):
    """Return a list of `n` brand-new index-name strings (e.g.
    ['\\mu_{0}', '\\mu_{1}']), guaranteed distinct from every name
    fresh_indices() has already handed out in this kernel session, and
    pre-declared via `{...}::Indices(position=...)` so they can be used
    immediately inside Ex()/substitute() strings.

    Always take dummy-pair names for a new substitute() RHS from this
    pool instead of hand-picking greek letters -- reusing a letter that
    is already a dummy elsewhere in the target expression is exactly
    what triggers cadabra2's dummy-collision bug (see module docstring,
    bug #3)."""
    start = _fresh_index_counter[0]
    names = [f'{letter}_{{{start + i}}}' for i in range(n)]
    _fresh_index_counter[0] = start + n
    Ex('{' + ','.join(names) + '}::Indices(position=' + position + ').')
    return names


def sub_copy(expr, rule_str):
    """Copy expr, apply a single substitute() rule, return the copy."""
    e = expr.copy()
    substitute(e, Ex(rule_str))
    return e


def replace_term(v, term_str, *replacement_strs):
    """Replace the single additive term `term_str` inside `v` with the
    sum of `replacement_strs`, each substituted in from its own isolated
    copy of `term_str` (so mismatched index patterns across the
    replacement pieces never appear in one parsed string)."""
    pieces = [sub_copy(Ex(term_str), f'{term_str} -> {r}') for r in replacement_strs]
    combined = pieces[0]
    for p in pieces[1:]:
        combined = combined + p
    v_without = sub_copy(v, f'{term_str} -> 0')
    return v_without + combined


def combine(*pieces):
    total = pieces[0].copy()
    for p in pieces[1:]:
        total = total + p
    return total
