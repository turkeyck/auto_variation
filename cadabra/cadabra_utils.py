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
"""
from cadabra2 import Ex, substitute


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
