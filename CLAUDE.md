# CLAUDE.md

Project-level guidance for Claude Code sessions working in this repo. Full
scope/architecture/acceptance criteria live in [docs/SPEC.md](docs/SPEC.md)
(the "規格整理 v 1.2.0" spec) — read it before making structural changes.
This file only holds the operational rules a fresh session needs before
its first tool call.

## Execution environment

This project spans two interpreters that cannot be merged into one:

- **`sympy_layer/` and everything under `conventions/`, `verification/`,
  `traceability/`, `tests/`** — run with the Windows Python launcher:
  `py -3 <script>.py` (confirmed: sympy 1.14.0). The system-default
  `python`/`python3` on this machine has neither `sympy` nor `cadabra2`
  installed — do not use them.
- **`cadabra/*.py`** — require `cadabra2` (2.5.14), which only exists
  inside WSL Ubuntu. Run via: `wsl.exe -- python3 cadabra/<file>.py`.
  WSL may show as "Stopped"; it starts automatically on first use.

`verification/run_all_checks.py` aggregates all `sympy_layer/` checks plus
the pure-Python `cadabra/check_no_rename_dummies.py` lint; it does **not**
run the cadabra2-dependent scripts (they're reported as explicitly
skipped, not silently omitted). Run it with `py -3 verification/run_all_checks.py`.

## Before writing or editing any `cadabra/*.py` script

Read `cadabra/cadabra_utils.py`'s module docstring first. It documents
three confirmed cadabra2 2.5.14 bugs (re-verified 2026-09-12; see
docs/SPEC.md §0 for the upstream version survey — no released or
in-development cadabra2 version fixes them) and three mandatory rules:

1. Never let two additive terms that each carry a `\delta{...}` (Accent)
   factor appear together in one `Ex()`-parsed string or one
   `substitute()` RHS string. Build each additive piece as its own
   single-term `Ex`/`substitute`, then combine the built trees with
   Python's `+` — use `cadabra_utils.sub_copy`/`replace_term`.
2. Never call `rename_dummies()` directly. Always use `canonicalise()`.
   Enforced by `cadabra/check_no_rename_dummies.py` — run it (or
   `verification/run_all_checks.py`, which includes it) before considering
   a new cadabra script done.
3. Never hand-pick index letters for a dummy pair introduced by a
   `substitute()` RHS. Always draw fresh names from
   `cadabra_utils.fresh_indices()`.

## Symbol conventions

`conventions/gp_conventions.py` is the single source of truth for shared
symbols (`t`, `k`, `Mpl`), the `X`/`F`/`Y` invariant definitions, the
signature `(-,+,+,+)`, the Riemann/Ricci sign convention, and the
`G2(X,F,Y)`..`G6(X)`/`g5(X)` theory registry. New `sympy_layer/` code
should `import` from it rather than redeclaring these symbols locally —
at least 6 existing files predate this module and still redeclare their
own copies; that duplication caused a real `G5,X` sign-error incident
(see `sympy_layer/proca_tensor.py`'s history). `cadabra/` scripts cannot
import this Python module (different interpreter) and must keep their
`Ex()` declarations consistent with it by hand — a known gap, not an
oversight.

## Testing

`py -3 verification/run_all_checks.py [--verbose]` is the one-shot check.
It prints known issues first (never buried under passing checks — a
direct fix for a past incident where a "6/7 passed" summary let the 1
open issue go unnoticed), then unexpected failures, then a per-module
summary, then the cadabra lint result. It also writes a machine-readable
`verification/last_run_summary.json` (gitignored — regenerate, don't
commit it).

`py -3 -m pytest tests/ -v` gives the same sympy_layer coverage at
per-check granularity (one pytest test per `run_checks()` entry across
every module), plus the cadabra lint, the H3 known-limits catalog
(`verification/known_limits.py`), and a static guard against
`schutz_sorkin_scalar.py` regressing back to the k-essence duality
shortcut. As of 2026-09-18 all 88 tests pass. cadabra's own symbolic
derivation scripts (`cadabra/*.py`, not `check_no_rename_dummies.py`)
are NOT pytest-wrapped — they need `cadabra2` (WSL-only) and are
print-based demonstrations, not `run_checks()`-returning functions; run
them individually via `wsl.exe -- python3 cadabra/<file>.py`.
