"""
Schema for the three registry entities defined in docs/SPEC.md's
"核心資料模型" section:
  - DerivedFormula   a single paper-equation reproduction task
  - ConventionEntry  a single symbol/sign convention (module A)
  - BugWorkaround    a single cadabra2 known-bug workaround (module D)

verification/registry.json stores a list of DerivedFormula records (as
plain dicts matching DerivedFormula.to_dict()). This module is the schema
those dicts must conform to -- not a database layer, not an ORM. All
"persistence" in this project is git-tracked source + this JSON file
(see docs/SPEC.md "State 管理與持久化": there is no runtime database).

State machine (docs/SPEC.md "狀態模型與揭露策略"):
    not_started -> derived_single_route -> self_consistent -> paper_matched
                                                             -> regression_locked
    (any state) -> known_issue -> derived_single_route (once resolved)

Invariants enforced by DerivedFormula.__post_init__:
  - state in {self_consistent, paper_matched, regression_locked} requires
    derivation_route_b to be non-null (H2: two independent routes).
  - state == known_issue requires known_issue_note to be non-empty.
"""
from dataclasses import dataclass, asdict
from enum import Enum
from typing import Optional, List


class Provenance(str, Enum):
    DERIVED = 'derived'
    LITERATURE_UNVERIFIED = 'literature_unverified'
    SUPERSEDED_NON_CONFORMING = 'superseded_non_conforming'


class DerivedFormulaState(str, Enum):
    NOT_STARTED = 'not_started'
    DERIVED_SINGLE_ROUTE = 'derived_single_route'
    SELF_CONSISTENT = 'self_consistent'
    PAPER_MATCHED = 'paper_matched'
    REGRESSION_LOCKED = 'regression_locked'
    KNOWN_ISSUE = 'known_issue'


_STATES_REQUIRING_ROUTE_B = {
    DerivedFormulaState.SELF_CONSISTENT,
    DerivedFormulaState.PAPER_MATCHED,
    DerivedFormulaState.REGRESSION_LOCKED,
}


@dataclass
class DerivedFormula:
    id: str
    paper_ref: str
    ansatz: str
    provenance: Provenance
    derivation_route_a: str
    state: DerivedFormulaState
    derivation_route_b: Optional[str] = None
    last_verified_date: Optional[str] = None
    known_issue_note: Optional[str] = None

    def __post_init__(self):
        if not isinstance(self.provenance, Provenance):
            self.provenance = Provenance(self.provenance)
        if not isinstance(self.state, DerivedFormulaState):
            self.state = DerivedFormulaState(self.state)
        if self.state in _STATES_REQUIRING_ROUTE_B and not self.derivation_route_b:
            raise ValueError(
                f'DerivedFormula {self.id!r}: state={self.state.value} requires a non-null '
                f'derivation_route_b (docs/SPEC.md H2: two independent derivation routes).')
        if self.state == DerivedFormulaState.KNOWN_ISSUE and not self.known_issue_note:
            raise ValueError(
                f'DerivedFormula {self.id!r}: state=known_issue requires a non-empty '
                f'known_issue_note (what step it is stuck on, what was tried).')

    def to_dict(self):
        d = asdict(self)
        d['provenance'] = self.provenance.value
        d['state'] = self.state.value
        return d

    @classmethod
    def from_dict(cls, d):
        return cls(**d)


@dataclass
class ConventionEntry:
    symbol: str
    definition: str
    sign_note: str = ''

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, d):
        return cls(**d)


@dataclass
class BugWorkaround:
    bug_id: str
    symptom: str
    workaround_rule: str
    cadabra_version_tested: List[str]
    last_reverified_date: str
    upstream_issue_url: Optional[str] = None

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, d):
        return cls(**d)


# ---------------------------------------------------------------------
# The three cadabra2 bug workarounds already documented (and re-verified
# 2026-09-12) in cadabra/cadabra_utils.py's module docstring, expressed
# in this schema. This is the initial, factual content of the
# BugWorkaround registry -- not a placeholder.
# ---------------------------------------------------------------------
KNOWN_CADABRA_BUGS = [
    BugWorkaround(
        bug_id='accent_boundary_free_index_checker',
        symptom=(
            "Free/dummy-index checker inside Ex()-string-parsing and substitute()'s RHS-string "
            "parsing cannot see a contraction that crosses an Accent (\\delta{...}) boundary; a "
            "sum of two terms each carrying \\delta{...} gets flagged as 'Free indices in "
            "different terms in a sum do not match' even when every index is a valid dummy "
            "contraction."
        ),
        workaround_rule=(
            "D1: never let two additive terms that each carry a \\delta{...} factor appear "
            "together in one Ex()-parsed string or one substitute() RHS string. Build each "
            "additive piece as its own single-term Ex/substitute, then combine the built Ex "
            "trees with Python's `+` (use cadabra_utils.sub_copy/replace_term)."
        ),
        cadabra_version_tested=['2.4.5.4', '2.5.14'],
        last_reverified_date='2026-09-12',
        upstream_issue_url=None,
    ),
    BugWorkaround(
        bug_id='rename_dummies_inverse_metric_contraction',
        symptom=(
            "rename_dummies() raises \"No index set for index ... known\" on a valid contraction "
            "of a plain (undeclared-property) one-index tensor against an InverseMetric/Symmetric "
            "two-index tensor (e.g. g^{mu nu} A_mu A_nu), regardless of what property (if any) is "
            "declared on the one-index tensor. canonicalise() handles the identical expression "
            "without complaint."
        ),
        workaround_rule=(
            "D2: never call rename_dummies() directly. Always use canonicalise() instead (it "
            "performs equivalent dummy-renaming internally without hitting the bug). Enforced by "
            "cadabra/check_no_rename_dummies.py lint."
        ),
        cadabra_version_tested=['2.4.5.4', '2.5.14'],
        last_reverified_date='2026-09-12',
        upstream_issue_url=None,
    ),
    BugWorkaround(
        bug_id='substitute_dummy_name_collision',
        symptom=(
            "substitute() with a replacement RHS that introduces a dummy pair whose name (e.g. "
            "\\mu) is already in use as an unrelated dummy pair elsewhere in the target "
            "expression raises \"Failed to find dummy property for $\\mu$ while renaming "
            "dummies\", even though the two pairs are logically unrelated."
        ),
        workaround_rule=(
            "D3: never hand-pick index letters for a dummy pair introduced by a substitute() RHS. "
            "Always draw fresh names from cadabra_utils.fresh_indices(), which guarantees "
            "distinctness from every name already handed out in the kernel session."
        ),
        cadabra_version_tested=['2.4.5.4', '2.5.14'],
        last_reverified_date='2026-09-12',
        upstream_issue_url=None,
    ),
]
