"""
Module A (docs/SPEC.md): single source of truth for symbol/sign conventions
used across this project (near-term target arXiv:1703.09573 generalized
Proca theory; background cross-reference arXiv:1603.05806).

docs/SPEC.md's module A requirement: this must be the ONLY place these
conventions are declared. Every sympy_layer/ script should import shared
symbols and theory registrations from here instead of redeclaring its own
copies of t, k, Mpl, X, G2..G6, etc. (docs/SPEC.md "開發應注意重點" item 1
names this exact duplication -- at least 6 existing files each redeclare
t/k/Mpl/X independently -- as the root cause of the G5,X sign-error
incident recorded in sympy_layer/proca_tensor.py's history.)

cadabra/ scripts run under a separate interpreter (cadabra2, only
reachable via WSL) and cannot `import` this SymPy module directly. Until
a cadabra-side mirror is built (not in Stage 0's scope per docs/SPEC.md),
cadabra scripts must keep their own Ex() declarations consistent with the
definitions documented here BY HAND -- this file remains the textual
source of truth for both layers even though only one of them can import it.
This is a known gap, not an oversight; tracked as follow-on work.

SIGNATURE: mostly-plus, (-,+,+,+). Matches
g = diag(-N(t)^2, a(t)^2, a(t)^2, a(t)^2), used throughout sympy_layer/
(proca_minisuperspace.py, adm_scalar.py, fR_gravity.py, ...). Do not
introduce a (+,-,-,-) convention anywhere in this project.

RIEMANN / RICCI SIGN CONVENTION: as implemented by
sympy_layer/tensor_utils.py's riemann_tensor()/ricci_tensor():
    R^rho_{sigma mu nu} = d_mu Gamma^rho_{nu sigma} - d_nu Gamma^rho_{mu sigma}
                         + Gamma^rho_{mu lambda} Gamma^lambda_{nu sigma}
                         - Gamma^rho_{nu lambda} Gamma^lambda_{mu sigma}
    R_{mu nu} = R^rho_{mu rho nu}
This is the ONE Riemann/Ricci sign convention for the whole project;
any covariant (cadabra) computation of curvature must agree with it.

BACKGROUND VECTOR FIELD ANSATZ: A^mu = (A0(t), 0, 0, 0). A0(t) is declared
`real=True`, NOT `positive=True` -- the papers do not restrict its sign,
and several existing files (sympy_layer/proca_minisuperspace.py) already
follow this; do not tighten the assumption elsewhere.
"""
import sympy as sp

# ---------------------------------------------------------------------
# Shared bare symbols. Import these instead of re-declaring
# `t = sp.Symbol('t', real=True)` etc. locally.
# ---------------------------------------------------------------------
t = sp.Symbol('t', real=True)
k = sp.Symbol('k', positive=True)
Mpl = sp.Symbol('M_pl', positive=True)

# ---------------------------------------------------------------------
# Invariant definitions -- symbolic, ansatz-independent. These are
# DEFINITIONS (the single textual source new code must match), not yet
# generic covariant-construction functions: module C (invariants.py,
# scheduled for a later stage per docs/SPEC.md, not built yet) owns
# building X/F/Y covariantly from g, A, nabla for an arbitrary ansatz.
# Until module C exists, any script that needs a concrete value of X/F/Y
# under a specific ansatz must derive it from these definitions in place
# and cite this docstring, rather than inventing an independent definition.
# ---------------------------------------------------------------------
INVARIANT_DEFINITIONS = {
    'X': 'X = -1/2 * A_mu * A^mu',
    'F': 'F = -1/4 * F_{mu nu} * F^{mu nu},  F_{mu nu} = nabla_mu A_nu - nabla_nu A_mu',
    'Y': 'Y = A^mu * A^nu * F_mu^{alpha} * F_{nu alpha}',
}

Xsym, Fsym, Ysym = sp.symbols('X F Y', real=True)

# ---------------------------------------------------------------------
# G_i(X,F,Y) theory registry -- arXiv:1703.09573 action:
#   S = int d^4x sqrt(-g) (L2 + L3 + L4 + L5 + L6 + L_M)
# Symbolic function registration only, NOT tied to any ansatz. Concrete
# X-derivatives (G_i,X etc.) are produced ad hoc per script via sp.diff()
# on these; this registry exists so every script asks for "the G2 of this
# project" from one shared place instead of re-declaring
# sp.Function('G2') independently.
# ---------------------------------------------------------------------
G2 = sp.Function('G2')(Xsym, Fsym, Ysym)   # L2 = G2(X,F,Y)
G3 = sp.Function('G3')(Xsym)               # L3 = G3(X) nabla_mu A^mu
G4 = sp.Function('G4')(Xsym)               # L4 = G4(X) R + G4,X[(nabla.A)^2 - nabla_mu A_nu nabla^nu A^mu]
G5 = sp.Function('G5')(Xsym)               # L5 = G5(X) G_{mu nu} nabla^mu A^nu - (1/6) G5,X[...]
g5 = sp.Function('g5')(Xsym)               # L5 extra term: -g5(X) Ftilde^{alpha mu} Ftilde^beta_mu nabla_alpha A_beta
G6 = sp.Function('G6')(Xsym)               # L6 = G6(X) L^{mu nu al be} nabla_mu A_nu nabla_al A_be + (1/2) G6,X Ftilde Ftilde nabla A nabla A

THEORY_REGISTRY = {
    'L2': G2, 'L3': G3, 'L4': G4, 'L5': G5, 'g5': g5, 'L6': G6,
}


def register_convention(symbol, definition, sign_note=''):
    """Build a ConventionEntry-shaped dict (see verification/registry_schema.py
    for the authoritative field set). Call sites should use this helper
    rather than constructing the dict by hand, so the two stay in sync."""
    return {'symbol': symbol, 'definition': definition, 'sign_note': sign_note}


CONVENTION_REGISTRY = [
    register_convention(
        'signature', '(-,+,+,+) mostly-plus',
        'Matches g=diag(-N^2,a^2,a^2,a^2) already used throughout sympy_layer/.'),
    register_convention(
        'riemann_ricci', 'R^rho_{sigma mu nu} per tensor_utils.riemann_tensor(); R_{mu nu}=R^rho_{mu rho nu}',
        'The one Riemann/Ricci sign convention for the whole project; see module docstring for the exact formula.'),
    register_convention(
        'X', INVARIANT_DEFINITIONS['X'],
        'X = A0^2/(2N^2) > 0 on the background ansatz A^mu=(A0(t),0,0,0); A0(t) itself carries no assumed sign.'),
    register_convention(
        'F', INVARIANT_DEFINITIONS['F'],
        'Vanishes identically on the pure-temporal background ansatz A^mu=(A0(t),0,0,0).'),
    register_convention(
        'Y', INVARIANT_DEFINITIONS['Y'],
        'Vanishes identically on the pure-temporal background ansatz A^mu=(A0(t),0,0,0).'),
    register_convention(
        'A0_sign', 'A0(t) declared sp.Function(\'A0\')(t) with no positivity assumption.',
        'Do not add positive=True anywhere; the papers do not restrict A0(t)\'s sign.'),
]
