"""
H1 ground-truth transcription: arXiv:1703.09573v2 (De Felice, Heisenberg,
Tsujikawa, "Observational constraints on generalized Proca theories",
Phys. Rev. D 95, 123540), Section II background equations.

SOURCE: downloaded directly from the arXiv e-print (`arxiv.org/e-print/
1703.09573v2`), single source file `astrophv2.tex`. NOT OCR'd from the
PDF, NOT transcribed from memory. Line numbers below refer to that file
as retrieved 2026-09-17.

Action (astrophv2.tex lines 198-254, \\label{LagProca}, L2, L3, L4, L5, L6):
    S = int d^4x sqrt(-g) (L2+L3+L4+L5+L6 + L_M)
    L2 = G2(X,F,Y)
    L3 = G3(X) nabla_mu A^mu
    L4 = G4(X) R + G4,X(X) [(nabla_mu A^mu)^2 - nabla_rho A_sigma nabla^sigma A^rho]
    L5 = G5(X) G_{mu nu} nabla^mu A^nu
         - (1/6) G5,X(X) [(nabla_mu A^mu)^3 - 3(nabla_mu A^mu)(nabla_rho A_sigma nabla^sigma A^rho)
                          + 2 nabla_rho A_sigma nabla^gamma A^rho nabla^sigma A_gamma]
         - g5(X) Ftilde^{alpha mu} Ftilde^beta_mu nabla_alpha A_beta
    L6 = G6(X) L^{mu nu alpha beta} nabla_mu A_nu nabla_alpha A_beta
         + (1/2) G6,X(X) Ftilde^{alpha beta} Ftilde^{mu nu} nabla_alpha A_mu nabla_beta A_nu
    X = -1/2 A_mu A^mu,  F = -1/4 F_{mu nu} F^{mu nu},  Y = A^mu A^nu F_mu^alpha F_{nu alpha}
    (lines 213-246)

NOTE ON L4/L5 vs. arXiv:1603.05806 (the background-theory companion paper
this project's sympy_layer/proca_minisuperspace.py was originally written
against): 1603.05806's L4, L5 carry free "intrinsic vector mode" constants
c2 (in L4) and d2 (in L5) (jcapresub1.tex lines 264-282). At c2=0, d2=0,
1603.05806's L4/L5 reduce EXACTLY (term-for-term, confirmed by direct
comparison of the two LaTeX sources) to 1703.09573's L4/L5 above minus the
g5(X) term. This is why proca_minisuperspace.py's pre-existing L4 (c2=0
already hard-coded) and L5 (with symbolic d2, whose background contribution
was independently proven d2-INDEPENDENT by that file's own check 7) are not
simply "the wrong theory" -- they already coincide with 1703.09573's L4/L5
on background contributions. What was missing, and what this module and
sympy_layer/proca_background_1703_09573.py add, is: (a) an explicit,
literal transcription of 1703.09573's own background equations (be1-be3)
to compare against, since the pre-existing file only ever checked internal
self-consistency (regression, Bianchi-identity reduction, d2-independence)
and NEVER did the H1 "paper is ground truth" comparison this module now
provides; (b) an explicit verification that G2's extra F,Y dependence and
the g5(X), G6(X) terms truly drop out of the background equations, rather
than assuming the paper's own claim.

BACKGROUND ANSATZ (astrophv2.tex lines 268-277): flat FLRW,
ds^2 = -dt^2 + a(t)^2 delta_ij dx^i dx^j (lapse already gauge-fixed to
N=1 in the paper's own presentation -- unlike sympy_layer/
proca_minisuperspace.py, which keeps a general N(t) and gauge-fixes to
N=1 only after taking the N-variation, exactly so that "vary w.r.t. N"
recovers the same Hamiltonian-constraint content the paper calls be1).
A^mu = (phi(t), 0, 0, 0).

KEY STRUCTURAL FACT, quoted verbatim (astrophv2.tex lines 323-326):
    "The functions g_5, G_6 and the additional dependence of F and Y in
    the function G_2, which correspond to intrinsic vector modes, do not
    contribute to the background equations of motion as expected."
This is why the background transcriptions below use G2(X) (not G2(X,F,Y))
and omit g5, G6 entirely -- not because this module assumes it, but
because it is what the paper's own text states. F and Y both vanish
identically for this purely-temporal-A, purely-time-dependent-a ansatz
(F_{mu nu} = nabla_mu A_nu - nabla_nu A_mu is identically zero for
A_mu=(phi(t),0,0,0) on FLRW -- every Christoffel-symbol term that could
make it nonzero, Gamma^0_{0i} and the symmetric Gamma^0_{ij} pair,
cancels in the antisymmetrization), which is the more fundamental
reason G2's F,Y-dependence and the g5/G6 terms (which are all built from
F_{mu nu} or its dual) cannot contribute here regardless of the paper's
own remark -- both facts are checked independently in
sympy_layer/proca_background_1703_09573.py rather than taken on faith.

BACKGROUND FIELD EQUATIONS (astrophv2.tex lines 292-321):
  Varying S w.r.t. g_{mu nu} gives the modified Einstein equations
  (\\label{be1}, \\label{be2}):

    be1:  G2 - G2,X phi^2 - 3 G3,X H phi^3 + 6 G4 H^2
          - 6(2 G4,X + G4,XX phi^2) H^2 phi^2
          + G5,XX H^3 phi^5 + 5 G5,X H^3 phi^3  =  rho_M

    be2:  G2 - phidot phi^2 G3,X + 2 G4 (3H^2 + 2 Hdot)
          - 2 G4,X phi (3H^2 phi + 2 H phidot + 2 Hdot phi)
          - 4 G4,XX H phidot phi^3 + G5,XX H^2 phidot phi^4
          + G5,X H phi^2 (2 Hdot phi + 2 H^2 phi + 3 H phidot)  =  -P_M

  Varying S w.r.t. phi gives (\\label{be3}):

    be3:  phi ( G2,X + 3 G3,X H phi + 6 G4,X H^2 + 6 G4,XX H^2 phi^2
                - 3 G5,X H^3 phi - G5,XX H^3 phi^3 )  =  0

  where H = adot/a, a dot denotes d/dt, and every G_i, G_i,X, G_i,XX is
  evaluated at the background value X = phi(t)^2/2 (from X=-1/2 A_mu A^mu
  at N=1, A^mu=(phi,0,0,0)).

These three are your original spec's "(2.11)-(2.13)" -- the paper's own
compiled section/equation numbering was not independently re-verified
against a compiled PDF in this pass (unlike arXiv:1603.05806's numbering,
which the project's schutz_sorkin_vector.py history records as having
been checked via a compiled .aux file after an earlier mislabeling
incident); the \\label{be1}/\\label{be2}/\\label{be3} anchors are used
here as the authoritative identifiers instead, with "as referenced by the
professor's spec as eq (2.11)-(2.13)" recorded as an assumption, not a
re-verified fact. Flagged in docs/SPEC.md "風險與未決事項" item 2.
"""
import sympy as sp


def background_field_equations(a_t, phi_t, G2, G3, G4, G5, Xsym, t):
    """Build (be1_LHS, be2_LHS, be3_LHS) exactly as transcribed above, as
    functions of a_t=a(t), phi_t=phi(t) (the Function(t) objects used by
    the caller), G2..G5 (sp.Function('Gi')(Xsym) objects; only their
    X-dependence is used -- G2's F,Y-dependence is irrelevant here per
    the paper's own statement quoted in this module's docstring), the
    caller's own X symbol (deliberately taken as a parameter rather than
    declared in this module: sympy_layer/proca_minisuperspace.py declares
    its X symbol with `positive=True` while a module-local declaration
    here would default to a plain/real symbol -- two Symbol('X', ...)
    objects with different assumptions do NOT compare equal to sympy, so
    Xsym must be threaded through from the caller for .subs() to work),
    and the time variable t.

    Returns three expressions with X already substituted to phi_t**2/2
    (background value) in G_i and its X-derivatives; RHS (rho_M, -P_M, 0)
    is NOT included -- callers compare LHS against their own derived
    gravity+Proca-sector result, which is the natural "= source term"
    split for a minimally-coupled matter action (see module docstring).
    """
    H = sp.diff(a_t, t) / a_t
    Hdot = sp.diff(H, t)
    phi = phi_t
    phidot = sp.diff(phi_t, t)
    Xbg = phi_t**2 / 2

    def at_bg(expr):
        return expr.subs(Xsym, Xbg)

    G2_ = at_bg(G2)
    G2X = at_bg(sp.diff(G2, Xsym))
    G3X = at_bg(sp.diff(G3, Xsym))
    G4_ = at_bg(G4)
    G4X = at_bg(sp.diff(G4, Xsym))
    G4XX = at_bg(sp.diff(G4, Xsym, 2))
    G5X = at_bg(sp.diff(G5, Xsym))
    G5XX = at_bg(sp.diff(G5, Xsym, 2))

    be1 = (G2_ - G2X * phi**2 - 3 * G3X * H * phi**3 + 6 * G4_ * H**2
           - 6 * (2 * G4X + G4XX * phi**2) * H**2 * phi**2
           + G5XX * H**3 * phi**5 + 5 * G5X * H**3 * phi**3)

    be2 = (G2_ - phidot * phi**2 * G3X + 2 * G4_ * (3 * H**2 + 2 * Hdot)
           - 2 * G4X * phi * (3 * H**2 * phi + 2 * H * phidot + 2 * Hdot * phi)
           - 4 * G4XX * H * phidot * phi**3 + G5XX * H**2 * phidot * phi**4
           + G5X * H * phi**2 * (2 * Hdot * phi + 2 * H**2 * phi + 3 * H * phidot))

    be3 = phi * (G2X + 3 * G3X * H * phi + 6 * G4X * H**2 + 6 * G4XX * H**2 * phi**2
                 - 3 * G5X * H**3 * phi - G5XX * H**3 * phi**3)

    return be1, be2, be3
