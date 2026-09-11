"""
Tensor sector of generalized Proca theory: no-ghost function G_T,
propagation-normalisation function F_T, and tensor sound speed c_T^2.

    G_T = 2 [ G4 - 2 X G4_X - X A0 H G5_X ]
    F_T = 2 [ G4 - X G5_X Adot0 ]
    c_T^2 = F_T / G_T

NOTE ON PROVENANCE: the G4-only part of these expressions (G_T = F_T =
2 G4) follows from the tensor quadratic action of a G4(X) R
non-minimal coupling and is straightforward to re-derive (the -2X G4_X
correction to G_T comes from expanding the inverse-metric contractions
inside the G4_X[(nabla.A)^2 - nabla A nabla A] term to quadratic order
in the TT metric perturbation h_ij). The G5-dependent terms are quoted
from the generalized-Proca cosmological-perturbation literature
(De Felice et al., JCAP 06 (2016) 048, arXiv:1603.05806, and the
companion arXiv:1605.05066) and have NOT been independently
re-derived from the covariant action inside this project yet -- that
re-derivation (full quadratic tensor action from L5) is flagged as
outstanding work, consistent with the plan's own Phase 3.2 note that
the L6 contribution to c_T^2 is likewise still to be added.
"""
import sympy as sp

Mpl = sp.Symbol('M_pl', positive=True)


def G_T(G4, G4X, X, A0, H, G5X):
    return 2 * (G4 - 2 * X * G4X - X * A0 * H * G5X)


def F_T(G4, X, G5X, A0dot):
    return 2 * (G4 - X * G5X * A0dot)


def c_T2(G4, G4X, X, A0, H, G5X, A0dot):
    gt = G_T(G4, G4X, X, A0, H, G5X)
    ft = F_T(G4, X, G5X, A0dot)
    return sp.simplify(ft / gt)


def run_checks():
    results = []
    X, A0, H, G5X, A0dot = sp.symbols('X A0 H G5X A0dot', real=True)

    # ---- Check 1 (plan 3.c): GR limit. G4 = Mpl^2/2 constant (G4_X=0),
    # G5 = 0 (G5_X = 0) -> G_T = F_T = Mpl^2, c_T^2 = 1
    G4_gr = Mpl**2 / 2
    G4X_gr = 0
    gt = G_T(G4_gr, G4X_gr, X, A0, H, 0)
    ft = F_T(G4_gr, X, 0, A0dot)
    ct2 = sp.simplify(ft / gt)
    ok1 = sp.simplify(gt - Mpl**2) == 0 and sp.simplify(ft - Mpl**2) == 0 and sp.simplify(ct2 - 1) == 0
    results.append(('1. GR limit: G_T = F_T = Mpl^2, c_T^2 = 1', ok1, (gt, ft, ct2)))

    # ---- Check 2: G4=const (X-independent but not necessarily Mpl^2/2),
    # G4_X=0, G5=0 -> still G_T = F_T = 2 G4 for ANY constant, c_T^2 = 1
    G4c = sp.Symbol('G4c', positive=True)
    gt2 = G_T(G4c, 0, X, A0, H, 0)
    ft2 = F_T(G4c, X, 0, A0dot)
    ok2 = sp.simplify(gt2 - 2 * G4c) == 0 and sp.simplify(ft2 - 2 * G4c) == 0
    results.append(('2. Constant-G4, no G5: G_T = F_T = 2 G4 identically', ok2, (gt2, ft2)))

    # ---- Check 3: turning on G4_X alone (still G5=0) breaks G_T = F_T
    # (this is the hallmark of beyond-GR gravity: G_T picks up -4 X G4_X,
    # F_T does not) -- verifies the two formulas are NOT accidentally
    # identical once G4_X != 0.
    G4Xsym = sp.Symbol('G4X', real=True)
    gt3 = G_T(G4c, G4Xsym, X, A0, H, 0)
    ft3 = F_T(G4c, X, 0, A0dot)
    ok3 = sp.simplify(gt3 - ft3) == sp.simplify(-4 * X * G4Xsym)
    results.append(('3. G4_X != 0 splits G_T from F_T by exactly -4 X G4_X', ok3, sp.simplify(gt3 - ft3)))

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
    print('\nNOTE: G5-dependent terms in G_T/F_T are literature values, not yet')
    print('independently re-derived from the covariant action in this project.')
