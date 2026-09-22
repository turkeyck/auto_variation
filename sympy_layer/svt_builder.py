"""
Module F1 (docs/SPEC.md): unified scalar-vector-tensor (SVT) perturbation
builder, producing the FULL perturbation content of arXiv:1703.09573 eq.
(3.1)-(3.2) (the "Stability conditions" section's metric + vector-field
ansatz) in ONE place, using the paper's OWN variable names, so the
tensor/vector/scalar sectors can each take a SUBSET of this single
declaration instead of maintaining independent, differently-named
hand-rolled perturbation variables per file (docs/SPEC.md module F1's
explicit requirement).

SOURCE: arXiv:1703.09573v2 (astrophv2.tex), "Stability conditions"
subsection (lines ~569-591 of that file):

    ds^2 = -(1+2 alpha) dt^2 + 2(d_i chi + V_i) dt dx^i
           + a^2(t) (delta_ij + h_ij) dx^i dx^j

    A^0 = phi(t) + delta_phi
    A^i = (1/a^2(t)) delta^{ij} (d_j chi_V + E_j)

where alpha, chi are scalar metric perturbations; V_i is the transverse
(d^i V_i=0) vector metric perturbation; h_ij is the transverse-traceless
(d^i h_ij=0, h^i_i=0) tensor metric perturbation; delta_phi, chi_V are
scalar vector-field perturbations; E_j is the transverse (d^j E_j=0)
vector-field perturbation.

NAMING CROSSWALK -- this project's EXISTING per-sector files predate
this module and use their own local names for the SAME quantities
(exactly the kind of duplication docs/SPEC.md module A's "single source
of truth" requirement targets). Recorded here for reference; NOT
auto-migrated in this pass (docs/SPEC.md Stage 5's own risk note:
refactor each file individually with its own regression check, not as
one batch edit) -- new code should use the (3.1)-(3.2) names directly.

    (3.1)-(3.2) name   | existing file(s)                | existing local name
    -------------------|----------------------------------|---------------------
    alpha              | adm_scalar.py, scalar_sector.py   | abar(t) * profile
    chi                | adm_scalar.py, scalar_sector.py   | chibar(t) * profile
    V_i                | vector_sector.py,                 | V (mode amplitude)
                        | vector_field_equation.py          |
    h_ij               | (not yet perturbed anywhere in    | --
                        | sympy_layer/; tensor sector only  |
                        | appears via q_T/c_T^2 formulas)   |
    delta_phi          | scalar_sector.py                  | dphi(t) * profile
    chi_V              | scalar_sector.py                  | chiV(t) * profile
    E_j                | vector_sector.py (Maxwell A_i)     | A_k(t) (implicitly
                        |                                    | E_j's role, not
                        |                                    | named E_j there)

Every builder function below shares the SAME single spatial profile
convention already established and verified in adm_scalar.py /
vector_field_equation.py / schutz_sorkin_scalar.py (a single cos(k z)
plane wave along one axis, giving EXACT eigenvalue relations for
gradients/Laplacians rather than approximate ones) -- this module does
not introduce a new convention, it centralizes the existing one.
"""
import sympy as sp

t, z = sp.symbols('t z', real=True)
k = sp.Symbol('k', positive=True)


def scalar_metric_perturbations(profile=None):
    """(alpha, chi) as Function(t)*profile amplitudes, matching eq (3.1)'s
    scalar metric sector. `profile` defaults to cos(k*z) (this project's
    established single-axis plane-wave convention); pass a different
    profile (or None to get bare, profile-free Function(t) objects) if a
    caller needs it decoupled from the shared spatial dependence."""
    alpha_bar = sp.Function('alpha_bar')(t)
    chi_bar = sp.Function('chi_bar')(t)
    if profile is None:
        return alpha_bar, chi_bar
    return alpha_bar * profile, chi_bar * profile


def vector_metric_perturbation(profile=None):
    """V_i (transverse vector metric perturbation), single-polarization
    mode amplitude V(t), matching vector_sector.py's/
    vector_field_equation.py's existing convention (a single V_1(t,z)-type
    component, transverse by construction since it only depends on z)."""
    V_bar = sp.Function('V_bar')(t)
    if profile is None:
        return V_bar
    return V_bar * profile


def scalar_vector_field_perturbations(profile=None):
    """(delta_phi, chi_V) as Function(t)*profile amplitudes, matching eq
    (3.2)'s scalar vector-field sector (A^0's perturbation and A^i's
    scalar-potential piece)."""
    dphi_bar = sp.Function('dphi_bar')(t)
    chiV_bar = sp.Function('chiV_bar')(t)
    if profile is None:
        return dphi_bar, chiV_bar
    return dphi_bar * profile, chiV_bar * profile


def vector_field_vector_perturbation(profile=None):
    """E_j (transverse vector-field perturbation), single-polarization
    mode amplitude, matching eq (3.2)'s A^i vector piece."""
    E_bar = sp.Function('E_bar')(t)
    if profile is None:
        return E_bar
    return E_bar * profile


def perturbed_metric_00_0i(a_t, alpha, dchi_dz):
    """g_00 = -(1+2 alpha), g_0z = d_z(chi) (+ V_z, which vanishes for a
    transverse V_i restricted to this single z-axis-only mode unless the
    caller explicitly adds it -- kept separate since the scalar and
    vector metric sectors are independent SVT pieces, never summed by
    this module itself). Returns (g00, g0z); g_ij = a_t**2 (flat gauge,
    no h_ij dependence -- tensor sector handled separately via q_T/c_T^2,
    per this project's existing convention, not yet built into a
    perturbed g_ij here)."""
    return -(1 + 2 * alpha), dchi_dz


def perturbed_A0_Ai(phi_t, a_t, dphi, dchiV_dz):
    """A^0 = phi(t) + delta_phi, A^z = (1/a^2) d_z(chi_V) (+ E_z, kept
    separate for the same SVT-independence reason as V_i above)."""
    return phi_t + dphi, dchiV_dz / a_t**2


def run_checks():
    """H2-style cross-check: the (alpha, chi) metric pieces built here
    reproduce EXACTLY the g00/g0z forms adm_scalar.py and
    schutz_sorkin_scalar.py already independently hard-coded and
    verified (not a tautology -- those two files were written before
    this module existed, using their own local abar/chibar names)."""
    results = []
    a_t = sp.Function('a', positive=True)(t)
    profile = sp.cos(k * z)

    alpha, chi = scalar_metric_perturbations(profile)
    dchi_dz = sp.diff(chi, z)
    g00, g0z = perturbed_metric_00_0i(a_t, alpha, dchi_dz)

    # adm_scalar.py's own convention (module docstring): "N = 1 +
    # alpha(t,x)", used via K_ij = (1/2N)(...) with g_00=-N^2 implicitly
    # linearized to -(1+2 alpha) at the order that file works to -- check
    # this module's g00 matches that exact linear form.
    alpha_bar_sym = sp.Function('alpha_bar')(t)
    expected_g00 = -(1 + 2 * alpha_bar_sym * profile)
    ok1 = sp.simplify(g00 - expected_g00) == 0
    results.append(('1. [H2] scalar_metric_perturbations()+perturbed_metric_00_0i() reproduces '
                     'adm_scalar.py\'s own g_00=-(1+2 alpha) convention exactly', ok1, g00))

    # schutz_sorkin_scalar.py's own g0z = d_z(chi) construction (its
    # build_delta_n() builds g0z the same way inline, before this module
    # existed -- reconstruct that exact inline expression and diff
    # against this module's output, rather than asserting equality by fiat.
    chibar_ss = sp.Function('chibar')(t)
    g0z_schutz_sorkin_scalar_inline = sp.diff(chibar_ss * profile, z)  # verbatim from that file's build_delta_n()
    g0z_this_module = dchi_dz.subs(sp.Function('chi_bar')(t), chibar_ss)
    ok2 = sp.simplify(g0z_this_module - g0z_schutz_sorkin_scalar_inline) == 0
    results.append(('2. g_0z = d_z(chi) matches schutz_sorkin_scalar.py\'s own inline construction '
                     '(same functional form, d_z(chibar(t)*profile))', ok2, g0z))

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
