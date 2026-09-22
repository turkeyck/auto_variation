"""
Module H3 (docs/SPEC.md): catalog of this project's standing known-limit
tests (GR limit, standard Proca limit, f(R) limit, Maxwell limit).

This is a CATALOG/INDEX, not a re-implementation -- per docs/SPEC.md
Stage 6 ("不刪除原本檔案裡的極限測試，作為 regression 來源，
known_limits.py 是彙總入口"), the actual check logic stays exactly where
it already lives (each sympy_layer/*.py file's own run_checks()); this
module records WHICH check, in WHICH file, covers WHICH known limit, and
provides a function to verify none of them have silently disappeared
(e.g. via a rename that broke the substring match below) -- that
verification is what tests/test_known_limits.py exercises.

Horndeski-correspondence limit (also named in docs/SPEC.md H3) is NOT
yet covered by any existing check -- listed as `covered=False` rather
than omitted, so its absence is visible instead of silently unrecorded.
"""

KNOWN_LIMITS = [
    {
        'name': 'GR limit (vacuum, no matter)',
        'module': 'proca_tensor',
        'check_substring': 'GR limit: G_T = F_T = Mpl^2, c_T^2 = 1',
        'covered': True,
    },
    {
        'name': 'GR + Maxwell limit (vector sector)',
        'module': 'vector_sector',
        'check_substring': 'GR+Maxwell (check 3.d): c_V^2 = 1',
        'covered': True,
    },
    {
        'name': 'GR + standard Proca limit (scalar sector, no-ghost positivity)',
        'module': 'scalar_sector',
        'check_substring': 'K > 0 for all physical',
        'covered': True,
    },
    {
        'name': 'Massless limit of standard Proca (m^2->0): longitudinal mode stops propagating',
        'module': 'scalar_sector',
        'check_substring': 'Massless limit',
        'covered': True,
    },
    {
        'name': 'GR + fluid-only limit (m^2->0, Proca sector off): canonical-scalar normalisation',
        'module': 'scalar_plus_fluid_sector',
        'check_substring': 'GR+fluid-only limit',
        'covered': True,
    },
    {
        'name': 'f(R)=R exact GR limit (covariant route)',
        'module': 'fR_gravity',
        'check_substring': "f(R)=R gives F=f'(R)=1",
        'covered': True,
    },
    {
        'name': 'f(R) GR limit (mini-superspace route, vacuum Friedmann)',
        'module': 'fR_gravity',
        'check_substring': 'GR limit (f=R)',
        'covered': True,
    },
    {
        'name': 'Standard Proca limit of the full GP background (be1-be3, G3=G4,X=G5=G6=g5=0)',
        'module': None,
        'check_substring': None,
        'covered': False,  # docs/SPEC.md 階段三 scope, not yet built
    },
    {
        'name': 'Horndeski-correspondence limit (L=F+L4+L6, constant G4,G6)',
        'module': None,
        'check_substring': None,
        'covered': False,  # not yet built
    },
]


def covered_limits():
    return [entry for entry in KNOWN_LIMITS if entry['covered']]


def uncovered_limits():
    return [entry for entry in KNOWN_LIMITS if not entry['covered']]
