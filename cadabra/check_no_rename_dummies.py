"""
Lint check enforcing cadabra_utils.py Rule 2: never call rename_dummies()
directly (it has a known bug on plain-tensor/InverseMetric contractions
that canonicalise() does not share -- see cadabra_utils.py docstring).

Run as: python3 cadabra/check_no_rename_dummies.py
Exits non-zero and lists offending files/lines if any call is found.
"""
import re
import sys
from pathlib import Path

CALL_RE = re.compile(r'\brename_dummies\s*\(')
HERE = Path(__file__).parent
SELF_NAME = Path(__file__).name


def find_violations():
    violations = []
    for path in sorted(HERE.glob('*.py')):
        if path.name in (SELF_NAME, 'cadabra_utils.py'):
            continue
        for lineno, line in enumerate(path.read_text(encoding='utf-8').splitlines(), start=1):
            if CALL_RE.search(line):
                violations.append((path.name, lineno, line.strip()))
    return violations


if __name__ == '__main__':
    violations = find_violations()
    if violations:
        print('rename_dummies() is banned outside cadabra_utils.py (Rule 2) -- '
              'use canonicalise() instead:')
        for name, lineno, line in violations:
            print(f'  {name}:{lineno}: {line}')
        sys.exit(1)
    print('OK: no standalone rename_dummies() calls found.')
