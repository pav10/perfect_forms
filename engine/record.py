"""JSON records (plan, section 1.2): one per field, with provenance."""
import os
import platform
import subprocess
import time
from collections import Counter
from fractions import Fraction

from . import __version__, certify, extreme
from .field import Order, member
from .pari import pari
from .unary import Voronoi

PRECISION_BITS = 256          # precision of the floating point embeddings (heuristics only)


def _git_hash():
    here = os.path.dirname(os.path.abspath(__file__))
    try:
        h = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=here, capture_output=True, text=True,
                           timeout=10).stdout.strip()
        dirty = subprocess.run(['git', 'status', '--porcelain', '--', 'engine'], cwd=here,
                               capture_output=True, text=True, timeout=10).stdout.strip()
        return (h or 'unknown') + ('+dirty' if dirty else '')
    except Exception:
        return 'unknown'


_GIT = None


def provenance(seconds):
    global _GIT
    if _GIT is None:
        _GIT = _git_hash()
    try:
        import cypari2
        from importlib.metadata import version
        cyv = version('cypari2')
    except Exception:
        cyv = 'unknown'
    return {
        'code': 'perfect_forms/engine %s' % __version__,
        'git': _GIT,
        'pari': '.'.join(str(int(t)) for t in pari('version()')[:3]),
        'cypari2': cyv,
        'python': platform.python_version(),
        'precision_bits': PRECISION_BITS,
        'seed': 1,
        'seconds': round(seconds, 3),
    }


def _q(x):
    """Exact rational as a string."""
    x = Fraction(x)
    return str(x.numerator) if x.denominator == 1 else '%d/%d' % (x.numerator, x.denominator)


def compute(pol, label=None, forms=True, relmin=False, certify_bnf=False):
    """Run Voronoi's algorithm for the maximal order of Q[x]/(pol) and return the record."""
    t0 = time.time()
    O = Order(pol, label=label)
    V = Voronoi(O).run()
    cert = certify.certificates(O, V)
    ext = V.extreme_square_classes()
    gam = extreme.summary(O, V)
    f = V.f_vector()
    m = O.m
    rec = {
        'label': label,
        'pol': O.polredabs(),
        'pol_input': str(O.pol),
        'degree': m,
        'signature': list(O.signature),
        'disc': O.disc,
        'order_index': O.index,
        'class_number': O.class_number(),
        'regulator': O.regulator(),
        'fundamental_units': [O.poly_str(u) for u in O.reduced_units()],
        'integral_basis': [str(pari.lift(a)) for a in O.basis_alg],
        'n': f[-1],
        'e': f[m - 2] if m >= 2 else 0,
        'mu': f[0],
        'f_vector': f,
        'degrees': dict(sorted(Counter(V.degrees()).items())),
        'facet_sizes': dict(sorted(Counter(V.facet_sizes()).items())),
        'extreme_square_norms': sorted(abs(O.norm(x)) for x in ext.values()),
        'gamma_pow': _q(gam['gamma_pow']),
        'gamma': gam['gamma'],
        'gamma_argmax': gam['argmax'],
        'n_extreme': gam['n_extreme'],
        'n_weakly_eutactic': gam['n_weakly_eutactic'],
        'leibak_ratio': _q(gam['leibak_ratio']) if 'leibak_ratio' in gam else None,
        'certificates': cert,
    }
    if certify_bnf:
        rec['bnfcertify'] = O.certify_bnf()
    if relmin:
        from .relmin import RelativeMinima
        R = RelativeMinima(O).run()
        rec['ell'] = R.count()
        rec['certificates']['relmin'] = certify.relmin_certificates(O, V, R)
        rec['certificates']['certified'] = rec['certificates']['certified'] and rec['certificates']['relmin']['ok']
    if forms:
        out = []
        for i, P in enumerate(V.classes):
            out.append({
                'alpha': O.poly_str(P.alpha),
                'alpha_coords': [_q(t) for t in P.alpha],
                'norm': _q(P.norm),
                'degree': P.degree,
                'eutaxy': gam['eutaxy'][i],
                'minvecs': [list(x) for x in P.minvecs],
                'extreme': [k in P.vertices for k in range(len(P.minvecs))],
                'minvec_norms': [abs(O.norm(x)) for x in P.minvecs],
                'neighbours': [j for j, _ in P.neighbours],
            })
        rec['forms'] = out
    rec['status'] = 'CERTIFIED' if rec['certificates']['certified'] else 'QUARANTINE'
    rec['stats'] = V.stats
    rec['provenance'] = provenance(time.time() - t0)
    return rec
