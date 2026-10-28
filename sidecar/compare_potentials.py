#!/usr/bin/env python3
"""Run the gradient test under both potentials and report the difference.

The clearest sensitivity test this project has: the stand-in and the
published AWSEM tables are different enough that agreement between them
means the result is not an artefact of either.

    python3 sidecar/compare_potentials.py
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from scipy import stats
from glueforge.cif import read_structure
from glueforge.frustration import FrustrationConfig, compute_frustration
from glueforge.interface import assign_roles, label_contacts
from glueforge.potential import builtin_potential, load_awsem_json

ENTRIES = ["6TD3", "8BU5", "8BU7", "8BU9", "8BUA", "8BUJ", "8BUK", "8BUL",
           "6Q0W", "6Q0R", "6Q0V", "6UD7", "6H0F", "5FQD"]
pots = {"builtin": builtin_potential(),
        "awsem": load_awsem_json(Path(__file__).resolve().parent / 'data' / 'AWSEM_2015.json')}

print(f"{'entry':6s} {'rho builtin':>12s} {'rho awsem':>10s} {'p awsem':>9s}")
print("-" * 42)
res = {k: [] for k in pots}
for e in ENTRIES:
    s = read_structure(Path(__file__).resolve().parent.parent / 'cache' / 'raw' / f'{e}.cif')
    a = assign_roles(s)
    rs = [r for r in s.residues if r.chain in (a.anchor_chain, a.partner_chain)]
    line = {}
    for name, pot in pots.items():
        fr = compute_frustration(rs, pot, FrustrationConfig(n_decoys=800))
        lb = label_contacts(rs, fr.index_i, fr.index_j, a)
        m = lb.mask_interface
        r = stats.spearmanr(lb.contact_min_dist[m], fr.frustration[m])
        res[name].append(r.statistic)
        line[name] = (r.statistic, r.pvalue)
    print(f"{e:6s} {line['builtin'][0]:+12.3f} {line['awsem'][0]:+10.3f} "
          f"{line['awsem'][1]:9.4f}")

print()
for name in pots:
    v = np.array(res[name])
    pos = int((v > 0).sum())
    p = stats.binomtest(pos, len(v), 0.5, alternative='greater').pvalue
    print(f"{name:8s}  mean rho {v.mean():+.3f}   positive {pos}/{len(v)}   sign test p={p:.4g}")

b, a_ = np.array(res['builtin']), np.array(res['awsem'])
print(f"\ncorrelation between the two potentials' per-system rho: "
      f"{np.corrcoef(b, a_)[0,1]:.3f}")
