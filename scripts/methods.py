"""One registry of which prediction columns should exist, and what makes them.

Without this, "is the benchmark up to date?" has no answer. `update.py` compared
the reference ids against the predictions files that EXIST, so a method whose
file had never been written was invisible: it reported "nothing to update" while
three rowan-cofold columns had never been produced at all. A missing file and a
file missing rows are different failures and both need naming.

The producer string must be the script that actually writes the file. Four PAE
columns were attributed to run_predictors.py and run_protenix_confidence.py when
interface_metrics.py writes them -- a registry that names the wrong producer is
worse than none, because it sends you to the wrong script with confidence. Check
with: grep -n 'write("<column>"' scripts/*.py

Adding a scorer means adding it here. That is the point: the registry is what
lets a tool say "this column is absent" rather than silently not looking for it.
"""

from __future__ import annotations

# column stem (predictions/<stem>_predictions.csv) -> what produces it
PRODUCERS: dict[str, str] = {
    "boltz-2": "run_predictors.py boltz",
    "boltz-2-iptm": "run_predictors.py boltz",
    "boltz-2-pae-mean": "interface_metrics.py",
    "boltz-2-pae-min": "interface_metrics.py",
    "boltz-2-pbind": "run_predictors.py boltz",
    "nesso-1": "run_predictors.py nesso",
    "nesso-1-pbind": "run_predictors.py nesso",
    "protenix-gpde": "run_protenix_confidence.py",
    "protenix-iptm": "run_protenix_confidence.py",
    "protenix-ligiptm": "run_protenix_confidence.py",
    "protenix-pae-mean": "interface_metrics.py",
    "protenix-pae-min": "interface_metrics.py",
    "protenix-plddt-binder": "run_protenix_confidence.py",
    "protenix-plddt-ligand": "run_protenix_confidence.py",
    "protenix-ptm": "run_protenix_confidence.py",
    "protenix-ranking": "run_protenix_confidence.py",
    "interface-ligburial": "run_interface.py",
    "interface-ncontacts": "run_interface.py",
    "clogp": "run_predictors.py baselines",
    "molecular-weight": "run_predictors.py baselines",
    "dtsfm-cosine": "run_dtsfm.py",
    "rowan-sqm": "run_rowan.py --submit (pose-based, needs a local co-fold)",
    "rowan-gnina": "run_rowan.py --submit (pose-based, needs a local co-fold)",
    "rowan-aevplig": "run_rowan.py --submit (pose-based, needs a local co-fold)",
    # hosted end to end: sequence + SMILES in, no GPU
    "rowan-cofold-affinity": "run_rowan_cofold.py --submit",
    "rowan-cofold-pbind": "run_rowan_cofold.py --submit",
    "rowan-cofold-iptm": "run_rowan_cofold.py --submit",
}

# What kind of method a column is. This lives here rather than in make_figure.py
# because it was a second, parallel dict there: adding the three rowan-cofold
# columns to the registry left the figure's map untouched and the figure crashed
# with a bare KeyError in CI. One registry, one place to add a method.
#   aff   dedicated affinity head        prox   structural proxy
#   base  trivial baseline               seq    sequence-native FM
#   phys  physics (SQM)                  dock   docking / ML rescore
#   iface interface geometry             combo  combination (LOO)
CATEGORY: dict[str, str] = {
    "boltz-2": "aff", "boltz-2-pbind": "aff",
    "nesso-1": "aff", "nesso-1-pbind": "aff",
    "rowan-cofold-affinity": "aff", "rowan-cofold-pbind": "aff",
    "boltz-2-iptm": "prox", "boltz-2-pae-min": "prox", "boltz-2-pae-mean": "prox",
    "rowan-cofold-iptm": "prox",
    "protenix-iptm": "prox", "protenix-ligiptm": "prox", "protenix-gpde": "prox",
    "protenix-ranking": "prox", "protenix-pae-min": "prox", "protenix-pae-mean": "prox",
    "protenix-ptm": "prox", "protenix-plddt-binder": "prox",
    "protenix-plddt-ligand": "prox",
    "clogp": "base", "molecular-weight": "base",
    "dtsfm-cosine": "seq",
    "rowan-sqm": "phys", "rowan-gnina": "dock", "rowan-aevplig": "dock",
    "interface-ncontacts": "iface", "interface-ligburial": "iface",
}

# Columns that are allowed to be absent: they are planned or optional, and their
# absence is not a gap to report every time anyone asks what needs running.
OPTIONAL: frozenset[str] = frozenset({
    "rowan-cofold-affinity", "rowan-cofold-pbind", "rowan-cofold-iptm",
})


def expected(include_optional: bool = True) -> list[str]:
    return sorted(PRODUCERS if include_optional else set(PRODUCERS) - OPTIONAL)
