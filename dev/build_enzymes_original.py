#!/usr/bin/env python3
"""Build an expanded Enzymes.csv containing every NEB-supplied restriction enzyme.

RECORD OF HOW THE FIRST VERSION OF inst/extdata/Enzymes.csv WAS MADE (2026-09-30). Kept for provenance, not used by the app.
Written by Claude Code at the authors' request. It reads enzymes from the REBASE snapshot bundled with Biopython (not from NEB's web site),
keeping those whose supplier list includes NEB. Later changes (DamDcm column dropped, Esp3I, PaqCI and the Common column added) were made
by hand; see docs/DATA.md and dev/make_enzy.R.

Usage:
    pip install biopython
    python build_enzymes.py Enzymes.csv Enzymes_expanded.csv [neb_methylation.csv]
"""
import csv
import sys

from Bio.Restriction import AllEnzymes

GOLDEN_GATE = {"BsaI", "BsmBI", "BbsI", "SapI", "BspQI", "BtgZI", "PaqCI"}


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    in_path, out_path = sys.argv[1], sys.argv[2]
    meth_path = sys.argv[3] if len(sys.argv) > 3 else None

    existing = {r["Enzyme"]: r for r in read_csv(in_path)}
    meth = {}
    if meth_path:
        meth = {r["Enzyme"]: r["DamDcm"].strip().upper() for r in read_csv(meth_path)}

    rows = dict(existing)
    unverified = []
    for enz in AllEnzymes:
        name = str(enz)
        if "N" not in enz.suppl or name in rows:
            continue
        dam = meth.get(name)
        if dam is None:
            dam = "FALSE"
            unverified.append(name)
        rows[name] = {
            "Enzyme": name,
            "Site": enz.site,
            "DamDcm": dam,
            "GG": "TRUE" if name in GOLDEN_GATE else "FALSE",
        }

    with open(out_path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["Enzyme", "Site", "DamDcm", "GG"])
        for name in sorted(rows, key=str.lower):
            r = rows[name]
            w.writerow([r["Enzyme"], r["Site"], r["DamDcm"], r["GG"]])

    with open("DamDcm_unverified.txt", "w") as f:
        f.write("\n".join(unverified) + "\n")

    print(f"{len(rows)} enzymes written to {out_path} "
          f"({len(existing)} kept from input, {len(rows) - len(existing)} added)")
    print(f"{len(unverified)} new enzymes have unverified DamDcm (see DamDcm_unverified.txt)")


if __name__ == "__main__":
    main()
