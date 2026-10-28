#!/usr/bin/env python3
"""Pull structures from RCSB into cache/raw.

The sandbox this was built in had no outbound network, so the fourteen
structures already in the cache were attached by hand. This script is how the
cache gets filled anywhere else -- run it once on a networked machine and the
rest of the pipeline works unchanged.

    python3 sidecar/fetch.py                  # the default benchmark set
    python3 sidecar/fetch.py 6TD3 5FQD        # named entries
    python3 sidecar/fetch.py --list           # show the default set and exit

Uses only the standard library, so there is nothing to install first.
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import sys
import urllib.error
import urllib.request
from pathlib import Path

RAW = Path(__file__).resolve().parent.parent / "cache" / "raw"
URL = "https://files.rcsb.org/download/{}.cif.gz"

#: The benchmark set, grouped by biological system.
DEFAULT_SET: dict[str, list[str]] = {
    "CDK12-cyclin K x DDB1": [
        "6TD3",  # CR8, the original
        "8BU5", "8BU7", "8BU9", "8BUA", "8BUJ", "8BUK", "8BUL",  # the SAR series
    ],
    "DCAF15 x RBM39": [
        "6Q0W",  # indisulam
        "6Q0V",  # tasisulam
        "6Q0R",  # E7820
        "6UD7",  # indisulam, the other crystal form
    ],
    "CRBN x C2H2 zinc finger": [
        "6H0F",  # pomalidomide + IKZF1 ZF2
        "5FQD",  # lenalidomide + CK1a
    ],
}

#: Not in the default set, but the obvious fourth system: a PPI stabiliser
#: outside the DDB1/CRL4 family, which removes the shared-machinery confound.
SUGGESTED_EXTRA = {
    "14-3-3 x client (fusicoccin-class)": [
        "see the 14-3-3 stabiliser literature; entries vary by client",
    ],
}


def fetch_one(entry: str, dest_dir: Path, force: bool = False) -> bool:
    entry = entry.upper()
    dest = dest_dir / f"{entry}.cif"
    if dest.exists() and not force:
        print(f"  {entry}: already present ({dest.stat().st_size // 1024} kB)")
        return True

    url = URL.format(entry)
    tmp = dest_dir / f".{entry}.cif.gz"
    try:
        with urllib.request.urlopen(url, timeout=90) as resp, open(tmp, "wb") as fh:
            shutil.copyfileobj(resp, fh)
    except urllib.error.HTTPError as exc:
        print(f"  {entry}: HTTP {exc.code} -- is that a real entry id?")
        tmp.unlink(missing_ok=True)
        return False
    except Exception as exc:  # noqa: BLE001
        print(f"  {entry}: {type(exc).__name__}: {exc}")
        tmp.unlink(missing_ok=True)
        return False

    try:
        with gzip.open(tmp, "rb") as src, open(dest, "wb") as out:
            shutil.copyfileobj(src, out)
    finally:
        tmp.unlink(missing_ok=True)

    print(f"  {entry}: {dest.stat().st_size // 1024} kB")
    return True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("entries", nargs="*", help="PDB ids; default is the benchmark set")
    ap.add_argument("--force", action="store_true", help="re-download existing files")
    ap.add_argument("--list", action="store_true", help="print the default set and exit")
    ap.add_argument("--out", default=str(RAW))
    args = ap.parse_args()

    if args.list:
        for system, ids in DEFAULT_SET.items():
            print(f"{system}:")
            print("  " + " ".join(ids))
        print("\nnot included, worth adding:")
        for system, note in SUGGESTED_EXTRA.items():
            print(f"  {system}: {note[0]}")
        return 0

    dest_dir = Path(args.out)
    dest_dir.mkdir(parents=True, exist_ok=True)

    entries = args.entries or [e for ids in DEFAULT_SET.values() for e in ids]
    print(f"fetching {len(entries)} entries into {dest_dir}")

    ok = sum(fetch_one(e, dest_dir, args.force) for e in entries)
    print(f"\n{ok}/{len(entries)} available")
    if ok:
        print("next:  python3 sidecar/run_oracle.py")
        print("       python3 sidecar/export_web.py && node build.js")
    return 0 if ok == len(entries) else 1


if __name__ == "__main__":
    raise SystemExit(main())
