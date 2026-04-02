"""
Command-line interface for terrahedral.

    terrahedral info structure.cif --metal Cu1
    terrahedral list-metals --pdb 6GTL
    terrahedral fetch --pdb 6GTL
    terrahedral fetch --cod 4070506
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _add_source_args(parser: argparse.ArgumentParser) -> None:
    """Add mutually exclusive CIF source arguments: file path, --pdb, or --cod."""
    group = parser.add_argument_group("CIF source (provide one)")
    group.add_argument("cif", nargs="?", default=None,
                       help="Path to a local .cif file")
    group.add_argument("--pdb", default=None, metavar="ID",
                       help="Fetch mmCIF from RCSB PDB (e.g. 6GTL)")
    group.add_argument("--cod", default=None, metavar="ID",
                       help="Fetch CIF from Crystallography Open Database (e.g. 4070506)")
    group.add_argument("--cache-dir", default=None, metavar="DIR",
                       help="Cache directory for fetched files (default: ~/.cache/terrahedral)")
    group.add_argument("--no-cache", action="store_true",
                       help="Re-download even if file is cached")


def _resolve_cif(args) -> str:
    """
    Resolve the CIF file path from args.  Downloads if --pdb or --cod given.
    Also sets args.mmcif = True when --pdb is used.
    """
    sources = sum(x is not None for x in [args.cif, args.pdb, args.cod])
    if sources == 0:
        print("Error: provide a CIF file path, --pdb ID, or --cod ID", file=sys.stderr)
        sys.exit(1)
    if sources > 1:
        print("Error: provide only one of: CIF file path, --pdb, --cod", file=sys.stderr)
        sys.exit(1)

    if args.pdb:
        from terrahedral.fetch import fetch_pdb
        path = fetch_pdb(
            args.pdb,
            cache_dir=args.cache_dir,
            force=getattr(args, "no_cache", False),
        )
        # PDB files are always mmCIF
        if hasattr(args, "mmcif"):
            args.mmcif = True
        print(f"Using {path}")
        return str(path)

    if args.cod:
        from terrahedral.fetch import fetch_cod
        path = fetch_cod(
            args.cod,
            cache_dir=args.cache_dir,
            force=getattr(args, "no_cache", False),
        )
        print(f"Using {path}")
        return str(path)

    return args.cif


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="terrahedral",
        description="CIF -> coordination geometry analysis",
    )
    sub = parser.add_subparsers(dest="command")

    # ── info ──
    p_info = sub.add_parser("info", help="Print coordination site summary")
    _add_source_args(p_info)
    p_info.add_argument("--metal", "-m", required=True, help="Metal atom label")
    p_info.add_argument("--mmcif", action="store_true")
    p_info.add_argument("--asym", default="B")
    p_info.add_argument("--seq", default=".")

    # ── list-metals ──
    p_list = sub.add_parser("list-metals", help="List metal atoms in a CIF")
    _add_source_args(p_list)
    p_list.add_argument("--mmcif", action="store_true")

    # ── fetch ──
    p_fetch = sub.add_parser("fetch", help="Download a CIF file without processing")
    p_fetch.add_argument("--pdb", default=None, metavar="ID",
                         help="Fetch mmCIF from RCSB PDB")
    p_fetch.add_argument("--cod", default=None, metavar="ID",
                         help="Fetch CIF from COD")
    p_fetch.add_argument("--cache-dir", default=None, metavar="DIR",
                         help="Cache directory (default: ~/.cache/terrahedral)")
    p_fetch.add_argument("--no-cache", action="store_true",
                         help="Re-download even if cached")
    p_fetch.add_argument("--output", "-o", default=None,
                         help="Copy fetched file to this path")

    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return 0

    if args.command == "fetch":
        return _cmd_fetch(args)
    elif args.command == "info":
        return _cmd_info(args)
    elif args.command == "list-metals":
        return _cmd_list_metals(args)

    return 0


def _load_site(args):
    from terrahedral import CoordinationSite
    filepath = _resolve_cif(args)
    if args.mmcif:
        return CoordinationSite.from_mmcif(
            filepath, metal_asym=args.asym, metal_seq=args.seq,
            site_label=args.metal,
        )
    else:
        return CoordinationSite.from_cif(filepath, metal=args.metal)


def _cmd_info(args) -> int:
    site = _load_site(args)
    site.classify_geometry()
    print(site.summary())

    t5 = site.tau5()
    t4 = site.tau4()
    if t5 is not None:
        print(f"  tau_5 = {t5:.3f}")
    if t4 is not None:
        print(f"  tau_4 = {t4:.3f}")

    return 0


def _cmd_list_metals(args) -> int:
    filepath = _resolve_cif(args)

    if args.mmcif:
        from terrahedral.parsers.mmcif_parser import parse, find_metals
        data = parse(filepath)
        metals = find_metals(data)
        if not metals:
            print("No metal HETATM atoms found.")
            return 0
        print(f"Found {len(metals)} metal atom(s):")
        for m in metals:
            print(f"  {m['comp_id']} asym={m['asym_id']} seq={m['seq_id']} "
                  f"({m['element']}) at ({m['x']:.3f}, {m['y']:.3f}, {m['z']:.3f})")
    else:
        from terrahedral.parsers.cif_parser import parse
        from terrahedral.core import METAL_ELEMENTS
        data = parse(filepath)
        metals = [
            a for a in data.atoms
            if a["element"] in METAL_ELEMENTS
        ]
        if not metals:
            print("No metal atoms found.")
            return 0
        print(f"Found {len(metals)} metal atom(s):")
        for m in metals:
            print(f"  {m['label']} ({m['element']}) at "
                  f"({m['frac_x']:.5f}, {m['frac_y']:.5f}, {m['frac_z']:.5f})")

    return 0


def _cmd_fetch(args) -> int:
    import shutil

    if not args.pdb and not args.cod:
        print("Error: provide --pdb ID or --cod ID", file=sys.stderr)
        return 1
    if args.pdb and args.cod:
        print("Error: provide only one of --pdb or --cod", file=sys.stderr)
        return 1

    if args.pdb:
        from terrahedral.fetch import fetch_pdb
        path = fetch_pdb(args.pdb, cache_dir=args.cache_dir, force=args.no_cache)
    else:
        from terrahedral.fetch import fetch_cod
        path = fetch_cod(args.cod, cache_dir=args.cache_dir, force=args.no_cache)

    if args.output:
        shutil.copy2(path, args.output)
        print(f"Saved to {args.output}")
    else:
        print(f"Downloaded to {path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
