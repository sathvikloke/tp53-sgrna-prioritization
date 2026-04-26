"""Unified command-line interface for the TP53 sgRNA pipeline.

Subcommands:
    fetch-cds         Download NM_000546.6 from NCBI.
    fetch-cohorts     Pull cBioPortal mutation calls for the configured cohorts.
    enumerate         Find all candidate sgRNAs around the hotspot codons.
    score-on-target   Score candidates with the on-target model.
    score-off-target  Score candidates against GRCh38 with Cas-OFFinder + CFD.
    rank              Compute per-cohort composite scores and rankings.
    tables            Build summary JSON + figures.
    all               Run every step in order.

Run `python -m src.cli <subcommand> --help` for per-step options.
"""

from __future__ import annotations

import argparse
import logging
import sys
from typing import Callable, Dict


def _set_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level, format="%(asctime)s %(levelname)s %(name)s: %(message)s"
    )


def _cmd_fetch_cds(args) -> int:
    from . import fetch_cds
    fetch_cds.run()
    return 0


def _cmd_fetch_cohorts(args) -> int:
    from . import fetch_cohort_mutations
    fetch_cohort_mutations.run()
    return 0


def _cmd_enumerate(args) -> int:
    from . import enumerate_sgrnas
    enumerate_sgrnas.run()
    return 0


def _cmd_score_on(args) -> int:
    from . import score_on_target
    score_on_target.run()
    return 0


def _cmd_score_off(args) -> int:
    from . import score_off_target
    score_off_target.run(device=args.device)
    return 0


def _cmd_rank(args) -> int:
    from . import compose_rankings
    compose_rankings.run()
    return 0


def _cmd_tables(args) -> int:
    from . import make_tables
    make_tables.run()
    return 0


def _cmd_all(args) -> int:
    """Run every step; abort on the first failure."""
    steps: Dict[str, Callable] = {
        "fetch-cds": _cmd_fetch_cds,
        "fetch-cohorts": _cmd_fetch_cohorts,
        "enumerate": _cmd_enumerate,
        "score-on-target": _cmd_score_on,
        "score-off-target": _cmd_score_off,
        "rank": _cmd_rank,
        "tables": _cmd_tables,
    }
    for name, fn in steps.items():
        logging.info("=== running %s ===", name)
        rc = fn(args)
        if rc != 0:
            logging.error("step %s failed (rc=%d); aborting.", name, rc)
            return rc
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="tp53-sgrna",
                                description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-v", "--verbose", action="store_true")
    sub = p.add_subparsers(dest="cmd", required=True)

    sub.add_parser("fetch-cds").set_defaults(func=_cmd_fetch_cds)
    sub.add_parser("fetch-cohorts").set_defaults(func=_cmd_fetch_cohorts)
    sub.add_parser("enumerate").set_defaults(func=_cmd_enumerate)
    sub.add_parser("score-on-target").set_defaults(func=_cmd_score_on)

    s = sub.add_parser("score-off-target")
    s.add_argument("--device", default="C", choices=["C", "G", "A"],
                   help="Cas-OFFinder device flag (CPU/GPU/Accelerator).")
    s.set_defaults(func=_cmd_score_off)

    sub.add_parser("rank").set_defaults(func=_cmd_rank)
    sub.add_parser("tables").set_defaults(func=_cmd_tables)

    a = sub.add_parser("all")
    a.add_argument("--device", default="C", choices=["C", "G", "A"])
    a.set_defaults(func=_cmd_all)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    _set_logging(args.verbose)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
