from __future__ import annotations

import argparse

from .commands.mut_clean import add_mutclean_parser


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="enzywizard-mut-clean",
        description="EnzyWizard-Mut-Clean: Clean both a wild-type protein structure and its corresponding mutant protein structure, remap cleaned amino acid substitution(s), generate multi-format cleaned protein files (CIF, PDB, and FASTA), and provide a detailed traceable cleaning report."
    )
    add_mutclean_parser(parser)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)