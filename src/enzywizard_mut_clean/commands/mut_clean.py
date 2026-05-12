from __future__ import annotations
from argparse import ArgumentParser, Namespace
from ..services.mut_clean_service import run_mutclean_service

def add_mutclean_parser(parser: ArgumentParser) -> None:
    parser.add_argument("-w","--wt_input_path", required=True, help="Path to wild-type protein structure file (CIF/PDB).")
    parser.add_argument("-m","--mut_input_path", required=True, help="Path to mutant protein structure file (CIF/PDB).")
    parser.add_argument("-s","--mutation", required=True, help="Amino acid substitution(s) describing mutations, e.g., A123V or A123V,G456D.")
    parser.add_argument("-o","--output_dir", required=True, help="Path to a directory for outputting cleaned CIF, PDB, and FASTA files and a JSON report.")
    parser.add_argument("--no_add_H",action="store_false",dest="add_H",help="Disable adding hydrogens using OpenMM (default: enabled).")
    parser.set_defaults(add_H=True)
    parser.add_argument("--pH",type=float,default=7.0,help="pH value for hydrogen addition (default: 7.0).")
    parser.set_defaults(func=run_mutclean)

def run_mutclean(args: Namespace) -> None:
    run_mutclean_service(wt_input_path=args.wt_input_path, mut_input_path=args.mut_input_path,mutation=args.mutation,output_dir=args.output_dir,add_H=args.add_H, pH=args.pH)





