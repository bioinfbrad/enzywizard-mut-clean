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


# ==============================
# Command: enzywizard-mut-clean
# ==============================

# brief introduction:
'''
EnzyWizard-Mut-Clean is a command-line tool for cleaning a wild-type protein 
structure and its corresponding mutant protein structure, generating multi-format 
cleaned protein files (CIF, PDB, and FASTA), remapping amino acid substitution(s), 
and providing a detailed traceable cleaning report that records residue mapping 
between the original and cleaned structures. Using PDBFixer APIs, it removes 
heterogens, repairs missing heavy atoms, optionally adds hydrogens, and converts 
both structures into continuous protein chains with standardized residue numbering.
The tool also remaps the input amino acid substitution(s) to the cleaned residue 
numbering system and outputs a JSON report summarizing residue mapping and cleaning 
statistics for both wild-type and mutant structures.

'''

# example usage:
'''
Example command:

enzywizard-mut-clean -w examples/input/1ZG4_WT.cif -m examples/input/1ZG6_S70G.cif -s S70G -o examples/output/

'''

# input parameters:
'''
-w, --wt_input_path
Required.
Path to the input wild-type protein structure file in CIF or PDB format.

-m, --mut_input_path
Required.
Path to the input mutant protein structure file in CIF or PDB format.

-s, --mutation
Required.
Amino acid substitution(s) describing the mutation relationship between the
wild-type and mutant structures.
Examples:
S70G
S70G,A123V
Multiple mutations should be separated by ','.

-o, --output_dir
Required.
Path to the output directory for saving cleaned structure files and the JSON report.

--no_add_H
Optional flag.
Disable hydrogen addition using OpenMM.
By default, hydrogens are added.

--pH
Optional.
pH value used for hydrogen addition.
Default: 7.0
Valid range: 0.0 to 14.0
'''

# output content:
'''
The program outputs the following files into the output directory:

1. Cleaned structure files for the wild-type structure
   - cleaned_{wt_name}.cif
   - cleaned_{wt_name}.pdb
   - cleaned_{wt_name}.fasta

2. Cleaned structure files for the mutant structure
   - cleaned_{mut_name}.cif
   - cleaned_{mut_name}.pdb
   - cleaned_{mut_name}.fasta

3. A JSON report
   - mut_clean_report_{wt_name}_to_{mut_name}.json

   The JSON report contains:

   - "output_type"
     A string identifying the report type:
     "enzywizard_mut_clean"

   - "amino_acid_substitution"
     The original amino acid substitution(s) provided by the user.

   - "cleaned_amino_acid_substitution"
     The amino acid substitution(s) after remapping from the original residue
     numbering system to the cleaned residue numbering system.

   - "wt_amino_acid_mapping_old_to_new"
     A list describing how residues in the original wild-type structure
     correspond to residues in the cleaned wild-type structure.

     Each entry contains:
     - "old_residue"
       Information for the original residue before cleaning:
       - "aa_id": original residue index
       - "aa_name": original residue one-letter amino acid code
       - "hydrogen_atom_count": number of hydrogen atoms found in the original residue

     - "new_residue"
       Information for the cleaned residue after cleaning:
       - "aa_id": new residue index after continuous renumbering
       - "aa_name": cleaned residue one-letter amino acid code
       - "hydrogen_atom_count": number of hydrogen atoms in the cleaned residue

   - "wt_clean_statistics"
     A dictionary summarizing the overall cleaning process for the wild-type structure.

     It includes:
     - "removed_heterogen"
       Number of non-protein residues removed.

     - "changed_resname"
       Number of residues identified as non-standard and replaced by PDBFixer.

     - "fixed_residues"
       Number of residues where missing atoms were repaired.

     - "added_heavy_atoms"
       Total number of heavy atoms added during missing atom reconstruction.

     - "added_hydrogen_atoms"
       Number of hydrogen atoms added (if hydrogen addition is enabled).

     - "kept_residues"
       Number of residues kept in the final cleaned wild-type structure.

   - "mut_amino_acid_mapping_old_to_new"
     A list describing how residues in the original mutant structure
     correspond to residues in the cleaned mutant structure.

     Each entry contains:
     - "old_residue"
       Information for the original residue before cleaning:
       - "aa_id": original residue index
       - "aa_name": original residue one-letter amino acid code
       - "hydrogen_atom_count": number of hydrogen atoms found in the original residue

     - "new_residue"
       Information for the cleaned residue after cleaning:
       - "aa_id": new residue index after continuous renumbering
       - "aa_name": cleaned residue one-letter amino acid code
       - "hydrogen_atom_count": number of hydrogen atoms in the cleaned residue

   - "mut_clean_statistics"
     A dictionary summarizing the overall cleaning process for the mutant structure.

     It includes:
     - "removed_heterogen"
       Number of non-protein residues removed.

     - "changed_resname"
       Number of residues identified as non-standard and replaced by PDBFixer.

     - "fixed_residues"
       Number of residues where missing atoms were repaired.

     - "added_heavy_atoms"
       Total number of heavy atoms added during missing atom reconstruction.

     - "added_hydrogen_atoms"
       Number of hydrogen atoms added (if hydrogen addition is enabled).

     - "kept_residues"
       Number of residues kept in the final cleaned mutant structure.

   This report helps users track:
   - how residue numbering changed in both wild-type and mutant structures,
   - whether residue names were standardized,
   - how hydrogen content changed after optional hydrogen addition,
   - and how the input mutation positions were translated into the cleaned
     residue numbering system.
'''

# Process:
'''
This command processes the input wild-type and mutant protein structures as follows:

1. Load the input structures
   - Read the wild-type CIF or PDB file using Biopython and PDBFixer.
   - Read the mutant CIF or PDB file using Biopython and PDBFixer.
   - Resolve the protein names from the input filenames.

2. Validate basic input conditions
   - Check that both input files exist.
   - Check that the wild-type and mutant filenames are valid and not identical.
   - Check that the pH value is within the valid range.
   - Extract a single chain from each original structure.
   - Check that the wild-type and mutant sequence lengths are valid and equal.
   - Check that the input amino acid substitution format is valid.
   - Check that mutation positions fall within the residue index ranges of both structures.

3. Clean both structures independently (PDBFixer-based processing)
   - Keep only the first chain.
   - Identify non-standard residues.
   - Remove all heterogens.
   - Replace non-standard residues with standard residues.
   - Disable missing residue reconstruction (do not add missing residues).
   - Detect missing atoms.
   - Add missing heavy atoms.
   - Optionally add hydrogens using OpenMM ForceField with specified pH.
   - Check for invalid coordinates (e.g., NaN values).

4. Renumber both structures
   - Rebuild each topology to ensure:
     - single chain (chain ID A),
     - continuous residue numbering starting from 1.

5. Remap amino acid substitution
   - Build residue mapping relationships from original residue numbering to cleaned
     residue numbering for both wild-type and mutant structures.
   - Ensure that each mutation position can be found in both mapping results.
   - Ensure that wild-type and mutant mappings produce the same cleaned residue
     position for each input mutation.
   - Convert the input amino acid substitution(s) into cleaned amino acid
     substitution(s).

6. Validate the cleaned structures
   - Load the cleaned wild-type and mutant structures from the saved CIF files.
   - Ensure each cleaned structure contains exactly one model and one chain.
   - Ensure each cleaned chain ID is A.
   - Ensure no hetero residues remain.
   - Ensure residue names are standardized.
   - Ensure residue numbering is continuous.
   - Ensure all required backbone and heavy atoms are present.
   - Check that the cleaned wild-type and mutant sequence lengths are equal.

7. Save outputs
   - Save the cleaned wild-type structure in CIF format.
   - Save the cleaned wild-type structure in PDB format.
   - Extract and save the cleaned wild-type amino acid sequence in FASTA format.
   - Save the cleaned mutant structure in CIF format.
   - Save the cleaned mutant structure in PDB format.
   - Extract and save the cleaned mutant amino acid sequence in FASTA format.
   - Generate and save the JSON report summarizing mutation remapping,
     residue mappings, and cleaning statistics for both structures.
'''

# dependencies:
'''
- Biopython
- OpenMM
- PDBFixer
- NumPy
'''

# references:
'''
- Biopython:
  https://biopython.org/

- OpenMM:
  https://openmm.org/

- PDBFixer:
  https://github.com/openmm/pdbfixer

- wwPDB Chemical Component Dictionary:
  https://www.wwpdb.org/data/ccd
'''


