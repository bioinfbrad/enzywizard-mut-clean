from __future__ import annotations
from pathlib import Path
from ..utils.logging_utils import Logger
from ..utils.IO_utils import file_exists, get_stem, check_filename_length, load_protein_structure, write_json_from_dict_inline_leaf_lists, write_fasta, load_pdbfixer, write_cif_from_pdbfixer, write_pdb_from_pdbfixer
from ..utils.mut_clean_utils import check_amino_acid_substitution
from ..utils.structure_utils import get_single_chain,get_chain_length
from ..algorithms.clean_algorithms import clean_pdbfixer_to_single_chain_A, check_cleaned_structure
from ..algorithms.mut_clean_algorithms import get_cleaned_amino_acid_substitution, generate_mutclean_report
from ..utils.common_utils import get_optimized_filename


def run_mutclean_service(wt_input_path: str | Path, mut_input_path: str | Path, mutation:str, output_dir: str | Path, add_H: bool = True, pH: float = 7.0, force_field_file: str = "charmm36.xml")->bool:
    # ---- logger ----
    logger = Logger(output_dir)
    logger.print(f"[INFO] Mutclean processing started: {wt_input_path} {mut_input_path} {mutation}")

    # ---- check input ----
    wt_input_path=Path(wt_input_path)
    mut_input_path=Path(mut_input_path)
    output_dir = Path(output_dir)

    if not (0.0 <= pH <= 14.0):
        logger.print(f"[ERROR] Invalid pH value: {pH}. Must be between 0 and 14.")
        return False

    if not file_exists(wt_input_path):
        logger.print(f"[ERROR] Input not found: {wt_input_path}")
        return False

    if not file_exists(mut_input_path):
        logger.print(f"[ERROR] Input not found: {mut_input_path}")
        return False

    output_dir.mkdir(parents=True, exist_ok=True)

    # ---- get name ----
    wt_name = get_stem(wt_input_path)
    if not check_filename_length(wt_name,logger):
        return False
    mut_name = get_stem(mut_input_path)
    if not check_filename_length(mut_name,logger):
        return False
    if wt_name==mut_name:
        logger.print(f"[ERROR] Wild-type and mutant names are the same: {wt_name} {mut_name}")
        return False
    logger.print(f"[INFO] Protein names resolved: {wt_name} {mut_name}")

    # ---- load structure ----
    wt_structure = load_protein_structure(wt_input_path,wt_name,logger)
    mut_structure = load_protein_structure(mut_input_path,mut_name,logger)
    if wt_structure is None:
        logger.print(f"[ERROR] Failed to load wild-type structure: {wt_input_path}")
        return False
    if mut_structure is None:
        logger.print(f"[ERROR] Failed to load mutant structure: {mut_input_path}")
        return False

    wt_fixer = load_pdbfixer(wt_input_path, logger)
    if wt_fixer is None:
        logger.print(f"[ERROR] Failed to load wild-type PDBFixer: {wt_input_path}")
        return False

    mut_fixer = load_pdbfixer(mut_input_path, logger)
    if mut_fixer is None:
        logger.print(f"[ERROR] Failed to load mutant PDBFixer: {mut_input_path}")
        return False

    logger.print("[INFO] Structures loaded")

    # ---- check amino acid substitution ----
    wt_chain = get_single_chain(wt_structure, logger)
    mut_chain = get_single_chain(mut_structure, logger)
    if wt_chain is None or mut_chain is None:
        return False

    wt_seq_length=get_chain_length(wt_chain,logger)
    mut_seq_length=get_chain_length(mut_chain,logger)
    if wt_seq_length is None or mut_seq_length is None:
        return False
    if not check_amino_acid_substitution(mutation,wt_seq_length,mut_seq_length,logger):
        return False

    # ---- check length ----
    if wt_seq_length!=mut_seq_length:
        logger.print(f"[ERROR] Wild-type and mutant sequence lengths are not equal: {wt_seq_length} {mut_seq_length}")
        return False

    # ---- run algorithm ----
    logger.print("[INFO] Cleaning wild-type and mutant structures started")

    wt_clean_result = clean_pdbfixer_to_single_chain_A(fixer=wt_fixer,add_H=add_H,logger=logger,pH=pH,force_field_file=force_field_file)
    if wt_clean_result is None:
        return False

    wt_cleaned_fixer, wt_mapping_old_to_new, wt_stats = wt_clean_result

    mut_clean_result = clean_pdbfixer_to_single_chain_A(fixer=mut_fixer,add_H=add_H,logger=logger,pH=pH,force_field_file=force_field_file)
    if mut_clean_result is None:
        return False

    mut_cleaned_fixer, mut_mapping_old_to_new, mut_stats = mut_clean_result

    cleaned_mutation = get_cleaned_amino_acid_substitution(wt_mapping_old_to_new,mut_mapping_old_to_new,mutation,logger)
    if cleaned_mutation is None:
        return False

    # save files

    cleaned_wt_cif_path = output_dir / get_optimized_filename(f"cleaned_{wt_name}.cif")
    cleaned_wt_pdb_path = output_dir / get_optimized_filename(f"cleaned_{wt_name}.pdb")

    cleaned_mut_cif_path = output_dir / get_optimized_filename(f"cleaned_{mut_name}.cif")
    cleaned_mut_pdb_path = output_dir / get_optimized_filename(f"cleaned_{mut_name}.pdb")

    write_cif_from_pdbfixer(wt_cleaned_fixer, cleaned_wt_cif_path)
    logger.print(f"[INFO] Cleaned wild-type CIF saved: {cleaned_wt_cif_path}")

    write_pdb_from_pdbfixer(wt_cleaned_fixer, cleaned_wt_pdb_path)
    logger.print(f"[INFO] Cleaned wild-type PDB saved: {cleaned_wt_pdb_path}")

    write_cif_from_pdbfixer(mut_cleaned_fixer, cleaned_mut_cif_path)
    logger.print(f"[INFO] Cleaned mutant CIF saved: {cleaned_mut_cif_path}")

    write_pdb_from_pdbfixer(mut_cleaned_fixer, cleaned_mut_pdb_path)
    logger.print(f"[INFO] Cleaned mutant PDB saved: {cleaned_mut_pdb_path}")

    # load cleaned structure

    wt_cleaned_structure = load_protein_structure(cleaned_wt_cif_path, wt_name, logger)
    if wt_cleaned_structure is None:
        logger.print(f"[ERROR] Failed to load cleaned wild-type structure: {cleaned_wt_cif_path}")
        return False

    mut_cleaned_structure = load_protein_structure(cleaned_mut_cif_path, mut_name, logger)
    if mut_cleaned_structure is None:
        logger.print(f"[ERROR] Failed to load cleaned mutant structure: {cleaned_mut_cif_path}")
        return False

    logger.print("[INFO] Cleaned structures loaded from saved CIF files")

    # ---- post check ----
    wt_cleaned_chain = get_single_chain(wt_cleaned_structure, logger)
    mut_cleaned_chain = get_single_chain(mut_cleaned_structure, logger)
    if wt_cleaned_chain is None or mut_cleaned_chain is None:
        return False

    wt_cleaned_length=get_chain_length(wt_cleaned_chain,logger)
    mut_cleaned_length=get_chain_length(mut_cleaned_chain,logger)
    if wt_cleaned_length is None or mut_cleaned_length is None:
        return False
    if wt_cleaned_length!=mut_cleaned_length:
        logger.print(f"[ERROR] Cleaned wild-type and mutant sequence lengths are not equal: {wt_cleaned_length} {mut_cleaned_length}")
        return False

    if not check_cleaned_structure(wt_cleaned_structure, logger):
        return False
    if not check_cleaned_structure(mut_cleaned_structure, logger):
        return False

    # ---- save results ----

    cleaned_wt_fasta_path = output_dir / get_optimized_filename(f"cleaned_{wt_name}.fasta")

    cleaned_mut_fasta_path = output_dir / get_optimized_filename(f"cleaned_{mut_name}.fasta")

    json_report_path = output_dir / get_optimized_filename(f"mut_clean_report_{wt_name}_to_{mut_name}.json")

    if write_fasta(wt_cleaned_structure,cleaned_wt_fasta_path, logger):
        logger.print(f"[INFO] Cleaned wild-type FASTA saved: {cleaned_wt_fasta_path}")
    else:
        return False

    if write_fasta(mut_cleaned_structure,cleaned_mut_fasta_path, logger):
        logger.print(f"[INFO] Cleaned mutant FASTA saved: {cleaned_mut_fasta_path}")
    else:
        return False

    report=generate_mutclean_report(mutation,cleaned_mutation,wt_structure,wt_cleaned_structure,mut_structure,mut_cleaned_structure,wt_mapping_old_to_new,wt_stats,mut_mapping_old_to_new,mut_stats,logger)
    if report is None:
        return False
    write_json_from_dict_inline_leaf_lists(report,json_report_path)
    logger.print(f"[INFO] Report JSON saved: {json_report_path}")

    logger.print("[INFO] Mutclean processing finished")

    return True
