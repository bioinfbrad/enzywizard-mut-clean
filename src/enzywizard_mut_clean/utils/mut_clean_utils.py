import re
from ..utils.logging_utils import Logger
from typing import List, Tuple, Any, Dict

def check_amino_acid_substitution(mutation: str,wt_length: int,mut_length: int,logger: Logger) -> bool:

    if not isinstance(mutation, str) or not mutation.strip():
        logger.print("[ERROR] Empty amino acid substitution.")
        return False

    muts = mutation.split(",")
    seen_positions = set()
    for single_mut in muts:
        single_mut = single_mut.strip()

        m = re.fullmatch(
            r"([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])",
            single_mut,
            flags=re.I,
        )
        if not m:
            logger.print(f"[ERROR] Invalid amino acid substitution format: {single_mut}.")
            return False

        wt_aa = m.group(1).upper()
        pos_str = m.group(2)
        mut_aa = m.group(3).upper()

        if wt_aa == mut_aa:
            logger.print(f"[ERROR] No amino acid mutated: {single_mut}.")
            return False

        pos = int(pos_str)

        if pos in seen_positions:
            logger.print(f"[ERROR] Duplicate mutation position detected: {pos}.")
            return False
        seen_positions.add(pos)

        if pos < 1 or pos > wt_length:
            logger.print(f"[ERROR] Mutation position out of wild-type range: {single_mut} (position={pos}, wt_length={wt_length})")
            return False

        if pos < 1 or pos > mut_length:
            logger.print(f"[ERROR] Mutation position out of mutant range: {single_mut} (position={pos}, mut_length={mut_length})")
            return False
    return True

def get_muts_from_aas(mutation: str)->List[Tuple[str, int, str]]:
    muts_list: List[Tuple[str, int, str]]=[]
    muts = mutation.split(",")

    for single_mut in muts:
        single_mut = single_mut.strip()

        m = re.fullmatch(
            r"([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])",
            single_mut,
            flags=re.I,
        )

        wt_aa = m.group(1).upper()
        pos_str = m.group(2)
        mut_aa = m.group(3).upper()
        muts_list.append((wt_aa,int(pos_str),mut_aa))
    return muts_list



def postprocess_mutclean_report_to_schema(
    raw_report: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Map the raw EnzyWizard-Mut-Clean report to the new JSON Schema field names.

    This function only renames/restructures fields and does not validate or convert values.
    """

    def map_residue(raw_residue: Dict[str, Any]) -> Dict[str, Any]:
        return {
            "residue_index": raw_residue.get("aa_id"),
            "residue_name": raw_residue.get("aa_name"),
            "hydrogen_atom_count": raw_residue.get("hydrogen_atom_count"),
        }

    def map_residue_mapping(raw_mapping: list[Dict[str, Any]]) -> list[Dict[str, Any]]:
        schema_mapping: list[Dict[str, Any]] = []

        for mapping_item in raw_mapping:
            schema_mapping.append(
                {
                    "old_residue": map_residue(mapping_item.get("old_residue", {})),
                    "new_residue": map_residue(mapping_item.get("new_residue", {})),
                }
            )

        return schema_mapping

    def map_clean_statistics(raw_statistics: Dict[str, Any]) -> Dict[str, Any]:
        return {
            "removed_heterogen_count": raw_statistics.get("removed_heterogen"),
            "standardized_residue_name_count": raw_statistics.get("changed_resname"),
            "repaired_residue_count": raw_statistics.get("fixed_residues"),
            "added_heavy_atom_count": raw_statistics.get("added_heavy_atoms"),
            "added_hydrogen_atom_count": raw_statistics.get("added_hydrogen_atoms"),
            "retained_residue_count": raw_statistics.get("kept_residues"),
        }

    schema_report: Dict[str, Any] = {
        "report_type": raw_report.get("output_type"),
        "amino_acid_substitution": raw_report.get("amino_acid_substitution"),
        "cleaned_amino_acid_substitution": raw_report.get("cleaned_amino_acid_substitution"),
        "wild_type_residue_mapping_old_to_new": map_residue_mapping(
            raw_report.get("wt_amino_acid_mapping_old_to_new", [])
        ),
        "wild_type_clean_statistics": map_clean_statistics(
            raw_report.get("wt_clean_statistics", {})
        ),
        "mutant_residue_mapping_old_to_new": map_residue_mapping(
            raw_report.get("mut_amino_acid_mapping_old_to_new", [])
        ),
        "mutant_clean_statistics": map_clean_statistics(
            raw_report.get("mut_clean_statistics", {})
        ),
    }

    return schema_report