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
    logger: Logger,
) -> Dict[str, Any] | None:
    """
    Convert the raw EnzyWizard-Mut-Clean report to the new JSON Schema format.

    """

    if not isinstance(raw_report, dict):
        logger.print("[ERROR] Invalid raw Mut-Clean report: expected a dictionary.")
        return None

    def convert_residue_info(raw_residue: Dict[str, Any]) -> Dict[str, Any] | None:
        if not isinstance(raw_residue, dict):
            logger.print("[ERROR] Invalid residue information in raw Mut-Clean report.")
            return None

        residue_index = raw_residue.get("aa_id")
        residue_name = raw_residue.get("aa_name")
        hydrogen_atom_count = raw_residue.get("hydrogen_atom_count")

        if not isinstance(residue_index, int):
            logger.print("[ERROR] Invalid or missing residue index in raw Mut-Clean report.")
            return None

        if not isinstance(residue_name, str):
            logger.print("[ERROR] Invalid or missing residue name in raw Mut-Clean report.")
            return None

        if not isinstance(hydrogen_atom_count, int):
            logger.print("[ERROR] Invalid or missing hydrogen atom count in raw Mut-Clean report.")
            return None

        return {
            "residue_index": residue_index,
            "residue_name": residue_name,
            "hydrogen_atom_count": hydrogen_atom_count,
        }

    def convert_residue_mapping(
        raw_mapping: List[Dict[str, Any]],
        mapping_name: str,
    ) -> List[Dict[str, Any]] | None:
        if not isinstance(raw_mapping, list):
            logger.print(f"[ERROR] Invalid {mapping_name}: expected a list.")
            return None

        schema_mapping: List[Dict[str, Any]] = []

        for mapping_item in raw_mapping:
            if not isinstance(mapping_item, dict):
                logger.print(f"[ERROR] Invalid item in {mapping_name}: expected a dictionary.")
                return None

            raw_old_residue = mapping_item.get("old_residue")
            raw_new_residue = mapping_item.get("new_residue")

            old_residue = convert_residue_info(raw_old_residue)
            if old_residue is None:
                logger.print(f"[ERROR] Failed to convert old_residue in {mapping_name}.")
                return None

            new_residue = convert_residue_info(raw_new_residue)
            if new_residue is None:
                logger.print(f"[ERROR] Failed to convert new_residue in {mapping_name}.")
                return None

            schema_mapping.append({
                "old_residue": old_residue,
                "new_residue": new_residue,
            })

        return schema_mapping

    def convert_clean_statistics(
        raw_statistics: Dict[str, Any],
        statistics_name: str,
    ) -> Dict[str, Any] | None:
        if not isinstance(raw_statistics, dict):
            logger.print(f"[ERROR] Invalid {statistics_name}: expected a dictionary.")
            return None

        raw_to_schema_key = {
            "removed_heterogen": "removed_heterogen_count",
            "changed_resname": "standardized_residue_name_count",
            "fixed_residues": "repaired_residue_count",
            "added_heavy_atoms": "added_heavy_atom_count",
            "added_hydrogen_atoms": "added_hydrogen_atom_count",
            "kept_residues": "retained_residue_count",
        }

        schema_statistics: Dict[str, Any] = {}

        for raw_key, schema_key in raw_to_schema_key.items():
            value = raw_statistics.get(raw_key)

            if not isinstance(value, int):
                logger.print(
                    f"[ERROR] Invalid or missing value in {statistics_name}: "
                    f"{raw_key}={value}"
                )
                return None

            schema_statistics[schema_key] = value

        return schema_statistics

    report_type = raw_report.get("output_type")
    amino_acid_substitution = raw_report.get("amino_acid_substitution")
    cleaned_amino_acid_substitution = raw_report.get("cleaned_amino_acid_substitution")

    if report_type != "enzywizard_mut_clean":
        logger.print(f"[ERROR] Invalid output_type in raw Mut-Clean report: {report_type}")
        return None

    if not isinstance(amino_acid_substitution, str) or not amino_acid_substitution.strip():
        logger.print("[ERROR] Invalid amino_acid_substitution in raw Mut-Clean report.")
        return None

    if not isinstance(cleaned_amino_acid_substitution, str) or not cleaned_amino_acid_substitution.strip():
        logger.print("[ERROR] Invalid cleaned_amino_acid_substitution in raw Mut-Clean report.")
        return None

    wild_type_residue_mapping_old_to_new = convert_residue_mapping(
        raw_report.get("wt_amino_acid_mapping_old_to_new"),
        "wt_amino_acid_mapping_old_to_new",
    )
    if wild_type_residue_mapping_old_to_new is None:
        return None

    mutant_residue_mapping_old_to_new = convert_residue_mapping(
        raw_report.get("mut_amino_acid_mapping_old_to_new"),
        "mut_amino_acid_mapping_old_to_new",
    )
    if mutant_residue_mapping_old_to_new is None:
        return None

    wild_type_clean_statistics = convert_clean_statistics(
        raw_report.get("wt_clean_statistics"),
        "wt_clean_statistics",
    )
    if wild_type_clean_statistics is None:
        return None

    mutant_clean_statistics = convert_clean_statistics(
        raw_report.get("mut_clean_statistics"),
        "mut_clean_statistics",
    )
    if mutant_clean_statistics is None:
        return None

    return {
        "report_type": report_type,
        "amino_acid_substitution": amino_acid_substitution,
        "cleaned_amino_acid_substitution": cleaned_amino_acid_substitution,
        "wild_type_residue_mapping_old_to_new": wild_type_residue_mapping_old_to_new,
        "wild_type_clean_statistics": wild_type_clean_statistics,
        "mutant_residue_mapping_old_to_new": mutant_residue_mapping_old_to_new,
        "mutant_clean_statistics": mutant_clean_statistics,
    }