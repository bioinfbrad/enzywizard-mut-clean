from __future__ import annotations
from pathlib import Path
from typing import List, Dict, Any

from Bio.Data.IUPACData import protein_letters_3to1





def normalize_aa_name_to_one_letter(aa_name: Any) -> str:

    aa_name_clean = aa_name.strip()

    if len(aa_name_clean) == 1:
        return aa_name_clean.upper()

    if len(aa_name_clean) == 3:
        aa_name_3 = aa_name_clean.upper().capitalize()
        aa_name_1 = protein_letters_3to1.get(aa_name_3)
        if isinstance(aa_name_1, str) and aa_name_1 != "":
            return aa_name_1.upper()
    return "X"