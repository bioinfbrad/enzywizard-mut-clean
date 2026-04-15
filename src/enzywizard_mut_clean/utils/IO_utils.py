from __future__ import annotations

from Bio.PDB import MMCIFParser, PDBParser, MMCIFIO, PDBIO
from Bio.PDB.Structure import Structure
from pathlib import Path

from ..utils.logging_utils import Logger
import json
import tempfile
from ..utils.common_utils import convert_to_json_serializable, InlineJSONEncoder, wrap_leaf_lists_as_rawjson, get_clean_filename, get_optimized_filename
from openmm.app import PDBFile,PDBxFile, Modeller
from ..utils.structure_utils import get_single_chain,get_residues_by_chain,get_sequence




def file_exists(path: str | Path) -> bool:
    p = Path(path)
    return p.exists() and p.is_file()

def get_stem(input_path: str | Path) -> str:
    return Path(input_path).stem

MAXFILENAME=150

def check_filename_length(name: str, logger: Logger) -> bool:
    if len(name) > MAXFILENAME:
        logger.print(f"[ERROR] Filename too long (>{MAXFILENAME}): {name}")
        return False
    return True

def load_protein_structure(path: str | Path, protein_name:str, logger: Logger) -> Structure | None:
    p = Path(path)

    try:
        if p.suffix.lower() in {".cif", ".mmcif"}:
            parser = MMCIFParser(QUIET=True)
        elif p.suffix.lower() == ".pdb":
            parser = PDBParser(QUIET=True)
        else:
            logger.print(f"[ERROR] Unsupported format: {p}")
            return None

        return parser.get_structure(protein_name, str(p))

    except Exception as e:
        logger.print(f"[ERROR] Exception in loading structure for {str(p)}: {e}")
        return None



def structure_to_pdbfile(struct: Structure, logger: Logger, protein_name: str = "structure") -> PDBFile | None:
    try:
        with tempfile.NamedTemporaryFile(suffix=".pdb") as tmp:
            pdb_path = Path(tmp.name)

            io = PDBIO()
            io.set_structure(struct)
            io.save(str(pdb_path))

            pdb = PDBFile(str(pdb_path))
            return pdb

    except Exception as e:
        logger.print(f"[ERROR] Failed to convert Structure to PDBFile: {e}")
        return None



def modeller_to_structure(modeller: Modeller, logger: Logger, protein_name: str = "structure") -> Structure | None:
    try:
        with tempfile.NamedTemporaryFile(suffix=".pdb") as tmp:
            pdb_path = Path(tmp.name)

            with open(pdb_path, "w") as f:
                PDBFile.writeFile(modeller.topology, modeller.positions, f)

            parser = PDBParser(QUIET=True)
            struct = parser.get_structure(protein_name, str(pdb_path))

            return struct

    except Exception as e:
        logger.print(f"[ERROR] Failed to convert Modeller to Structure: {e}")
        return None



def write_cif(struct: Structure, output_path: str | Path) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    io = MMCIFIO()
    io.set_structure(struct)
    io.save(str(output_path))

def write_pdb(struct: Structure, output_path: str | Path) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    io = PDBIO()
    io.set_structure(struct)
    io.save(str(output_path))



def write_json_from_dict_inline_leaf_lists(dict_data: dict, output_path: str | Path) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    dict_data = convert_to_json_serializable(dict_data)
    dict_data = wrap_leaf_lists_as_rawjson(dict_data)

    with output_path.open("w", encoding="utf-8") as f:
        json.dump(
            dict_data,
            f,
            cls=InlineJSONEncoder,
            indent=2,
            ensure_ascii=False
        )

def write_fasta(struct: Structure, output_path: str | Path, logger: Logger) -> bool:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    chain = get_single_chain(struct, logger)
    if chain is None:
        return False

    residues = get_residues_by_chain(chain, logger)
    if residues is None:
        return False

    seq = get_sequence(residues, logger)
    if seq is None:
        return False

    header = str(struct.id).strip() if getattr(struct, "id", None) else None
    if header is None:
        logger.print("[ERROR] Structure header is None")
        return False

    try:
        with output_path.open("w", encoding="utf-8") as f:
            f.write(f">{header}\n")
            f.write(f"{seq}\n")
        return True
    except Exception as e:
        logger.print(f"[ERROR] Failed to write FASTA to {output_path}: {e}")
        return False

