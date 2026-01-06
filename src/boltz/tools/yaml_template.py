from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import yaml
from rdkit import Chem


def _chain_id_from_index(index: int) -> str:
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    if index < len(letters):
        return letters[index]
    return f"{letters[index % len(letters)]}{index // len(letters)}"


def _extract_chain_id(header: str) -> Optional[str]:
    if not header:
        return None
    token = header.split("|")[0].strip()
    if not token:
        parts = header.split()
        token = parts[0].strip() if parts else ""
    return token or None


def _parse_fasta(path: Path) -> list[tuple[str, str]]:
    headers: list[str] = []
    sequences: list[str] = []
    current: list[str] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    sequences.append("".join(current))
                    current = []
                headers.append(line[1:].strip())
            else:
                current.append(line)
        if current:
            sequences.append("".join(current))
    if not sequences:
        raise ValueError(f"No sequences found in {path}")
    if len(headers) < len(sequences):
        headers.extend([""] * (len(sequences) - len(headers)))
    return list(zip(headers, sequences))


def _load_ligand_smiles(path: Path) -> str:
    suffix = path.suffix.lower()
    mol = None
    if suffix in {".sdf", ".sd"}:
        supplier = Chem.SDMolSupplier(str(path), sanitize=True, removeHs=False)
        mol = next((item for item in supplier if item is not None), None)
    elif suffix == ".mol":
        mol = Chem.MolFromMolFile(str(path), sanitize=True, removeHs=False)
    elif suffix == ".mol2":
        mol = Chem.MolFromMol2File(str(path), sanitize=True, removeHs=False)
    elif suffix == ".pdb":
        mol = Chem.MolFromPDBFile(str(path), sanitize=True, removeHs=False)
    else:
        msg = "Ligand file must be .sdf, .mol, .mol2, or .pdb"
        raise ValueError(msg)
    if mol is None:
        raise ValueError(f"Failed to read ligand file {path}")
    smiles = Chem.MolToSmiles(mol)
    if not smiles:
        raise ValueError(f"Failed to generate SMILES from {path}")
    return smiles


def build_yaml_config(
    protein_path: Union[str, Path],
    ligand_path: Optional[Union[str, Path]] = None,
    output_path: Union[str, Path] = "config.yaml",
) -> Path:
    protein_path = Path(protein_path).expanduser()
    output_path = Path(output_path).expanduser()
    ligand_path = Path(ligand_path).expanduser() if ligand_path else None

    entries = []
    used_ids: set[str] = set()
    fallback_index = 0

    for header, sequence in _parse_fasta(protein_path):
        candidate = _extract_chain_id(header)
        if candidate and candidate not in used_ids:
            chain_id = candidate
        else:
            while True:
                chain_id = _chain_id_from_index(fallback_index)
                fallback_index += 1
                if chain_id not in used_ids:
                    break
        used_ids.add(chain_id)
        entries.append({"protein": {"id": chain_id, "sequence": sequence}})

    if ligand_path is not None:
        ligand_id = "L" if "L" not in used_ids else None
        if ligand_id is None:
            while True:
                ligand_id = _chain_id_from_index(fallback_index)
                fallback_index += 1
                if ligand_id not in used_ids:
                    break
        used_ids.add(ligand_id)
        smiles = _load_ligand_smiles(ligand_path)
        entries.append({"ligand": {"id": ligand_id, "smiles": smiles}})

    config = {
        "version": 1,
        "sequences": entries,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        yaml.safe_dump(config, sort_keys=False, default_flow_style=False)
    )
    return output_path
