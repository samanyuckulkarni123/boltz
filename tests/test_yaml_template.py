from pathlib import Path

import yaml
from rdkit import Chem

from boltz.tools.yaml_template import build_yaml_config


def _write_fasta(path: Path, sequence: str) -> None:
    path.write_text(f">A|protein\n{sequence}\n")


def _write_sdf(path: Path) -> None:
    mol = Chem.MolFromSmiles("CCO")
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()


def test_build_yaml_config_with_ligand(tmp_path: Path) -> None:
    fasta_path = tmp_path / "protein.fasta"
    ligand_path = tmp_path / "ligand.sdf"
    output_path = tmp_path / "config.yaml"

    _write_fasta(fasta_path, "MKT")
    _write_sdf(ligand_path)

    build_yaml_config(fasta_path, ligand_path, output_path)

    assert output_path.exists()
    data = yaml.safe_load(output_path.read_text())
    assert data["version"] == 1
    assert data["sequences"][0]["protein"]["sequence"] == "MKT"
    assert "ligand" in data["sequences"][1]
    assert data["sequences"][1]["ligand"]["smiles"]


def test_build_yaml_config_without_ligand(tmp_path: Path) -> None:
    fasta_path = tmp_path / "protein.fasta"
    output_path = tmp_path / "config.yaml"

    _write_fasta(fasta_path, "MAAA")

    build_yaml_config(fasta_path, None, output_path)

    assert output_path.exists()
    data = yaml.safe_load(output_path.read_text())
    assert data["version"] == 1
    assert len(data["sequences"]) == 1
    assert data["sequences"][0]["protein"]["sequence"] == "MAAA"
