"""Reusable WS8 P1 chemistry, deliberately outside the dependency-free package."""
import hashlib
import json

mcs_params = {
    "timeout": 10, "ringMatchesRingOnly": True, "completeRingsOnly": True,
    "atomCompare": "CompareElements", "bondCompare": "CompareOrder",
    "matchValences": False, "matchChiralTag": False, "heavy_atoms_only": True,
}
mcs_params_sha256 = hashlib.sha256(
    json.dumps(mcs_params, sort_keys=True, separators=(",", ":")).encode()
).hexdigest()


class RejectedStructure(ValueError):
    """The message is a stable, non-identifying rejection reason."""


def standardise(smiles):
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize
    mol = Chem.MolFromSmiles(smiles or "")
    if mol is None or not mol.GetNumAtoms():
        raise RejectedStructure("parse")
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    mol = max(fragments, key=lambda m: (m.GetNumHeavyAtoms(), Chem.MolToSmiles(m)))
    mol = rdMolStandardize.Uncharger().uncharge(mol)
    Chem.RemoveStereochemistry(mol)
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)
    mol = Chem.RemoveHs(mol)
    Chem.SanitizeMol(mol)
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        raise RejectedStructure("radical")
    if any(a.GetSymbol() not in {"C", "H", "N", "O", "P", "S", "F", "Cl", "Br", "I", "Si", "B", "Se"}
           for a in mol.GetAtoms()):
        raise RejectedStructure("element")
    if not 4 <= mol.GetNumHeavyAtoms() <= 80:
        raise RejectedStructure("heavy_atom_range")
    return mol


def hill_formula(mol):
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    return CalcMolFormula(mol)


def first_block(mol):
    from rdkit.Chem import MolToInchiKey
    return MolToInchiKey(mol).split("-")[0]


def mcs_over_set(mols):
    from rdkit import Chem
    from rdkit.Chem import rdFMCS
    members = [Chem.RemoveHs(m) for m in mols]
    if len(members) < 2 or any(not m.GetNumHeavyAtoms() for m in members):
        return "", 0, False
    result = rdFMCS.FindMCS(
        members, timeout=mcs_params["timeout"], ringMatchesRingOnly=True,
        completeRingsOnly=True, atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareOrder, matchValences=False,
        matchChiralTag=False,
    )
    # A canceled search is never a measurement, even if RDKit returns a partial MCS.
    if result.canceled:
        return "", 0, True
    return result.smartsString, result.numAtoms, False


def contained(truth, smarts):
    from rdkit import Chem
    query = Chem.MolFromSmarts(smarts)
    return bool(query and query.GetNumAtoms() and truth.HasSubstructMatch(query, useChirality=False))
