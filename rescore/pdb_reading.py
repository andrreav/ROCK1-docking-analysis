from rdkit import Chem


def protein_sanitizer(protein: Chem.Mol):
    Chem.SanitizeMol(protein, Chem.SANITIZE_SETAROMATICITY)
    Chem.SanitizeMol(protein, Chem.SANITIZE_SETCONJUGATION)
    Chem.SanitizeMol(protein, Chem.SANITIZE_SETHYBRIDIZATION)
    Chem.SanitizeMol(protein, Chem.SANITIZE_SYMMRINGS)
    # Chem.SanitizeMol(prot, Chem.SANITIZE_PROPERTIES)
    Chem.SanitizeMol(protein, Chem.SANITIZE_CLEANUP)
    # prot = Chem.AddHs(prot, addCoords=True)
    return protein


def read_pdb_block(pdb_block):
    prot = Chem.MolFromPDBBlock(
        pdb_block, sanitize=False, removeHs=False, proximityBonding=True
    )
    assert prot is not None, "Could not read protein from PDB block: \n" + pdb_block
    prot = protein_sanitizer(prot)
    return prot


def read_pdb_file(pdb_file):
    prot = Chem.MolFromPDBFile(
        pdb_file, removeHs=False, sanitize=False, proximityBonding=True
    )
    assert prot is not None, "Could not read protein from PDB file: " + pdb_file
    prot = protein_sanitizer(prot)
    return prot
