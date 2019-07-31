# Split a PDB file into protein and named ligands

Simple code that takes a protein PDB, and splits it into a PDB of just the protein, just the excipient and then the protein combined with each excipient.
The ligands must be passed as csv.

## Usage

    pdb_split_protein_and_ligands.py --PDB file.pdb --ligands "NAG", "BMA"

## Output

1. PDB for the protein alone
2. Individual PDBs for each specified excipient
3. A combined PDB of the protein and the individual excipient
