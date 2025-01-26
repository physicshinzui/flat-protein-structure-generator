# Flat Peptide Structure Generation using PyMOL

## Overview
This script generates flat peptide structures using PyMOL's fab command and extracts atomic coordinates into a structured NumPy dictionary. The generated data is stored in a .npz file for further analysis.

## An example entry in the .npz file:
```
{
  'ACDEFGHIK_coordinates': np.array([...]),
  'ACDEFGHIK_atom_names': np.array([...]),
  'ACDEFGHIK_residue_ids': np.array([...]),
  'ACDEFGHIK_residue_names': np.array([...]),
  'sequence_order': np.array([...])
}
```