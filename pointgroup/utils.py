from ase.io import write
from ase.data import chemical_symbols
from ase import Atoms

def get_cif(cell, pos, species, name="vis.cif"):
    symbols = [chemical_symbols[Z] if isinstance(Z, int) else Z for Z in species]
    atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)
    write(name, atoms, format='cif')
    return
