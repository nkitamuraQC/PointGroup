import numpy as np
import spglib

lattice = np.array([[0.0, 0.5, 0.5],
                    [0.5, 0.0, 0.5],
                    [0.5, 0.5, 0.0]]) * 5.4
positions = [[0.875, 0.875, 0.875],
            [0.125, 0.125, 0.125]]
numbers= [1,] * 2
cell = (lattice, positions, numbers)
# print("space group:", spglib.get_spacegroup(cell, symprec=1e-5))
mesh = [8, 8, 8]
print("symmetry:", spglib.spglib.get_symmetry(cell))
# print("symmetry_dataset:", spglib.spglib.get_symmetry_dataset(cell))
