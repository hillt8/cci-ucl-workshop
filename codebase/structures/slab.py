import numpy as np
import sys
import os
sys.path.append(rf"{os.getcwd()}")
from codebase.io.klmc_parser import gin_files_to_db
from codebase.rematch.prescreen import match_cluster_size 

def get_random_samples(arr, size=10):
    return arr[np.random.randint(0, arr.shape[0], size=10)]

from ase.io import read
metal_oxide = read(r'codebase\data\prescreened_structures.db', index=0)
metal_oxide = metal_oxide[[atom.index for atom in metal_oxide if atom.symbol in ['Ce', 'O']]]
metal_oxide.wrap()

from ase.visualize import view
from ase.build import fcc100
overslab = fcc100('Au', size=(4,4,3), vacuum=0)
for atom in overslab:
    atom.position[0] += 3.0
    atom.position[1] += 3.0
    atom.position[2] += 9.5
metal_oxide += overslab
metal_oxide.center(vacuum=10.0, axis=2)
# view(metal_oxide)

indices = np.array([atom.index for atom in metal_oxide if atom.symbol in ['Au']])
indices_metal_oxide = np.array([atom.index for atom in metal_oxide if atom.symbol in ['Ce', 'O']])

random_indices = get_random_samples(indices, size=10)

target_size = 12
atoms_list = []
for i in range(10000):
    random_indices = get_random_samples(indices, size=target_size)
    cluster = metal_oxide[random_indices].copy()
    if match_cluster_size(slab=cluster, size=target_size, species=["Au"]):
        atoms_list.append(cluster)

view(atoms_list)

"""
# Recreate full slab with metal oxide
full_slab = np.append(random_indices, indices_metal_oxide)
cluster = metal_oxide[full_slab]
view(cluster)
"""