from ase.io import read
import sys
import os
sys.path.append(rf"{os.getcwd()}") # This is specific to the VSCode project to run code as modules
from ase.db import connect
import numpy as np
from dscribe.kernels import REMatchKernel
from sklearn.preprocessing import normalize
import argparse
'''
parser = argparse.ArgumentParser(description="similarity testing.")
parser.add_argument('--db', type=str, default='structures.db', help='Name of the database file.')
parser.add_argument('--csv', type=str, default='output.csv', help='output csv filename.')
args = parser.parse_args()

db_name = args.db
output_csv = args.csv
'''
db_name = "codebase/2-layer-set-1.db"
output_csv = "similarity_test_output.csv"

_kernel_cache = {}

def _get_kernel(gamma=1.0, alpha=0.01, threshold=1e-6):
    """
    Reuse REMatch kernels for identical hyperparameters since constructing them
    repeatedly is relatively expensive.
    """
    key = (gamma, alpha, threshold)
    kernel = _kernel_cache.get(key)
    if kernel is None:
        kernel = REMatchKernel(metric="rbf", gamma=gamma, alpha=alpha, threshold=threshold)
        _kernel_cache[key] = kernel
    return kernel

def get_soap(atoms):
    from dscribe.descriptors import SOAP

    species = ["Ce", "O", "Au", "Cu"]
    r_cut = 6.0
    n_max = 8
    l_max = 6

    # Setting up the SOAP descriptor
    soap = SOAP(
        species=species,
        periodic=True,
        r_cut=r_cut,
        n_max=n_max,
        l_max=l_max,
    )

    return soap.create(atoms, n_jobs=1)

def soap_rematch_similarity(soap1, soap2, gamma=1.0, alpha=0.01, threshold=1e-6, kernel=None):
    """
    input two SOAP, return REMatch
    """
    soap1_norm = normalize(soap1)
    soap2_norm = normalize(soap2)

    if kernel is None:
        kernel = _get_kernel(gamma=gamma, alpha=alpha, threshold=threshold)
    kernel_matrix = kernel.create([soap1_norm, soap2_norm])

    return kernel_matrix[0, 1]

def shave_slab(atoms, threshold=3.0, fix=["Ce", "O"]):
    ''' Remove atoms below a certain z-coordinate threshold, only for specified species to fix
    '''
    positions = atoms.get_positions()
    oxide_z_coordinate = float(np.max([[position[2] for _, position in enumerate(positions) if atoms[_].symbol in fix]]))

    threshold = oxide_z_coordinate - threshold
    symbols = np.array(atoms.get_chemical_symbols())
    fix_mask = np.isin(symbols, fix)
    bottom_mask = positions[:, 2] <= threshold
    keep_indices = np.nonzero(~(fix_mask & bottom_mask))[0].tolist()
    if not keep_indices:
        return
    temp_atoms = atoms.copy()[keep_indices]
    return keep_indices, temp_atoms


with connect(db_name) as db:
    clusters = [row.toatoms() for row in db.select(selection="id<3000")]
    files = [row.filename for row in db.select(selection="id<3000")]
    energies = [row.energy_eV for row in db.select(selection="id<3000")]

#find the lowest energy structure
idx_min = np.argmin(energies)
'''
# shaved slab similarities
clusters_shaved = [shave_slab(cluster, threshold=3)[1] for cluster in clusters]

soaps = [get_soap(cluster) for cluster in clusters_shaved]
similarities_shaved = [soap_rematch_similarity(soaps[idx_min], soap) for soap in soaps]
similarities_shaved = [round(sim, 6) for sim in similarities_shaved]
'''

# shaved slab similarities
clusters_cluster_only = [shave_slab(cluster, threshold=-1)[1] for cluster in clusters]

soaps = [get_soap(cluster) for cluster in clusters_cluster_only]

soaps = [get_soap(cluster) for cluster in clusters_cluster_only]
similarities_cluster_only = [soap_rematch_similarity(soaps[idx_min], soap) for soap in soaps]
similarities_cluster_only = [round(sim, 6) for sim in similarities_cluster_only]

'''
import csv

rows = zip(files, energies, similarities_shaved, similarities_cluster_only)

with open(output_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["filename", "energy (eV)", "sim_shaved", "sim_cluster_only"])  # header
    writer.writerows(rows)


import csv
import numpy as np

rows = list(zip(files, energies, similarities_cluster_only))
rows_sorted = sorted(rows, key=lambda x: x[1])


sim_threshold = 0.99      
energy_threshold = 0.2 

filtered_rows = []

for row in rows_sorted:
    fname, E, sim_cluster = row

    is_duplicate = False

    for kept in filtered_rows:
        fname_k, E_k, sim_cluster_k = kept

        # Check if within energy window
        if abs(E - E_k) <= energy_threshold:
            # Check if similarity indicates same structure
            if round(sim_cluster, 2) == round(sim_cluster_k, 2):
                is_duplicate = True
                break

    if not is_duplicate:
        filtered_rows.append(row)

with open(output_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["filename", "energy (eV)", "sim_cluster_only"])
    writer.writerows(filtered_rows)

new_db = "unique_structures.db"

from ase.db import connect

with connect(new_db) as db_new:
    for atoms, fname, E, sim in filtered_rows:
        db_new.write(
            atoms,
            filename=fname,
            energy_eV=E,
            sim_cluster_only=sim
        )

'''
import csv
import numpy as np
from ase.db import connect

# ----------------------------
# Build raw rows
# ----------------------------
rows = list(zip(clusters, files, energies, similarities_cluster_only))
rows_sorted = sorted(rows, key=lambda x: x[2])   # x[2] = energy

sim_threshold = 0.99      
energy_threshold = 0.2 

filtered_rows = []
unique_structures = []   # keep atoms + metadata for DB

for atoms, fname, E, sim_cluster in rows_sorted:

    is_duplicate = False

    for kept_atoms, fname_k, E_k, sim_cluster_k in unique_structures:

        # Check if within energy window
        if abs(E - E_k) <= energy_threshold:
            # Check if similarity indicates same structure
            if round(sim_cluster, 2) == round(sim_cluster_k, 2):
                is_duplicate = True
                break

    if not is_duplicate:
        filtered_rows.append((fname, E, sim_cluster))
        unique_structures.append((atoms, fname, E, sim_cluster))

# ----------------------------
# Save filtered CSV
# ----------------------------
with open(output_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["filename", "energy (eV)", "sim_cluster_only"])
    writer.writerows(filtered_rows)

# ----------------------------
# Save unique structures to new ASE DB
# ----------------------------
new_db = "unique_structures.db"

with connect(new_db) as db_new:
    for atoms, fname, E, sim in unique_structures:
        db_new.write(
            atoms,
            filename=fname,
            energy_eV=E,
            sim_cluster_only=sim
        )
