from ase.io import read
from ase.visualize import view
from ase.build import bulk
from ase.build import molecule
from ase.build import fcc100, add_adsorbate

"""
# Create a bulk structure with a defect
atoms = bulk('Au', 'fcc', a=4.08, cubic=True)
atoms = atoms * (2, 2, 2)
atoms[0].symbol = 'Ag'  # introduce a defect
"""
# Create a surface with an adsorbate
atoms = fcc100('Cu', size=(3, 3, 3), vacuum=10.0)
adsorbate = molecule('H2O')
adsorbate.rotate(180, 'x')
add_adsorbate(atoms, adsorbate, height=2.5, position='ontop')

view(atoms)

# Calculate the energy of each structure and save in a database
# for image in images:
from mace.calculators import mace_mp
from ase.optimize import BFGS

PARAMS = {'model': "medium",
          'dispersion': True,
          'default_dtype': 'float64',
          'device': 'cpu'}

samples = [(atoms, 0)]

def evaluate_structure(atoms, index):
    atoms.calc = mace_mp(**PARAMS)
    opt = BFGS(atoms, trajectory=f"{index}.traj", logfile=f"{index}.log")
    opt.run(fmax=0.05)

evaluate_structure(atoms, 0)

from ase.io import read
view(read("0.traj@:"))


"""
from multiprocessing import Pool
with Pool(2) as p:
    p.starmap(evaluate_structure, images)
"""
