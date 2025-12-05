from ase.db import connect
import sys
import os
from ase.io import read
import argparse

parser = argparse.ArgumentParser(description="Populate a database with POSCAR files in the current directory.")
parser.add_argument('--output', type=str, default='structures.db', help='Name of the database file to create or append to.')
parser.add_argument('--input', type=str, default='.', help='Directory containing POSCAR files.')
parser.add_argument('--energies', type=str, default='energies.csv', help='CSV file containing energies for the structures.')
args = parser.parse_args()

db_name = args.output
path = args.input
energy_file = args.energies
with connect(db_name) as db:
    list = [f for f in os.listdir(path) if f.startswith('POSCAR')]
    for name in list:
        db.write(read(f"{path}/{name}"), filename=name, format='vasp')


import csv
energy_dict = {}   # filename → energy
with open(energy_file) as f:
    reader = csv.reader(f)
    next(reader)   # skip header
    for filename, energy in reader:
        energy_dict[filename] = float(energy)

with connect(db_name) as db:
    for row in db.select():
        fname = row.filename

        if fname in energy_dict:
            db.update(row.id, energy_eV=energy_dict[fname])
            print(f"Updated {fname} with energy {energy_dict[fname]}")