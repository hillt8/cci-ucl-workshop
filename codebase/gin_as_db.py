import sys
sys.path.append(r"C:\Users\c1528354\GitHub\ucl-cci\cci-ucl-workshop")
from codebase.io.klmc_parser import gin_files_to_db
from codebase.rematch.prescreen import match_cluster_size

# gin_files_to_db("codebase/data")

from ase.db import connect

with connect("codebase/data/structures.db") as db:
    with connect("codebase/data/prescreened_structures.db") as db_out:
        for row in db.select():
            is_single_cluster = match_cluster_size(slab=row.toatoms(), size=12, species=["Au"])
            print(f"ID: {row.id}, Formula: {row.formula}, Single cluster {is_single_cluster}")
            
            if is_single_cluster:
                db_out.write(row.toatoms())
           
from ase.visualize import view
from ase.io import read
view(read("codebase/data/prescreened_structures.db@:"))