import sys
import os
sys.path.append(rf"{os.getcwd()}") # This is specific to the VSCode project to run code as modules
from codebase.rematch.memory_rematch import filter_similar_structures
from ase.db import connect

keep_indices = filter_similar_structures(r"C:\Users\User\OneDrive - Cardiff University\Documents\python\matt\cci-ucl-workshop\codebase\rematch\2-layer-set-1.db")
print(f"Keeping {len(keep_indices)} unique structures.")


#MAKE SURE TO SAVE THE KEEP INDICES SOMEWHERE SAFE
with connect("codebase/data/prescreened_structures.db") as db:
    with connect("codebase/data/unique_structures.db") as db_out:
        for i in keep_indices:
            row = db.get(i + 1)  # ASE DB indices are 1-based
            db_out.write(row.toatoms())