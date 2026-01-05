import MDAnalysis as mda
from tqdm import tqdm
from pathlib import Path

# Check if the pdf and trj file exist

def file_check(pdb_file, traj_file):
    pdb_path = Path(pdb_file)
    traj_path = Path(traj_file)

    pdb_check = pdb_path.exists()
    traj_check = traj_path.exists()

    if not pdb_check:
        print(f"File path {pdb_path} does not exists")
        return False
    
    elif not traj_check:
        print(f"File path {traj_path} does not exists")
        return False
    else:
        return True

# First I need to create a function that can convert .dcd file to .xtc

def dcd_to_xtc(pdb_file, dcd_file):
    u = mda.Universe(pdb_file, dcd_file)
    ag = u.select_atoms('all')

    with mda.Writer('traj.xtc', ag.n_atoms) as W:

        for ts in tqdm(u.trajectory):
            W.write(ag)

