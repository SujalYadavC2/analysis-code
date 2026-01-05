import MDAnalysis as mda
from tqdm import tqdm

# First I need to create a function that can convert .dcd file to .xtc

def dcd_to_xtc(pdb_file, dcd_file):
    u = mda.Universe(pdb_file, dcd_file)
    ag = u.select_atoms('all')

    with mda.Writer('traj.xtc', ag.n_atoms) as W:

        for ts in tqdm(u.trajectory):
            W.write(ag)

