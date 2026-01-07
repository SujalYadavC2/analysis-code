import MDAnalysis as mda
from tqdm import tqdm
from pathlib import Path

# Check if the pdf and trj file exist

def file_check(pdb_file:str, traj_file:str) -> bool:
    """
    This functions checks if the file exists or not
    
    :param pdb_file (str): File name or Path
    :param traj_file (str): File name or Path

    Return:
        A bool (True or False)
    """

    pdb_path = Path(pdb_file)
    traj_path = Path(traj_file)

    pdb_check = pdb_path.exists()
    traj_check = traj_path.exists()

    if not pdb_check:
        return False
    
    elif not traj_check:
        return False
    else:
        return True

# First I need to create a function that can convert .dcd file to .xtc

def dcd_to_xtc(pdb_file:str, dcd_file:str) -> None:
    """
    Fuction for converting .dcd file to .xtc
    
    :param pdb_file (str): File name or path
    :param dcd_file (str): File name or path
    """

    # file check
    check = file_check(pdb_file, dcd_file)

    if check:
        u = mda.Universe(pdb_file, dcd_file)
        ag = u.select_atoms('all')

        with mda.Writer('traj.xtc', ag.n_atoms) as W:

            for ts in tqdm(u.trajectory):
                W.write(ag)

        return None
    else:
        print("One of the file does not exits. Please check again.")
        return None


# Check if simulation has produced the number of frames that you were hoping for

def frame_num_check(mda_universe:mda.Universe, expecting_num:int) -> bool:
    """
    Checks the number of frames in the simulation
    
    :param mda_universe: Universe from MDAnalysis
    :type mda_universe: mda.Universe

    :param expecting_num: Number of frames
    :type expecting_num: int
    
    :return: If the simulation has same of number of frames as you were expecting
    :rtype: bool
    """
    
    u = mda_universe
    num_frames = u.trajectory.n_frames

    if num_frames == expecting_num:
        return True
    else:
        return False
