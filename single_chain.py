import MDAnalysis as mda
import numpy as np
from tqdm import tqdm

# Radius of Gyration

def Rg(chain:mda.AtomGroup, start:int=None, stop:int=None, step:int=1) -> list:
    """
    Calculates radius of gyration over the trajectory for a chain
    
    :param chain: molecular chain
    :type chain: mda.AtomGroup

    :param start: start of the traj
    :type start: int

    :param stop: end of the traj
    :type stop: int

    :param step: step to take over traj
    :type step: int

    :return: A list with first index as np.array of time and second index np.array of radius of gyration
    :rtype: list
    """
    u = chain.universe

    rg_list = []

    for ts in tqdm(u.trajectory[start:stop:step]):
        rg = chain.radius_of_gyration()
        rg_list.append(rg)

    time = np.array([ts.time for ts in u.trajectory[start:stop:step]])

    rg_list = np.array(rg_list)

    return [time, rg_list]


# End-to-End Distance
def Ree(chain:mda.AtomGroup, start:int=None, stop:int=None, step:int=1) -> list:
    """
    Calculate end-to-end distance over trajectory
    
    :param chain: molecular chain
    :type chain: mda.AtomGroup

    :param start: start of the traj
    :type start: int

    :param stop: end of the traj
    :type stop: int

    :param step: step to take over traj
    :type step: int

    :return: A list with first index as np.array of time and second index np.array of end-to-end distance
    :rtype: list
    """

    u = chain.universe

    ree_list = []
    for ts in u.trajectory[start:stop:step]:
        mda.lib.nsgrid.FastNS(u.dimensions).search(chain.positions)

        pos_first = chain.positions[0]
        pos_last = chain.positions[-1]

        dist = np.linalg.norm(pos_last-pos_first)
        ree_list.append(dist)
    
    time = np.array([ts.time for ts in u.trajectory[start:stop:step]])

    ree_list = np.array(ree_list)

    return [time, ree_list]
