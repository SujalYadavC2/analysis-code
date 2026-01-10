import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def density_profile(targets:dict, start:int=None, stop:int=None, step:int=1):

    u = list(targets.values())[0].universe

    z_dim_index = 2
    box_z = u.dimensions[z_dim_index]*0.1 #nm
    n_bins = 500
    dz = box_z/n_bins #nm

    bins = np.linspace(0, box_z, n_bins + 1)
    z_centers = ((bins[:-1] + bins[1:]) /2) #nm

    store_hists = {}
    for name, target in targets.items():
        store_hists[name] = np.zeros(n_bins)
    
    n_frames = 0
    for ts in tqdm(u.trajectory[start:stop:step]):

        n_frames +=1

        for name, target in targets.items():
            target.wrap(compound='atoms')
            
            target_z = target.positions[:,z_dim_index]*0.1 #nm
            h_target, _ = np.histogram(target_z, bins=bins)
            store_hists[name] += h_target

    box_x = u.dimensions[0]*0.1 #nm
    box_y = u.dimensions[1]* 0.1 #nm
    bin_vol_nm3 = box_x * box_y * dz


    na = 6.022e23
    molar_factor = 1e27 / na

    store_densities = {}
    for name, hists in store_hists.items():
        store_densities[name] = ((hists/n_frames)/bin_vol_nm3)*molar_factor*1e-3 #mM
    
    #store_wrap_densities = {}
    for name, density in store_densities.items():

        temp_density = np.concatenate((density, density))
        temp_density = temp_density[int(n_bins/2): int(n_bins+n_bins/2)]
        
        np.save(f"{name}_density_profile.npy", np.array([z_centers,temp_density]))
    