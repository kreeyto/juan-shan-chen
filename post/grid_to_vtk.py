from pyevtk.hl import gridToVTK
from get_sim_info import *
import numpy as np

def save_vtk_2d(macrs_dict, path, filename_write, points=True, norm_val=1):
    info = retrieve_sim_info(path)

    if norm_val == 0:
        norm_val = info['NX']
        if points:
            norm_val -= 1

    dx, dy = 1.0 / norm_val, 1.0 / norm_val
    dz = 1.0  

    if info['Prc'] == 'double':
        prc = 'float64'
    else:
        prc = 'float32'

    if not points:
        x = np.arange(0, info['NX'] / norm_val + 0.1 * dx, dx, dtype=prc)
        y = np.arange(0, info['NY'] / norm_val + 0.1 * dy, dy, dtype=prc)
        z = np.array([0], dtype=prc)
        gridToVTK(path + filename_write, x, y, z, cellData=macrs_dict)
    else:
        x = np.arange(0, (info['NX'] - 1) / norm_val + 0.1 * dx, dx, dtype=prc)
        y = np.arange(0, (info['NY'] - 1) / norm_val + 0.1 * dy, dy, dtype=prc)
        z = np.array([0], dtype=prc)
        gridToVTK(path + filename_write, x, y, z, pointData=macrs_dict)
