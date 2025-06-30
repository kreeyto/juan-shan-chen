import os
import glob
import numpy as np

__macr_names__ = ['rho_water', 'rho_air', 'u1_water', 'u1_air']
__info__ = dict()

def get_filenames_macr(macr_name, path):
    return sorted(glob.glob(path + "/*" + macr_name + "*.bin"))

def get_macr_steps(path):
    file_list = get_filenames_macr(__macr_names__[0], path)
    step_set = set()
    for file in file_list:
        step_str = file.split(__macr_names__[0])[-1]
        step_str = step_str[:-4]  # remove ".bin"
        step_set.add(int(step_str))
    return sorted(step_set)

def retrieve_sim_info(path):
    if len(__info__) == 0:
        filename = glob.glob(path + "/*info*.txt")[0]
        with open(filename, "r") as f:
            lines = [line.strip() for line in f.readlines()]

            def extract(keyword, cast=str):
                try:
                    return cast([line.split()[-1] for line in lines if keyword in line][0])
                except Exception:
                    print(f"Could not extract '{keyword}' from info file.")
                    return None

            __info__['ID'] = extract('Simulation ID')
            __info__['Prc'] = extract('Precision')
            __info__['NX'] = extract('NX', int)
            __info__['NY'] = extract('NY', int)
            __info__['Tau'] = extract('Tau', float)
            __info__['Umax'] = extract('Umax', float)
            __info__['Nsteps'] = extract('Nsteps', int)
    return __info__

def read_file_macr_2d(macr_filename, path):
    info = retrieve_sim_info(path)
    dtype = 'd' if info['Prc'] == 'double' else 'f'
    with open(macr_filename, "rb") as f:
        vec = np.fromfile(f, dtype)
        vec_2d = np.reshape(vec, (info['NY'], info['NX']), order='C')
        return vec_2d

def get_macrs_from_step(step, path):
    macr = dict()
    for macr_name in __macr_names__:
        filename = os.path.join(path, f"{__info__['ID']}_{macr_name}{step:06d}.bin")
        if os.path.exists(filename):
            macr[macr_name] = read_file_macr_2d(filename, path)
    return macr if macr else None

def get_all_macrs(path):
    macr = dict()
    filenames = {macr_name: get_filenames_macr(macr_name, path) for macr_name in __macr_names__}
    min_length = min(len(v) for v in filenames.values())

    for i in range(min_length):
        step_str = filenames[__macr_names__[0]][i].split(__macr_names__[0])[-1][:-4]
        step = int(step_str)
        macr[step] = {macr_name: read_file_macr_2d(filenames[macr_name][i], path)
                      for macr_name in __macr_names__}
    return macr
