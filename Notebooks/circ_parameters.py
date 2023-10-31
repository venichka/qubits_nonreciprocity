import argparse
import ast
import qucat as qc
import numpy as np
import h5py

def write_h5(file_name, data_dict):
    with h5py.File(file_name, "w") as f:
        for dset_name, dset_data in data_dict.items():
            f.create_dataset(dset_name, data=dset_data)
            print(dset_name, dset_data)
    return "File .h5 wrote successfully"

def read_h5(file_name):
    with h5py.File(file_name, "r") as f:
        dict_data = {dset_name: f[dset_name][:] for dset_name in f.keys()}
        for dset_name, dset_data in dict_data.items():
            print(dset_name, dset_data)
    return dict_data

def validate_input(L_js):
    if not isinstance(L_js, list) or not all(isinstance(x, float) for x in L_js):
        raise ValueError("L_js should be a list of floats")

def get_input_args():
    parser = argparse.ArgumentParser(description="Process input from a file or command line")
    parser.add_argument("--file", type=str, help="Input filename")
    parser.add_argument("--L_js", type=float, nargs='+', help="List of floats")
    return parser.parse_args()

def load_config_from_file(file_name):
    try:
        with open(file_name, 'r') as file:
            data = file.read()
            config = ast.literal_eval(data)  # Safely evaluate the content as a Python literal

            if not isinstance(config, dict):
                raise ValueError("File should contain a dictionary")

            return config.get("L_js")
    except FileNotFoundError:
        raise FileNotFoundError(f"Error: File '{file_name}' not found.")
    except ValueError as e:
        raise ValueError(f"Error: {e}")

def main():
    args = get_input_args()

    if args.file:
        L_js = load_config_from_file(args.file)
    else:
        L_js = args.L_js

    try:
        validate_input(L_js)
        
        DIRECTION = 'L'

        a_mag = 0.1 * np.sqrt(1.5e+9)
        a_R = a_mag if DIRECTION == 'R' else 0e-2
        a_L = a_mag if DIRECTION == 'L' else 0e-2
        a_inc = {'R': a_R, 'L': a_L}

        cir_1 = qc.GUI('circuits/1_transmon.txt', edit=False, plot=False, print_network=True)
        cir_2 = qc.GUI('circuits/1_transmon.txt', edit=False, plot=False, print_network=True)
        cirs = [cir_1, cir_2]

        p = cir_parameters(cirs, L_js)

        FILENAME = "../Data/params.h5"
        write_h5(FILENAME, p)

        # Now, you can use L_js in your script
        print("L_js:", L_js)
    except ValueError as e:
        print(f"Error: {e}")
        exit(1)

def cir_parameters(cirs, L_js):
    gammas = 2.0 * np.pi * np.array([cirs[i].loss_rates(Lj=L_js[i])[0] for i in range(len(cirs))])
    w_m = 2.0 * np.pi * np.array([cirs[i].eigenfrequencies(Lj=L_js[i])[0] for i in range(len(cirs))])
    w_ext = w_m[0]
    delta = (w_m[0] - w_m[1]) / np.mean(gammas)
    t_j = np.array([0.0, (np.pi - delta) / w_ext])
    t_ij = np.array([[t_j[i] - t_j[j] for j in range(len(cirs))] for i in range(len(cirs))])
    
    return {
        'w_m': w_m,
        'w_ext': w_ext,
        'A_m': 2.0 * np.pi * np.array([cirs[i].anharmonicities(Lj=L_js[i])[0] for i in range(len(cirs))]),
        'gammas': gammas,
        'gamma_nr': 190.0e3 * 2.0 * np.pi,
        'gamma_phi': 200.0e3 * 2.0 * np.pi,
        'delta': delta,
        't_j': t_j,
        't_ij': t_ij,
        'Omega': np.array([[0.5 * np.sqrt(gammas[i] * gammas[j]) *
                            np.sin(w_ext * np.abs(t_ij[i, j]))
                            for j in range(len(cirs))]
                           for i in range(len(cirs))]),
        'Gamma': np.array([[np.sqrt(gammas[i] * gammas[j]) *
                            np.cos(w_ext * np.abs(t_ij[i, j]))
                            for j in range(len(cirs))]
                           for i in range(len(cirs))]),
    }

if __name__ == "__main__":
    main()
