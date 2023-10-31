import argparse
import ast
import qucat as qc
import qutip as qt
import matplotlib.pyplot as plt
import numpy as np
import data_gen_time as dt


def validate_input(N_exc, nmax):
    if not isinstance(N_exc, int):
        raise ValueError("N_exc should be an integer")
    if not isinstance(nmax, int):
        raise ValueError("nmax should be an integer (default = 100)")


def get_input_args():
    parser = argparse.ArgumentParser(description="Process input from a file or command line")
    parser.add_argument("--file", type=str, help="Input filename")
    parser.add_argument("--N_exc", type=int, help="Integer: number of excitations")
    parser.add_argument("--test", action='store_true', help="Show test figs: True or False")
    parser.add_argument("--nmax", type=int, default=100, help="Integer: max number of L_j values")
    parser.add_argument("--version", type=str, default="v1.2", help="Version of the output file")
    return parser.parse_args()


def load_config_from_file(file_name):
    try:
        with open(file_name, 'r') as file:
            data = file.read()
            config = ast.literal_eval(data)  # Safely evaluate the content as a Python literal

            if not isinstance(config, dict):
                raise ValueError("File should contain a dictionary")

            return (config.get("L_js"), config.get("N_exc"),
                    config.get("t_max"), config.get("tau_max"),
                    config.get("test"), config.get("nmax"),
                    config.get("nmax_time"), config.get("time"))
    except FileNotFoundError:
        raise FileNotFoundError(f"Error: File '{file_name}' not found.")
    except ValueError as e:
        raise ValueError(f"Error: {e}")


args = get_input_args()

if args.file:
    N_exc, TEST, nmax, nmax_time, TIME, POWER_LJ = load_config_from_file(args.file)
else:
    N_exc = args.N_exc
    TEST = args.test
    NMAX = args.nmax
    VERSION = args.version

try:
    PATH_DATA, PATH_FIGS = dt.path()
    # Setting the system
    validate_input(N_exc, NMAX)

    cir_1 = qc.GUI('circuits/1_transmon.txt', edit=False, plot=False,
                   print_network=False)
    cir_2 = qc.GUI('circuits/1_transmon.txt', edit=False, plot=False,
                   print_network=False)
    cirs = [cir_1, cir_2]
    N_qubits = len(cirs)
    DIM = N_exc**N_qubits

    a_ops, i_op = dt.operators(N_exc, N_qubits)

    # Power dependent functions
    def transmission_power_dependence(a_in, L_j1, L_j2, direction='R'):
        N_qubits = len(cirs)
        if direction == 'R':
            a_inc = {'R': a_in, 'L': 0.0}
        elif direction == 'L':
            a_inc = {'R': 0.0, 'L': a_in}
        L_js = [L_j1, L_j2]
        p = dt.cir_parameters(cirs, L_js, a_inc, a_ops, i_op)
        out_ops = p['out_ops']
        rho_D, rho_B = dt.dark_bright_states(a_ops, p)

        H, _, _, _ = dt.hamiltonian_transmon(a_ops, p)
        L_diss = sum([sum([p['Gamma'][i, j] * qt.lindblad_dissipator(a_ops[i],
                                                                     a_ops[j])
                           for j in range(N_qubits)]) for i in range(N_qubits)])
        c_ops_nr = [np.sqrt(p["gamma_nr"]) * a_ops[j] for j in range(len(a_ops))]
        c_ops_deph = [np.sqrt(p["gamma_phi"]) * (a_ops[j].dag()*a_ops[j] -
                                                 a_ops[j] * a_ops[j].dag())
                      for j in range(len(a_ops))]
        L_0 = qt.liouvillian(H, c_ops=c_ops_nr + c_ops_deph)
        L = L_0 + L_diss

        rho_ss = qt.steadystate(L, method='direct', maxiter=1e9, tol=1e-24,
                                use_precond=True)
        return (dt.transmission(direction, rho_ss, a_in, out_ops),
                np.real((rho_ss * rho_D).tr()),
                np.real((rho_ss * rho_B).tr()),
                np.real((rho_ss * rho_ss).tr()))  # intensity transmission

    # Power dependence: time independent
    # field amplitudes
    L_js = [3.31e-9, 3.3e-9]
    a_inc_list = 10.0**np.linspace(-6, 1, NMAX)*np.sqrt(1.5e+9)

    # Write to file
    p_for_plot = dt.cir_parameters(cirs, L_js, {'L': 3873.0, 'R': 0.0},
                                   a_ops, i_op)

    # Power L_j dependence
    from IPython.display import clear_output
    L_j_list = np.linspace(3.15e-9, 3.45e-9, NMAX)
    np.insert(L_j_list, NMAX // 2, L_js[1])
    eff = np.zeros((NMAX, NMAX))
    trans_R = np.zeros((NMAX, NMAX))
    trans_L = np.zeros((NMAX, NMAX))
    for j in range(len(L_j_list)):
        tr_R = np.array(qt.parallel_map(transmission_power_dependence,
                                        a_inc_list, task_args=(L_j_list[j],
                                                               L_js[1], ),
                                        task_kwargs=dict(direction='R')))[:, 0]
        tr_L = np.array(qt.parallel_map(transmission_power_dependence,
                                        a_inc_list, task_args=(L_j_list[j],
                                                               L_js[1], ),
                                        task_kwargs=dict(direction='L')))[:, 0]
        eff[j] = (np.maximum(tr_L, tr_R)*np.abs((tr_L - tr_R) /
                                                (tr_L + tr_R)))
        trans_R[j] = tr_R
        trans_L[j] = tr_L
        clear_output(wait=True)  # Clear the output of the current cell
        print(j, ' out of ', len(L_j_list), "\n")

    # Write to file
    dict_power_Lj = {'a_inc_list': a_inc_list, 'L_j_list': L_j_list,
                     'transmission_R': trans_R,
                     'transmission_L': trans_L,
                     'efficiency': eff}
    dt.write_h5(PATH_DATA+"power_Lj_dependence_Nexc" + str(N_exc)
                + "_" + VERSION + ".h5", dict_power_Lj)

    if TEST:
        X = L_j_list / 1e-9
        Y = a_inc_list**2 / np.mean(p_for_plot["gammas"])
        fig_eff, ax = plt.subplots(1, 1, constrained_layout=True)
        cm0 = ax.contourf(X, Y, eff.T, cmap="BuPu", levels=30)
        ax.set_xscale('function',
                      functions=(lambda x: (x)**2,
                                 lambda x: (x)**(1/2)))
        ax.set_yscale('log')
        ax.set_xlabel(r"$L_j$ (nH)")
        ax.set_ylabel(r'$|a_\mathrm{inc}|^2 / \bar{\gamma}$')
        ax.set_xlim(3.22, 3.38)
        ax.set_ylim(1e-6, 1e1)
        fig_eff.colorbar(cm0, ax=ax, label=r"$\mathcal{M}$")

    if TEST:
        plt.show()

except ValueError as e:
    print(f"Error: {e}")
    exit(1)
