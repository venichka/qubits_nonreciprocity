import argparse
import ast
import qucat as qc
import qutip as qt
import scipy as sc
import matplotlib.pyplot as plt
import numpy as np
import data_gen_time as dt


def validate_input(L_js, t_max, tau_max, nmax, nmax_time):
    if not isinstance(L_js, list) or not all(isinstance(x, float) for x in L_js):
        raise ValueError("L_js should be a list of floats")
    if not isinstance(t_max, float):
        raise ValueError("t_max of t_D should be a float (usually 1e1)")
    if not isinstance(tau_max, float):
        raise ValueError("tau_max of t_D should be a float (usually 1e1 or 5e1)")
    if not isinstance(nmax, int):
        raise ValueError("nmax should be an integer (default = 100)")
    if not isinstance(nmax_time, int):
        raise ValueError("nmax_time should be an integer (default = 20)")


def get_input_args():
    parser = argparse.ArgumentParser(description="Process input from a file or command line")
    parser.add_argument("--file", type=str, help="Input filename")
    parser.add_argument("--L_js", type=float, nargs='+', help="List of floats")
    parser.add_argument("--t_max", type=float, help="Integration time (number of t_D): float")
    parser.add_argument("--tau_max", type=float, help="Integration time for g1 and g2 (number of t_D): float")
    parser.add_argument("--test", action='store_true', help="Show test figs: True or False")
    parser.add_argument("--nmax", type=int, default=100, help="Integer: max number of L_j values")
    parser.add_argument("--nmax_time", type=int, default=20, help="Integer: max number of L_j values for time dependent H")
    parser.add_argument("--time", action='store_true', help="Compute time dependent H: True or False")
    parser.add_argument("--power_lj", action='store_true', help="Compute power, L_j dependence: True or False")
    return parser.parse_args()


def load_config_from_file(file_name):
    try:
        with open(file_name, 'r') as file:
            data = file.read()
            config = ast.literal_eval(data)  # Safely evaluate the content as a Python literal

            if not isinstance(config, dict):
                raise ValueError("File should contain a dictionary")

            return (config.get("L_js"),
                    config.get("t_max"), config.get("tau_max"),
                    config.get("test"), config.get("nmax"),
                    config.get("nmax_time"), config.get("time"))
    except FileNotFoundError:
        raise FileNotFoundError(f"Error: File '{file_name}' not found.")
    except ValueError as e:
        raise ValueError(f"Error: {e}")


def DB_basis(a_ops, args, RWA=True):
    N_exc = a_ops[0].dims[0][0]
    gg = qt.basis(N_exc, 0)  # ground state
    ee = qt.basis(N_exc, N_exc - 1)  # excited state
    for j in range(N_qubits - 1):
        gg = qt.tensor(gg, qt.basis(N_exc, 0))
        ee = qt.tensor(ee, qt.basis(N_exc, N_exc - 1))
    if RWA:
        Gamma = args["Gamma"]
    else:
        Gamma = args["Gamma_nonRWA"]
    J_ops = []
    lam, b = sc.linalg.eig(Gamma)
    J_ops = [sum([b[j, i] * a_ops[j] for j in range(N_qubits)])
             for i in range(N_qubits)]
    # dark and bright operators
    for j in range(len(lam)):
        if lam[j] == np.min(lam):
            J_D = J_ops[j]
        elif lam[j] == np.max(lam):
            J_B = J_ops[j]
    dd = J_D.dag() * gg
    bb = J_B.dag() * gg
    basis = np.eye(N_exc**2, dtype=complex)
    basis[:, 1] = ((dd.data).toarray()).flatten()
    basis[:, N_exc] = ((bb.data).toarray()).flatten()
    return basis, gg, dd, bb, ee


def DB_basis_full(a_ops, args, RWA=True):
    N_exc = a_ops[0].dims[0][0]
    gg = qt.basis(N_exc, 0)  # ground state
    ee = qt.basis(N_exc, 1)  # excited state
    for j in range(N_qubits - 1):
        gg = qt.tensor(gg, qt.basis(N_exc, 0))
        ee = qt.tensor(ee, qt.basis(N_exc, 1))
    if RWA:
        Gamma = args["Gamma"]
    else:
        Gamma = args["Gamma_nonRWA"]
    J_ops = []
    lam, b = sc.linalg.eig(Gamma)
    J_ops = [sum([b[j, i] * a_ops[j] for j in range(N_qubits)])
             for i in range(N_qubits)]
    # dark and bright operators
    for j in range(len(lam)):
        if lam[j] == np.min(lam):
            J_D = J_ops[j]
        elif lam[j] == np.max(lam):
            J_B = J_ops[j]
    dd = J_D.dag() * gg
    bb = J_B.dag() * gg
    dd2 = (dd.dag()*qt.basis([N_exc, N_exc], [1, 0]) *
           qt.basis([N_exc, N_exc], [2, 0]) +
           dd.dag()*qt.basis([N_exc, N_exc], [0, 1]) *
           qt.basis([N_exc, N_exc], [0, 2]))  # J_D.dag() * J_D.dag() * gg
    bb2 = (bb.dag()*qt.basis([N_exc, N_exc], [1, 0]) *
           qt.basis([N_exc, N_exc], [2, 0]) +
           bb.dag()*qt.basis([N_exc, N_exc], [0, 1]) *
           qt.basis([N_exc, N_exc], [0, 2]))  # J_B.dag() * J_B.dag() * gg
    basis = np.eye(N_exc**2, dtype=complex)
    basis[:, 1] = ((dd.data).toarray()).flatten()
    basis[:, N_exc] = ((bb.data).toarray()).flatten()
    if N_exc > 2:
        basis[:, 2] = ((dd2.data).toarray()).flatten()
        basis[:, 2*N_exc] = ((bb2.data).toarray()).flatten()
    return basis, gg, dd, bb, ee, dd2, bb2


args = get_input_args()

if args.file:
    L_js, t_max, tau_max, TEST, nmax, nmax_time, TIME, POWER_LJ = load_config_from_file(args.file)
else:
    L_js = args.L_js
    t_max = args.t_max
    tau_max = args.tau_max
    TEST = args.test
    NMAX = args.nmax
    NMAX_TIME = args.nmax_time
    TIME = args.time
    POWER_LJ = args.power_lj

try:
    PATH_DATA, PATH_FIGS = dt.path()
    # Setting the system
    validate_input(L_js, t_max, tau_max, NMAX, NMAX_TIME)

    cir_1 = qc.GUI('circuits/1_transmon.txt', edit=False, plot=False,
                   print_network=False)
    cir_2 = qc.GUI('circuits/1_transmon.txt', edit=False, plot=False,
                   print_network=False)
    cirs = [cir_1, cir_2]
    N_qubits = len(cirs)
    a_ops, i_op = dt.operators(3, N_qubits)

    # Power dependent functions
    def transmission_power_dependence_Nexc(a_in, L_j1, L_j2, N_exc,
                                           direction='R'):
        N_qubits = len(cirs)
        a_ops, i_op = dt.operators(N_exc, N_qubits)
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
        if N_exc == 2:
            basis, _, _, _, _ = DB_basis(a_ops, p, RWA=True)
        elif N_exc > 2:
            basis, _, _, _, _, _, _ = DB_basis_full(a_ops, p, RWA=True)
        return (dt.transmission(direction, rho_ss, a_in, out_ops),
                np.real((rho_ss * rho_D).tr()),
                np.real((rho_ss * rho_B).tr()),
                np.real((rho_ss * rho_ss).tr()),
                (rho_ss.data).toarray(),
                (rho_ss.transform(basis).data).toarray())  # intensity transmission

    def transmission_power_dependence_time_dependent(a_in, L_j1, L_j2, N_exc,
                                                     direction: str = 'R'):
        N_qubits = len(cirs)
        a_ops, i_op = dt.operators(N_exc, N_qubits)
        if direction == 'R':
            a_inc = {'R': a_in, 'L': 0.0}
        elif direction == 'L':
            a_inc = {'R': 0.0, 'L': a_in}
        L_js = [L_j1, L_j2]
        p = dt.cir_parameters(cirs, L_js, a_inc, a_ops, i_op)
        out_ops = p['out_ops']
        rho_D, rho_B = dt.dark_bright_states_nonRWA(a_ops, p)

        H_t, H_args = dt.hamiltonian_transmon_time_drive(a_ops, p)
        J_ops = dt.jump_operators(a_ops, p)
        c_ops_nr = [np.sqrt(p["gamma_nr"]) * a_ops[j] for j in range(len(a_ops))]
        c_ops_deph = [np.sqrt(p["gamma_phi"]) * (a_ops[j].dag()*a_ops[j] -
                                                 a_ops[j]*a_ops[j].dag())
                      for j in range(len(a_ops))]
        # Master equation
        rho0 = dt.initial_state(N_qubits, N_exc)
        w = np.real(sc.linalg.eigvals(p['Gamma_nonRWA']))
        tlist = np.linspace(0.0, t_max/np.min(w), 2)
        result_tt = qt.mesolve(H_t, rho0, tlist,
                               c_ops=J_ops + c_ops_nr + c_ops_deph,
                               e_ops=[], args=H_args,
                               options=qt.Options(nsteps=1e7, rhs_reuse=False))
        rho_ss = result_tt.states[-1]
        basis, _, _, _, _ = DB_basis(a_ops, p, RWA=False)

        return (dt.transmission(direction, rho_ss, a_in, out_ops),
                np.real((rho_ss * rho_D).tr()),
                np.real((rho_ss * rho_B).tr()),
                np.real((rho_ss * rho_ss).tr()),  # intensity transmission
                (rho_ss.data).toarray(),
                (rho_ss.transform(basis).data).toarray())  # intensity transmission

    # Power and N_exc dependence: time independent
    a_inc_list = 10.0**np.linspace(-6, 1, NMAX)*np.sqrt(1.5e+9)
    N_exc_list = np.arange(2, 9)
    data_dict = {}

    for j in N_exc_list:
        result_R = qt.parallel_map(transmission_power_dependence_Nexc,
                                   a_inc_list, task_args=(L_js[0],
                                                          L_js[1],
                                                          j,
                                                          ),
                                   task_kwargs=dict(direction='R'),
                                   )
        result_L = qt.parallel_map(transmission_power_dependence_Nexc,
                                   a_inc_list, task_args=(L_js[0],
                                                          L_js[1],
                                                          j,
                                                          ),
                                   task_kwargs=dict(direction='L'),
                                   )
        transmission_R = np.array([result_R[i][0]
                                   for i in range(len(a_inc_list))])
        pop_D_R = np.array([result_R[i][1]
                            for i in range(len(a_inc_list))])
        pop_B_R = np.array([result_R[i][2]
                            for i in range(len(a_inc_list))])
        rho2_R = np.array([result_R[i][3]
                           for i in range(len(a_inc_list))])
        rhoss_R = np.array([result_R[i][4]
                            for i in range(len(a_inc_list))])
        rhoss_DB_R = np.array([result_R[i][5]
                               for i in range(len(a_inc_list))])
        transmission_L = np.array([result_L[i][0]
                                   for i in range(len(a_inc_list))])
        pop_D_L = np.array([result_L[i][1]
                            for i in range(len(a_inc_list))])
        pop_B_L = np.array([result_L[i][2]
                            for i in range(len(a_inc_list))])
        rho2_L = np.array([result_L[i][3]
                           for i in range(len(a_inc_list))])
        rhoss_L = np.array([result_L[i][4]
                            for i in range(len(a_inc_list))])
        rhoss_DB_L = np.array([result_L[i][5]
                               for i in range(len(a_inc_list))])
        efficiency = (np.maximum(transmission_L, transmission_R) *
                      np.abs((transmission_L - transmission_R) /
                             (transmission_L + transmission_R)))
        data_dict.update({
            'transmission_R_'+str(j): transmission_R,
            'transmission_L_'+str(j): transmission_L,
            'population_B_R_'+str(j): pop_B_R,
            'population_B_L_'+str(j): pop_B_L,
            'population_D_R_'+str(j): pop_D_R,
            'population_D_L_'+str(j): pop_D_L,
            'rho2_R_'+str(j): rho2_R,
            'rho2_L_'+str(j): rho2_L,
            'rho_ss_R_'+str(j): rhoss_R,
            'rho_ss_L_'+str(j): rhoss_L,
            'rho_ss_DB_R_'+str(j): rhoss_DB_R,
            'rho_ss_DB_L_'+str(j): rhoss_DB_L,
            'efficiency_'+str(j): efficiency
            })
        print('N_exc = ', j, '\n')

    # Time dependent
    if TIME:
        # TODO
        a_inc_list_t = 10.0**np.linspace(-6, 1, NMAX_TIME)*np.sqrt(1.5e+9)

        result_R_t = np.array(qt.parallel_map(
            transmission_power_dependence_time_dependent,
            a_inc_list_t, task_args=(L_js[0], L_js[1], ),
            task_kwargs=dict(direction='R'),
            progress_bar=True))
        result_L_t = np.array(qt.parallel_map(
            transmission_power_dependence_time_dependent,
            a_inc_list_t, task_args=(L_js[0], L_js[1], ),
            task_kwargs=dict(direction='L'),
            progress_bar=True))
        transmission_R_t = result_R_t[:, 0]
        pop_D_R_t = result_R_t[:, 1]
        pop_B_R_t = result_R_t[:, 2]
        rho2_R_t = result_R_t[:, 3]
        transmission_L_t = result_L_t[:, 0]
        pop_D_L_t = result_L_t[:, 1]
        pop_B_L_t = result_L_t[:, 2]
        rho2_L_t = result_L_t[:, 3]
        efficiency_t = (np.maximum(transmission_L_t, transmission_R_t) *
                        np.abs((transmission_L_t - transmission_R_t) /
                               (transmission_L_t + transmission_R_t)))

    # Write to file
    p_for_plot = dt.cir_parameters(cirs, L_js, {'L': 3873.0, 'R': 0.0},
                                   a_ops, i_op)
    data_dict.update({'a_inc_list': a_inc_list, 'L_js': [L_js[0], L_js[1]],
                      'gammas': p_for_plot["gammas"],
                      'delta': p_for_plot["delta"]})
    dt.write_h5(PATH_DATA+"power_dependence_Nexc" + "_delta_" +
                str(np.round(p_for_plot['delta'],
                             decimals=dt.determine_num_decimals_file_name(
                                 p_for_plot['delta']
                                 )
                             )
                    ) + "_v1.1.h5",
                data_dict)
    if TIME:
        # TODO
        dict_power_t = {'a_inc_list': a_inc_list_t, 'L_js': [L_js[0], L_js[1]],
                        'gammas': p_for_plot["gammas"], 'delta': p_for_plot["delta"],
                        'transmission_R': transmission_R_t,
                        'transmission_L': transmission_L_t,
                        'efficiency': efficiency_t,
                        'population_B_L': pop_B_L_t, 'population_D_L': pop_D_L_t,
                        'population_B_R': pop_B_R_t, 'population_D_R': pop_D_R_t,
                        'rho2_L': rho2_L_t, 'rho2_R': rho2_R_t}
        dt.write_h5(PATH_DATA+"power_dependence_t_Nexc" +
                    "_delta_" +
                    str(np.round(p_for_plot['delta'],
                                 decimals=dt.determine_num_decimals_file_name(
                                     p_for_plot['delta']
                                     )
                                 )
                        ) + "_v1.1.h5",
                    dict_power_t)

    if TEST:
        fig_1, ax = plt.subplots(2, 1)
        p_for_plot = dt.cir_parameters(cirs, L_js, {'L': 3873.0, 'R': 0.0},
                                       a_ops, i_op)

        ax[0].plot(a_inc_list**2 / np.mean(p_for_plot["gammas"]),
                   transmission_R, color='r', label='R')
        ax[0].plot(a_inc_list**2 / np.mean(p_for_plot["gammas"]),
                   transmission_L, color='b', label='L')
        ax[0].plot(a_inc_list**2 / np.mean(p_for_plot["gammas"]),
                   efficiency, color='black', label='M')
        ax[0].plot(a_inc_list**2 / np.mean(p_for_plot["gammas"]),
                   [0.67 for j in a_inc_list], '--', color='black',
                   label='M: 2 level best')
        if TIME:
            ax[0].plot(a_inc_list_t**2 / np.mean(p_for_plot["gammas"]),
                       transmission_R_t, 'o', color='r', label='R time')
            ax[0].plot(a_inc_list_t**2 / np.mean(p_for_plot["gammas"]),
                       transmission_L_t, 'o', color='b', label='L time')
            ax[0].plot(a_inc_list_t**2 / np.mean(p_for_plot["gammas"]),
                       efficiency_t, '+', color='black', label='M time')
        ax[0].set_xscale('log')
        ax[0].set_xlabel(r'$|a_\mathrm{inc}|^2 / \bar{\gamma}$')
        ax[0].set_ylabel('Transmission')
        ax[0].legend()

        ax[1].plot(a_inc_list**2 / np.mean(p_for_plot["gammas"]), pop_D_R,
                   color='r', label='|D>: R')
        ax[1].plot(a_inc_list**2 / np.mean(p_for_plot["gammas"]), pop_D_L,
                   color='b', label='|D>: L')
        ax[1].plot(a_inc_list**2 / np.mean(p_for_plot["gammas"]), pop_B_R,
                   '--', color='r', label='|B>: R')
        ax[1].plot(a_inc_list**2 / np.mean(p_for_plot["gammas"]), pop_B_L,
                   '--', color='b', label='|B>: L')
        ax[1].plot(a_inc_list**2 / np.mean(p_for_plot["gammas"]), rho2_L,
                   color='black', label=r'$\mathrm{Tr}(\rho^2)$: L')
        ax[1].plot(a_inc_list**2 / np.mean(p_for_plot["gammas"]), rho2_R,
                   color='grey', label=r'$\mathrm{Tr}(\rho^2)$: R')
        if TIME:
            ax[1].plot(a_inc_list_t**2 / np.mean(p_for_plot["gammas"]),
                       pop_D_R_t, 'o',
                       color='r', label='|D>: R time')
            ax[1].plot(a_inc_list_t**2 / np.mean(p_for_plot["gammas"]),
                       pop_D_L_t, 'o',
                       color='b', label='|D>: L time')
            ax[1].plot(a_inc_list_t**2 / np.mean(p_for_plot["gammas"]),
                       pop_B_R_t, 'x',
                       color='r', label='|B>: R time')
            ax[1].plot(a_inc_list_t**2 / np.mean(p_for_plot["gammas"]),
                       pop_B_L_t, 'x',
                       color='b', label='|B>: L time')
            ax[1].plot(a_inc_list_t**2 / np.mean(p_for_plot["gammas"]),
                       rho2_L_t, '+',
                       color='black', label=r'$\mathrm{Tr}(\rho^2)$: L time')
            ax[1].plot(a_inc_list_t**2 / np.mean(p_for_plot["gammas"]),
                       rho2_R_t, '+',
                       color='grey', label=r'$\mathrm{Tr}(\rho^2)$: R time')

        ax[1].set_xscale('log')
        ax[1].set_xlabel(r'$|a_\mathrm{inc}|^2 / \bar{\gamma}$')
        ax[1].set_ylabel('Population')
        ax[1].legend()

    if TEST:
        plt.show()

except ValueError as e:
    print(f"Error: {e}")
    exit(1)
