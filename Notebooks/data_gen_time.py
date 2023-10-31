import argparse
import ast
import qucat as qc
import qutip as qt
import scipy as sc
import matplotlib.pyplot as plt
import numpy as np
import h5py
import os


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


def path():
    home = os.path.expanduser("~")
    if home == "C:\\Users\\nnefedkin":
        PATH_FIGS = "D:/nnefedkin/Google_Drive/Work/In process/Projects/Qubits_nonreciprocity/Figs/"
        PATH_DATA = "D:/nnefedkin/Google_Drive/Work/In process/Projects/Qubits_nonreciprocity/Data/"
    elif home == "/home/nikita":
        PATH_FIGS = "/home/nikita/Documents/Work/Projects/Qubits_nonreciprocity/Figs/"
        PATH_DATA = "/home/nikita/Documents/Work/Projects/Qubits_nonreciprocity/Data/"
    elif home == "/Users/jimi":
        PATH_DATA = '/Users/jimi/Google Drive/Work/In process/Projects/Qubits_nonreciprocity/Data/'
        PATH_FIGS = '/Users/jimi/Google Drive/Work/In process/Projects/Qubits_nonreciprocity/Figs/'
    else:
        raise ValueError("Unknown home directory: {}".format(home))
    return PATH_DATA, PATH_FIGS


def determine_num_decimals_file_name(a0):
    abs_a0 = abs(a0)  # Get the absolute value of a0
    k = 0
    while abs_a0 < 1:
        abs_a0 *= 10
        k += 1
    return k + 1


def validate_input(L_js, N_exc, direction, a_inc, t_max, tau_max, nmax):
    if not isinstance(L_js, list) or not all(isinstance(x, float) for x in L_js):
        raise ValueError("L_js should be a list of floats")
    if not isinstance(N_exc, int):
        raise ValueError("N_exc should be an integer")
    if direction not in ['L', 'R']:
        raise ValueError("direction should be 'L' or 'R'")
    if not isinstance(a_inc, float):
        raise ValueError("a_inc should be a float (around 0.1*np.sqrt(1.5e+9))")
    if not isinstance(t_max, float):
        raise ValueError("t_max of t_D should be a float (usually 1e1)")
    if not isinstance(tau_max, float):
        raise ValueError("tau_max of t_D should be a float (usually 1e1 or 5e1)")
    if not isinstance(nmax, int):
        raise ValueError("nmax should be an integer (default = 100)")


def get_input_args():
    parser = argparse.ArgumentParser(description="Process input from a file or command line")
    parser.add_argument("--file", type=str, help="Input filename")
    parser.add_argument("--L_js", type=float, nargs='+', help="List of floats")
    parser.add_argument("--N_exc", type=int, help="Integer: number of excitations")
    parser.add_argument("--direction", type=str, choices=['L', 'R'], help="Direction ('L' or 'R')")
    parser.add_argument("--a_inc", type=float, default=3873.0, help="Input amplitude: float")
    parser.add_argument("--t_max", type=float, help="Integration time (number of t_D): float")
    parser.add_argument("--tau_max", type=float, help="Integration time for g1 and g2 (number of t_D): float")
    parser.add_argument("--test", action='store_true', help="Show test figs: True or False")
    parser.add_argument("--nmax", type=int, default=100, help="Integer: max number of L_j values")
    return parser.parse_args()


def load_config_from_file(file_name):
    try:
        with open(file_name, 'r') as file:
            data = file.read()
            config = ast.literal_eval(data)  # Safely evaluate the content as a Python literal

            if not isinstance(config, dict):
                raise ValueError("File should contain a dictionary")

            return config.get("L_js"), config.get("N_exc"), config.get("direction"), config.get("a_inc"), config.get("t_max"), config.get("tau_max"), config.get("test"), config.get("nmax")
    except FileNotFoundError:
        raise FileNotFoundError(f"Error: File '{file_name}' not found.")
    except ValueError as e:
        raise ValueError(f"Error: {e}")


def cir_parameters(cirs, L_js, a_inc, a_ops, i_op):
    gammas = 2.0*np.pi*np.array([cirs[i].loss_rates(Lj = L_js[i])[0] for i in range(len(cirs))])
    w_m = 2.0*np.pi*np.array([cirs[i].eigenfrequencies(Lj = L_js[i])[0] for i in range(len(cirs))])
    w_ext = w_m[0]
    delta = (w_m[0] - w_m[1]) / np.mean(gammas)
    t_j = np.array([0.0, (np.pi - delta) / w_ext])
    t_ij = np.array([[t_j[i] - t_j[j] for j in range(len(cirs))] for i in range(len(cirs))])
    return {
        'a_inc': a_inc,
        'w_m': w_m,
        'w_ext': w_ext,
        'A_m': 2.0*np.pi*np.array([cirs[i].anharmonicities(Lj = L_js[i])[0] for i in range(len(cirs))]),
        'gammas': gammas,
        'gamma_nr': 0*190.0e3 * 2.0*np.pi,
        'gamma_phi': 0*200.0e3 * 2.0*np.pi,
        'delta': delta,
        't_j': t_j,
        't_ij': t_ij,
        'Omega_R': np.array([1.0j*np.sqrt(0.5*gammas[j])*
                      (a_inc['R'] * np.exp(-1j * w_ext * t_j[j]) + 
                       a_inc['L'] * np.exp(1j * w_ext * t_j[j])) 
                      for j in range(len(cirs))]),
        'Omega': np.array([[0.5*np.sqrt(gammas[i]*gammas[j])*
                       np.sin(w_ext * np.abs(t_ij[i,j]))
                       for j in range(len(cirs))] 
                      for i in range(len(cirs))]),
        'Gamma': np.array([[np.sqrt(gammas[i]*gammas[j])*
                       np.cos(w_ext * np.abs(t_ij[i,j]))
                       for j in range(len(cirs))] 
                      for i in range(len(cirs))]),
        'Omega_R_nonRWA': np.array([1.0j*np.sqrt(0.5*gammas[j])*np.sqrt(w_ext/w_m[j])*
                      (a_inc['R'] * np.exp(-1j * w_ext * t_j[j]) + 
                       a_inc['L'] * np.exp(1j * w_ext * t_j[j])) 
                      for j in range(len(cirs))]),
        'Omega_nonRWA': np.array([[-1j*0.25*np.sqrt(gammas[i]*gammas[j])*
                       (
                           np.sqrt(w_m[i]/w_m[j])*np.exp(1j*w_m[i]*np.abs(t_ij[j,i])) - 
                           np.sqrt(w_m[j]/w_m[i])*np.exp(-1j*w_m[j]*np.abs(t_ij[i,j]))                              
                       )
                       for j in range(len(cirs))] 
                      for i in range(len(cirs))]),
        'Gamma_nonRWA': np.array([[0.5*np.sqrt(gammas[i]*gammas[j])*
                       (
                           np.sqrt(w_m[i]/w_m[j])*np.exp(1j*w_m[i]*np.abs(t_ij[j,i])) + 
                           np.sqrt(w_m[j]/w_m[i])*np.exp(-1j*w_m[j]*np.abs(t_ij[i,j]))                              
                       )
                       for j in range(len(cirs))] 
                      for i in range(len(cirs))]),
        # output operators
        'out_ops': {'R': a_inc['R']*i_op + sum([np.exp(-1j*w_m[j]*t_j[j])*
                                             np.sqrt(0.5*gammas[j])*a_ops[j] for j in range(len(cirs))]),
                    'L': a_inc['L']*i_op + sum([np.exp(1j*w_m[j]*t_j[j])*
                                             np.sqrt(0.5*gammas[j])*a_ops[j] for j in range(len(cirs))])},
        # for non-rotating frame, a_in * exp(-i w_ext t) as an input signal
        'out_ops_t': {'R': sum([np.exp(-1j*w_m[j]*t_j[j])*
                                             np.sqrt(0.5*gammas[j])*a_ops[j] for j in range(len(cirs))]),
                      'L': sum([np.exp(1j*w_m[j]*t_j[j])*
                                             np.sqrt(0.5*gammas[j])*a_ops[j] for j in range(len(cirs))])}
    }


def operators(N_excitations, N_qubits):
    # transition operators
    a_ops = [qt.destroy(N_excitations) if i == 0 else qt.qeye(N_excitations)
             for i in range(N_qubits)]

    for i in range(N_qubits):
        for j in range(1, N_qubits):
            if i == j:
                j_op = qt.destroy(N_excitations)
            else:
                j_op = qt.qeye(N_excitations)
            a_ops[i] = qt.tensor(a_ops[i], j_op)

    # unity operator
    i_op = qt.qeye(N_excitations)
    for j in range(1, N_qubits):
        i_op = qt.tensor(i_op, qt.qeye(N_excitations))

    return a_ops, i_op


def initial_state(N_qubits, N_excitations):
    psi0 = qt.basis(N_excitations, 0)  # start without excitations
    for j in range(N_qubits - 1):
        psi0 = qt.tensor(psi0, qt.basis(N_excitations, 0))
    rho0 = qt.ket2dm(psi0)
    return rho0


# Hamiltonians, jump operators, states
def hamiltonian_transmon(a_ops, parameters):
    p = parameters
    N_qubits = len(p["w_m"])
    H_0 = sum([(p['w_m'][j] - p['w_ext'])*a_ops[j].dag()*a_ops[j]
               - 0.5*p['A_m'][j]*a_ops[j].dag()*a_ops[j].dag()*a_ops[j]*a_ops[j]#(a_ops[j].dag() + a_ops[j])**4
               for j in range(N_qubits)])
    H_f = sum([p['Omega_R'][j]*a_ops[j] + np.conj(p['Omega_R'][j])*a_ops[j].dag() for j in range(N_qubits)])
    H_d = sum(sum([[p['Omega'][i,j]*a_ops[i].dag()*a_ops[j]#(a_ops[i] - a_ops[i].dag())*(a_ops[j] - a_ops[j].dag())
                    for j in range(N_qubits)] for i in range(N_qubits)], []))
    return H_0 + H_f + H_d, H_0, H_f, H_d


def hamiltonian_transmon_time_drive(a_ops, args):
    p = args
    N_qubits = len(p["w_m"])
    H_0 = sum([(p['w_m'][j] - p['w_ext'])*a_ops[j].dag()*a_ops[j] 
               - 0.5*p['A_m'][j]*a_ops[j].dag()*a_ops[j].dag()*a_ops[j]*a_ops[j]#(a_ops[j].dag() + a_ops[j])**4
               for j in range(N_qubits)])
    H_d = sum(sum([[p['Omega_nonRWA'][i,j]*a_ops[i].dag()*a_ops[j]#(a_ops[i] - a_ops[i].dag())*(a_ops[j] - a_ops[j].dag())
                    for j in range(N_qubits)] for i in range(N_qubits)], []))
    time_dep_coeff_m = [
        f'-2.0*sqrt(0.5*g[{j}])*sqrt(w_ext/w[{j}]) * \
        (a_L*sin(w_ext*(t+t_j[{j}])) + a_R*sin(w_ext*(t-t_j[{j}]))) * \
        exp(-1j*w_ext*t)'
        for j in range(N_qubits)
        ]
    time_dep_coeff_p = [
        f'-2.0*sqrt(0.5*g[{j}])*sqrt(w_ext/w[{j}]) * \
        (a_L*sin(w_ext*(t+t_j[{j}])) + a_R*sin(w_ext*(t-t_j[{j}]))) * \
        exp(1j*w_ext*t)'
        for j in range(N_qubits)
        ]
    H = [H_0, H_d,
        *[ [a_ops[i], time_dep_coeff_m[i]] for i in range(N_qubits) ],
        *[ [a_ops[i].dag(), time_dep_coeff_p[i]] for i in range(N_qubits) ]
        ]
    arguments_dict = {
        f'g': p['gammas'],
        'w_ext': p['w_ext'],
        f'w': p['w_m'],
        'a_L': p['a_inc']['L'],
        f't_j': p['t_j'],
        'a_R': p['a_inc']['R']
    }
    return H, arguments_dict


def jump_operators(a_ops, args):
    J_ops = []
    Gamma = args["Gamma_nonRWA"]
    N_qubits = len(args["w_m"])
    lam, b = sc.linalg.eig(Gamma)
    J_ops = [ np.sqrt(lam[i])*sum([b[j, i] * a_ops[j] for j in range(N_qubits) ]) for i in range(N_qubits)]
    return J_ops


def dark_bright_states(a_ops, args):
    J_ops = []
    Gamma = args["Gamma"]
    N_qubits = len(args["w_m"])
    N_excitations = a_ops[0].dims[0][0]  # size of the 1 subspace dimension
    lam, b = sc.linalg.eig(Gamma)
    J_ops = [ sum([b[j, i] * a_ops[j] for j in range(N_qubits) ]) for i in range(N_qubits)]
    # dark and bright operators
    for j in range(len(lam)):
        if lam[j] == np.min(lam):
            J_D = J_ops[j]
        elif lam[j] == np.max(lam):
            J_B = J_ops[j]
    # ground state
    psi0 = qt.basis(N_excitations, 0) # start without excitations
    for j in range(N_qubits - 1):
        psi0 = qt.tensor(psi0, qt.basis(N_excitations, 0))
    return qt.ket2dm(J_D.dag()*psi0), qt.ket2dm(J_B.dag()*psi0)


def dark_bright_states_nonRWA(a_ops, args):
    J_ops = []
    Gamma = args["Gamma_nonRWA"]
    N_qubits = len(args["w_m"])
    N_excitations = a_ops[0].dims[0][0]  # size of the 1 subspace dimension
    lam, b = sc.linalg.eig(Gamma)
    J_ops = [ sum([b[j, i] * a_ops[j] for j in range(N_qubits) ]) for i in range(N_qubits)]
    # dark and bright operators
    for j in range(len(lam)):
        if lam[j] == np.min(lam):
            J_D = J_ops[j]
        elif lam[j] == np.max(lam):
            J_B = J_ops[j]
    # ground state
    psi0 = qt.basis(N_excitations, 0) # start without excitations
    for j in range(N_qubits - 1):
        psi0 = qt.tensor(psi0, qt.basis(N_excitations, 0))
    return qt.ket2dm(J_D.dag()*psi0), qt.ket2dm(J_B.dag()*psi0)


# Hamiltonian eigenenergies
def hamiltonian_eigenenergies(a_in, L_j1, L_j2, cirs, a_ops, i_op,
                              direction: str = 'R', cutoff=None):
    if direction == 'R':
        a_inc = {'R': a_in, 'L': 0.0}
    elif direction == 'L':
        a_inc = {'R': 0.0, 'L': a_in}
    L_js = [L_j1, L_j2]
    p = cir_parameters(cirs, L_js, a_inc, a_ops, i_op)
    H, _, _, _ = hamiltonian_transmon(a_ops, p)
    return H.eigenenergies() if cutoff == None else H.eigenenergies()[cutoff[0]:cutoff[1]]


# Gamma eigenvalues
def gamma_eigvals(a_in, L_j1, L_j2, cirs, a_ops, i_op, direction: str = 'R'):
    if direction == 'R':
        a_inc = {'R': a_in, 'L': 0.0}
    elif direction == 'L':
        a_inc = {'R': 0.0, 'L': a_in}
    L_js = [L_j1, L_j2]
    p = cir_parameters(cirs, L_js, a_inc, a_ops, i_op)
    Gamma = p["Gamma"]
    return sc.linalg.eigvals(Gamma)


# Correlation functions for time dependent H
def g12_time_dep(H, rho0, taus, c_ops, a_op, args={}, func='g1'):
    # first calculate the occupation number as a function of time
    n = qt.mesolve(H, rho0, taus, c_ops, e_ops=[a_op.dag()*a_op], args=args,
                   progress_bar=True,
                   options=qt.Options(nsteps=1e9, rhs_reuse=False)).expect[0]

    # calculate the correlation function G1(2) and normalize with n(0)n(t) to
    # obtain g1(2)
    if func == 'g1':
        G1 = qt.correlation_2op_1t(H, rho0, taus, c_ops,
                                   a_op.dag(), a_op,
                                   solver='me',
                                   reverse=False, args=args,
                                   options=qt.Options(nsteps=1e9))
        return n, G1 / np.sqrt((n[0] * n))
    elif func == 'g2':
        G2 = qt.correlation_3op_1t(H, rho0, taus, c_ops,
                                   a_op.dag(), a_op.dag()*a_op, a_op,
                                   solver='me',
                                   args=args,
                                   options=qt.Options(nsteps=1e9))
        return n, G2 / (n[0] * n)


# Transmissions
def transmission(direction: str, rho, a_inc, out_ops):
    return qt.expect(out_ops[direction].dag()*out_ops[direction],
                     rho) / np.abs(a_inc)**2


def transmission_f(direction: str, rho, a_inc, out_ops):
    return np.abs(qt.expect(out_ops[direction], rho) / a_inc)


def main():
    args = get_input_args()

    if args.file:
        L_js, N_exc, DIRECTION, a_mag, t_max, tau_max, TEST, nmax = load_config_from_file(args.file)
    else:
        L_js = args.L_js
        N_exc = args.N_exc
        DIRECTION = args.direction
        a_mag = args.a_inc
        t_max = args.t_max
        tau_max = args.tau_max
        TEST = args.test
        NMAX = args.nmax

    try:
        PATH_DATA, PATH_FIGS = path()
        # Setting the system
        validate_input(L_js, N_exc, DIRECTION, a_mag, t_max, tau_max, NMAX)

        # a_mag = 0.1 * np.sqrt(1.5e+9)
        a_R = a_mag if DIRECTION == 'R' else 0e-2
        a_L = a_mag if DIRECTION == 'L' else 0e-2
        a_inc = {'R': a_R, 'L': a_L}

        cir_1 = qc.GUI('circuits/1_transmon.txt', edit=False, plot=False,
                       print_network=False)
        cir_2 = qc.GUI('circuits/1_transmon.txt', edit=False, plot=False,
                       print_network=False)
        cirs = [cir_1, cir_2]
        N_qubits = len(cirs)
        DIM = N_exc**N_qubits

        a_ops, i_op = operators(N_exc, N_qubits)

        p = cir_parameters(cirs, L_js, a_inc, a_ops, i_op)  # system parameters
        out_ops = p["out_ops"]

        # Computing H eigenvalues and g_B, g_D vs L_j
        L_j_list = np.linspace(2.80e-9, 3.80e-9, NMAX)
        H_eigs = np.array([hamiltonian_eigenenergies(a_mag, j, L_js[1],
                                                     cirs, a_ops, i_op)
                           for j in L_j_list])
        gammas_B = np.array([np.max(gamma_eigvals(a_mag, j, L_js[1],
                                                  cirs, a_ops, i_op))
                             for j in L_j_list])  # TODO: potential slow
        gammas_D = np.array([np.min(gamma_eigvals(a_mag, j, L_js[1],
                                                  cirs, a_ops, i_op))
                             for j in L_j_list])
        dict_Lj_dep = {'g_D(L_j)': np.real(gammas_D),
                       'g_B(L_j)': np.real(gammas_B),
                       'H_eig(L_j)': H_eigs, 'L_j_list': L_j_list}
        # TEST
        if TEST:
            fig_eg, ax = plt.subplots(1, 2)
            ax[0].plot(L_j_list/1e-9, np.real(gammas_B)/(2*np.pi*1e9), label=r'$\gamma_B$')
            ax[0].plot(L_j_list/1e-9, np.real(gammas_D)/(2*np.pi*1e9), label=r'$\gamma_D$')
            ax[0].plot([L_js[1]/1e-9, L_js[1]/1e-9], np.linspace(0.0, 0.17, 2), '--', color='black')
            ax[0].set_xlabel(r"$L_j$ (nH)")
            ax[0].set_ylabel(r"$\gamma_{B,D} / 2\pi$ (GHz)")
            ax[0].legend()

            ax[1].plot(L_j_list/1e-9, np.real(H_eigs)/(2*np.pi*1e9), 'o', ms=2)
            ax[1].plot(L_j_list/1e-9, [0.0 for j in L_j_list], '--', color='black')
            ax[1].plot(L_js[1]/1e-9, 0.0, 'o', color='r')
            ax[1].set_xlabel(r"$L_j$ (nH)")
            ax[1].set_ylabel(r"$\omega / 2\pi$ (GHz)")
            ax[1].legend()
            fig_eg.tight_layout()

        print('Computed: H eigenvalues \n')

        # Computing time evolution
        # ODE solver options
#        opts = qt.Options(atol=1e-16, rtol=1e-14, method='adams', order=12,
#                          nsteps=10000, first_step=0, max_step=0,min_step=0,
#                          average_expect=True, average_states=False, tidy=True,
#                          num_cpus=0, norm_tol=0.001, norm_t_tol=1e-06,
#                          norm_steps=5, rhs_reuse=False, rhs_filename=None,
#                          ntraj=500, gui=False, rhs_with_state=False,
#                          store_final_state=False, store_states=False,
#                          steady_state_average=False, seeds=None,
#                          normalize_output=True, use_openmp=None,
#                          openmp_threads=None)
        # Solver params
        w = np.real(sc.linalg.eigvals(p['Gamma']))
        tlist = np.linspace(0.0, t_max/np.min(w), 2000)
        rho0 = initial_state(N_qubits, N_exc)
        # Hamiltonian and Liiouvillian
        H, H_0, H_f, H_i = hamiltonian_transmon(a_ops, p)
        L_diss = sum([sum([p['Gamma'][i, j] * qt.lindblad_dissipator(a_ops[i],
                                                                     a_ops[j])
                           for j in range(N_qubits)]) for i in range(N_qubits)])
        c_ops_nr = [np.sqrt(p["gamma_nr"]) * a_ops[j] for j in range(len(a_ops))]
        c_ops_deph = [np.sqrt(p["gamma_phi"]) * (a_ops[j].dag()*a_ops[j] -
                                                 a_ops[j]*a_ops[j].dag())
                      for j in range(len(a_ops))]
        L_0 = qt.liouvillian(H, c_ops=c_ops_nr + c_ops_deph)
        L = L_0 + L_diss

        H_t, H_args = hamiltonian_transmon_time_drive(a_ops, p)
        J_ops = jump_operators(a_ops, p)

        # Solving H time independent
        rho_ss = qt.steadystate(L, method='direct', maxiter=1e9, tol=1e-24,
                                use_precond=True)
        result_t = qt.mesolve(L, rho0, tlist, c_ops=[], e_ops=[],
                              options=qt.Options(nsteps=1e9))
        # Solving time dependent
        result_tt = qt.mesolve(H_t, rho0, tlist,
                               c_ops=J_ops + c_ops_nr + c_ops_deph, e_ops=[],
                               args=H_args, progress_bar=True,
                               options=qt.Options(nsteps=1e9, rhs_reuse=False))
        rho_ss_t = result_tt.states[-1]

        # TEST
        if TEST:
            fig_0, ax = plt.subplots(3, 1)
            ax[0].plot(tlist, [result_t.states[i].diag()[1:DIM] for i in range(len(tlist))])
            ax[0].set_ylabel(r'$\rho_{jj}$')
            ax[1].plot(tlist, np.abs(qt.expect(a_ops[0], result_t.states)), color='r', label='qubit 1')
            ax[1].plot(tlist, np.abs(qt.expect(a_ops[1], result_t.states)), color='b', label='qubit 2')
            ax[1].set_ylabel('Dipole moment')
            ax[1].legend()
            ax[2].plot(tlist, np.abs(qt.expect(out_ops['R'], result_t.states)), color='r', label='R')
            ax[2].plot(tlist, np.abs(qt.expect(out_ops['L'], result_t.states)), color='b', label='L')
            ax[2].plot(tlist, [np.abs(a_inc[DIRECTION]) for i in tlist], '--', color='black', 
                       label='inc '+ DIRECTION)
            ax[2].set_xlabel('Time')
            ax[2].set_ylabel('Output')
            ax[2].legend()
            fig_0.suptitle(r"$a_R^2 = %s$, $a_L^2 = %s$" % (a_inc['R']**2/np.mean(p['gammas']), 
                                                        a_inc['L']**2/np.mean(p['gammas'])))
            fig_0.tight_layout()

        print('Computed: time evolution \n')

        # Coherence functions and spectra: time independent
        taulist = np.linspace(0.0, tau_max/np.min(w), 10000)
        g1_tau_R, _ = qt.coherence_function_g1(L, None, taulist, c_ops=[],
                                               a_op=out_ops['R'],
                                               options=qt.Options(nsteps=1e9))
        g1_tau_L, _ = qt.coherence_function_g1(L, None, taulist, c_ops=[],
                                               a_op=out_ops['L'],
                                               options=qt.Options(nsteps=1e9))
        g2_tau_R, _ = qt.coherence_function_g2(L, None, taulist, c_ops=[],
                                               a_op=out_ops['R'],
                                               options=qt.Options(nsteps=1e9))
        g2_tau_L, _ = qt.coherence_function_g2(L, None, taulist, c_ops=[],
                                               a_op=out_ops['L'],
                                               options=qt.Options(nsteps=1e9))
        # Spectrum
        _, spec_R = qt.spectrum_correlation_fft(taulist, g1_tau_R-g1_tau_R[-1])
        w_spec, spec_L = qt.spectrum_correlation_fft(taulist,
                                                     g1_tau_L - g1_tau_L[-1])
        print('Computed: time independent spectra \n')

        # Time dependent coherence functions and spectra
        _, g1_tau_R_t = g12_time_dep(H_t, result_tt.states[-1], taulist,
                                     J_ops + c_ops_nr + c_ops_deph,
                                     a_op=out_ops['R'], args=H_args)
        _, g1_tau_L_t = g12_time_dep(H_t, result_tt.states[-1], taulist,
                                     J_ops + c_ops_nr + c_ops_deph,
                                     a_op=out_ops['L'], args=H_args)
        _, g2_tau_R_t = g12_time_dep(H_t, result_tt.states[-1], taulist,
                                     J_ops + c_ops_nr + c_ops_deph,
                                     a_op=out_ops['R'], args=H_args, func='g2')
        _, g2_tau_L_t = g12_time_dep(H_t, result_tt.states[-1], taulist,
                                     J_ops + c_ops_nr + c_ops_deph,
                                     a_op=out_ops['L'], args=H_args, func='g2')
        # Spectrum
        _, spec_R_t = qt.spectrum_correlation_fft(taulist,
                                                  np.abs(g1_tau_R_t) - np.abs(np.mean(g1_tau_R_t[-(len(taulist)//10):-1])))
        w_spec_t, spec_L_t = qt.spectrum_correlation_fft(taulist,
                                                         np.abs(g1_tau_L_t) - np.abs(np.mean(g1_tau_L_t[-(len(taulist)//10):-1])))

        # TEST
        if TEST:
            fig_coh_t, ax = plt.subplots(2, 2)
            gs = ax[1, 0].get_gridspec()
# remove the underlying axes
            for ax_i in ax[1, :]:
                ax_i.remove()
            axbig = fig_coh_t.add_subplot(gs[1, :])
            ax[0, 0].plot(taulist, np.abs(g1_tau_L_t), color='b', label='L')
            ax[0, 0].plot(taulist, np.abs(g1_tau_L), '--',color='black', label='L')
            ax[0, 0].plot(taulist, np.abs(g1_tau_R_t), color='r', label='R')
            ax[0, 0].plot(taulist, np.abs(g1_tau_R), '-.', color='black', label='R')
            ax[0, 0].set_ylabel(r'$g^{(1)}(\tau)$')
            ax[0, 0].legend()
            ax[0, 1].plot(taulist, np.abs(g2_tau_L_t), color='b', label='L')
            ax[0, 1].plot(taulist, np.abs(g2_tau_R_t), color='r', label='R')
            ax[0, 1].plot(taulist, np.abs(g2_tau_L), '--', color='black', label='L')
            ax[0, 1].plot(taulist, np.abs(g2_tau_R), '-.', color='black', label='R')
            ax[0, 1].set_ylabel(r'$g^{(2)}(\tau)$')
            ax[0, 1].legend()
            axbig.plot(w_spec_t/(2*np.pi*1e6), spec_L_t, color='b', label='L')
            axbig.plot(w_spec/(2*np.pi*1e6), spec_L, '--', color='black', label='L')
            axbig.plot(w_spec/(2*np.pi*1e6), spec_R, '-.', color='black', label='R')
            axbig.plot(w_spec_t/(2*np.pi*1e6), spec_R_t, color='r', label='R')
            axbig.set_xlim(-5e1, 5e1)
            axbig.set_ylabel(r'$S$')
            axbig.set_xlabel(r'$\omega / 2\pi$ (MHz)')
            axbig.legend()
            fig_coh_t.suptitle(r"$a_R^2 = %s$, $a_L^2 = %s$" % (a_inc['R']**2/np.mean(p['gammas']),
                                                        a_inc['L']**2/np.mean(p['gammas'])))
            fig_coh_t.tight_layout()

        print('Computed: time dependent spectra \n')

        dict_time_dep = {'tlist': tlist, 'rho_ss': (rho_ss.data).toarray(),
                         'rho_ss_t': (rho_ss_t.data).toarray(),
                         'direction': DIRECTION, 'a_inc': a_mag,
                         'N_excitations': N_exc, 'delta': p["delta"],
                         'gammas': p['gammas'],
                         'transmission': transmission(DIRECTION, rho_ss, a_mag,
                                                      out_ops),
                         'transmission_f': transmission_f(DIRECTION, rho_ss,
                                                          a_mag, out_ops),
                         'transmission_t': transmission(DIRECTION, rho_ss_t,
                                                        a_mag, out_ops),
                         'transmission_f_t': transmission_f(DIRECTION, rho_ss_t,
                                                            a_mag, out_ops),
                         'dipole_1_t': qt.expect(a_ops[0],
                                                 result_tt.states),
                         'dipole_2_t': qt.expect(a_ops[1],
                                                 result_tt.states),
                         'output_R_t': qt.expect(p['out_ops']['R'],
                                                 result_tt.states),
                         'output_L_t': qt.expect(p['out_ops']['L'],
                                                 result_tt.states),
                         'dipole_1': qt.expect(a_ops[0],
                                               result_t.states),
                         'dipole_2': qt.expect(a_ops[1],
                                               result_t.states),
                         'output_R': qt.expect(p['out_ops']['R'],
                                               result_t.states),
                         'output_L': qt.expect(p['out_ops']['L'],
                                               result_t.states),
                         'g1_R': g1_tau_R, 'g1_L': g1_tau_L,
                         'g2_R': g2_tau_R, 'g2_L': g2_tau_L,
                         'g1_R_t': g1_tau_R_t, 'g1_L_t': g1_tau_L_t,
                         'g2_R_t': g2_tau_R_t, 'g2_L_t': g2_tau_L_t,
                         'taulist': taulist, 'w_spec': w_spec,
                         'spec_R': spec_R, 'spec_L': spec_L,
                         'spec_R_t': spec_R_t, 'spec_L_t': spec_L_t,
                         'w_spec_t': w_spec_t,
                         'L_js': L_js}
        dict_Lj_dep.update(dict_time_dep)

        # Write data to file
        write_h5(PATH_DATA+"time_dependence_Nexc" + str(N_exc) + "_" +
                 DIRECTION + "_ainc_" +
                 str(np.round(a_inc[DIRECTION]**2/np.mean(p['gammas']),
                              decimals=2))+"_delta_" +
                 str(np.round(p['delta'],
                              decimals=determine_num_decimals_file_name(
                                  p['delta']))) + "_v1.1.h5", dict_Lj_dep)

        if TEST:
            plt.show()

    except ValueError as e:
        print(f"Error: {e}")
        exit(1)


if __name__ == "__main__":
    main()
