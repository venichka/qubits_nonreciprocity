#!/bin/bash

# Activate the Python virtual environment
source ~/py_envs/qiskit/bin/activate

python data_gen_power.py --L_js 3.301e-9 3.3e-9 --N_exc 2 --t_max 1e1 --tau_max 1e1 --nmax 500 --nmax_time 50
python data_gen_power.py --L_js 3.3001e-9 3.3e-9 --N_exc 2 --t_max 1e1 --tau_max 1e1 --nmax 500 --nmax_time 50

python data_gen_power.py --L_js 3.301e-9 3.3e-9 --N_exc 3 --t_max 1e1 --tau_max 1e1 --nmax 500 --nmax_time 50
python data_gen_power.py --L_js 3.3001e-9 3.3e-9 --N_exc 3 --t_max 1e1 --tau_max 1e1 --nmax 500 --nmax_time 50

python data_gen_power.py --L_js 3.301e-9 3.3e-9 --N_exc 5 --t_max 1e1 --tau_max 1e1 --nmax 500 --nmax_time 50
python data_gen_power.py --L_js 3.3001e-9 3.3e-9 --N_exc 5 --t_max 1e1 --tau_max 1e1 --nmax 500 --nmax_time 50
