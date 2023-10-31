@echo off
setlocal enabledelayedexpansion

echo Script started.

call python data_gen_power.py --L_js 3.301e-9 3.3e-9 --N_exc 2 --t_max 1e1 --tau_max 1e1 --nmax 500 --nmax_time 50
call python data_gen_power.py --L_js 3.3001e-9 3.3e-9 --N_exc 2 --t_max 1e1 --tau_max 1e1 --nmax 500 --nmax_time 50

call python data_gen_power.py --L_js 3.301e-9 3.3e-9 --N_exc 3 --t_max 1e1 --tau_max 1e1 --nmax 500 --nmax_time 50
call python data_gen_power.py --L_js 3.3001e-9 3.3e-9 --N_exc 3 --t_max 1e1 --tau_max 1e1 --nmax 500 --nmax_time 50

call python data_gen_power.py --L_js 3.301e-9 3.3e-9 --N_exc 5 --t_max 1e1 --tau_max 1e1 --nmax 500 --nmax_time 50
call python data_gen_power.py --L_js 3.3001e-9 3.3e-9 --N_exc 5 --t_max 1e1 --tau_max 1e1 --nmax 500 --nmax_time 50


call python data_gen_time.py --L_js 3.3001e-9 3.3e-9 --N_exc 2 --direction L --a_inc 3873e0 --t_max 1e1 --tau_max 1e1 --nmax 500
call python data_gen_time.py --L_js 3.3001e-9 3.3e-9 --N_exc 2 --direction R --a_inc 3873e0 --t_max 1e-2 --tau_max 1e1 --nmax 500

call python data_gen_time.py --L_js 3.301e-9 3.3e-9 --N_exc 5 --direction L --a_inc 3873e0 --t_max 1e1 --tau_max 1e1 --nmax 500
call python data_gen_time.py --L_js 3.3001e-9 3.3e-9 --N_exc 5 --direction L --a_inc 3873e0 --t_max 1e1 --tau_max 1e1 --nmax 500

call python data_gen_time.py --L_js 3.301e-9 3.3e-9 --N_exc 5 --direction R --a_inc 3873e0 --t_max 1e-2 --tau_max 1e1 --nmax 500
call python data_gen_time.py --L_js 3.3001e-9 3.3e-9 --N_exc 5 --direction R --a_inc 3873e0 --t_max 1e-2 --tau_max 1e1 --nmax 500

call python data_gen_time.py --L_js 3.301e-9 3.3e-9 --N_exc 3 --direction L --a_inc 3873e0 --t_max 1e1 --tau_max 1e1 --nmax 500
call python data_gen_time.py --L_js 3.3001e-9 3.3e-9 --N_exc 3 --direction L --a_inc 3873e0 --t_max 1e1 --tau_max 1e1 --nmax 500

call python data_gen_time.py --L_js 3.301e-9 3.3e-9 --N_exc 3 --direction R --a_inc 3873e0 --t_max 1e-2 --tau_max 1e1 --nmax 500
call python data_gen_time.py --L_js 3.3001e-9 3.3e-9 --N_exc 3 --direction R --a_inc 3873e0 --t_max 1e-2 --tau_max 1e1 --nmax 500


REM Print completion message
echo Script completed.

