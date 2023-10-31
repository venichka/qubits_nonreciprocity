#!/bin/bash

# Activate the Python virtual environment
source ~/py_envs/qiskit/bin/activate

# Define common parameters
L_range=("3.301e-9 3.3e-9" "3.3001e-9 3.3e-9")
N_exc=("3" "5")
direction=("L" "R")
a_inc="3873e0"
t_max="1e1"
tau_max="1e1"
nmax="500"
nmax_time="50"

# Total number of iterations
total_iterations=$(( ${#L_range[@]} * ${#N_exc[@]} * ${#direction[@]} ))

# Initialize a counter for progress
progress=0

# Loop through combinations of parameters
for L in "${L_range[@]}"; do
  for exc in "${N_exc[@]}"; do
    for dir in "${direction[@]}"; do

      # Increment progress counter
      ((progress++))

      # Print progress message
      echo "Progress: $progress out of $total_iterations"

      # Run data generation
      python data_gen_time.py --L_js $L --N_exc $exc --direction $dir --a_inc $a_inc --t_max $t_max --tau_max $tau_max --nmax $nmax

      # Run power dependence if applicable
      if [ "$dir" == "L" ]; then
        python data_gen_power.py --L_js $L --N_exc $exc --t_max $t_max --tau_max $tau_max --nmax $nmax
      fi
    done
  done
done

# Print completion message
echo "Script completed."
