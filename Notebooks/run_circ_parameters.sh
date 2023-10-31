#!/bin/bash

# Activate the Python virtual environment
source ~/py_envs/qiskit/bin/activate

# Pass the command line arguments to the Python script
python circ_parameters.py --L_js "$1" "$2"
