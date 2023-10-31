import argparse
import ast

def validate_input(N_excitations, DIRECTION, L_js):
    if not isinstance(N_excitations, int):
        raise ValueError("N_excitations should be an integer")
    if DIRECTION not in ['L', 'R']:
        raise ValueError("DIRECTION should be 'L' or 'R'")
    if not isinstance(L_js, list) or not all(isinstance(x, float) for x in L_js):
        raise ValueError("L_js should be a list of floats")

parser = argparse.ArgumentParser(description="Process input from a file or command line")
parser.add_argument("--file", type=str, help="Input filename")
parser.add_argument("--N_excitations", type=int, help="Number of excitations")
parser.add_argument("--DIRECTION", type=str, choices=['L', 'R'], help="Direction ('L' or 'R')")
parser.add_argument("--L_js", type=float, nargs='+', help="List of floats")

args = parser.parse_args()

if args.file:
    try:
        with open(args.file, 'r') as file:
            data = file.read()
            config = ast.literal_eval(data)  # Safely evaluate the content as a Python literal

            if not isinstance(config, dict):
                raise ValueError("File should contain a dictionary")

            N_excitations = config.get("N_excitations")
            DIRECTION = config.get("DIRECTION")
            L_js = config.get("L_js")

            validate_input(N_excitations, DIRECTION, L_js)

    except FileNotFoundError:
        print(f"Error: File '{args.file}' not found.")
        exit(1)
    except ValueError as e:
        print(f"Error: {e}")
        exit(1)
else:
    try:
        validate_input(args.N_excitations, args.DIRECTION, args.L_js)
        N_excitations = args.N_excitations
        DIRECTION = args.DIRECTION
        L_js = args.L_js
    except ValueError as e:
        print(f"Error: {e}")
        exit(1)

# Now, you can use N_excitations, DIRECTION, and L_js in your script
print("N_excitations:", N_excitations)
print("DIRECTION:", DIRECTION)
print("L_js:", L_js)
