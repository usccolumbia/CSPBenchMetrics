import pandas as pd
import glob
import os
import sys
import argparse

parser = argparse.ArgumentParser(
        description=("groundtruth and predicted formulas name, they should be same")
    )

parser.add_argument(
    "--formula",
    type=str,
    default="Ca4S4",
    help="Formula name",
)

args = parser.parse_args(sys.argv[1:])

if not os.path.exists("../results"):
    os.mkdir("../results")

if os.path.exists(f"../results/distance_table_{args.formula}.csv"):
    os.remove(f"../results/distance_table_{args.formula}.csv")

groundtruth_cif = f'../data/ground_truth_structures/{args.formula}.cif'
predicted_cif = f'../data/predicted_structures/{args.formula}_*.cif'

files = glob.glob(predicted_cif)

for p in files:
    print('------------------------------------------')
    print('It is computing structure: ' + p)
    print('------------------------------------------')
    op1 = f'python3 distance_mul.py --cif {groundtruth_cif} --predicted {p} --formula {args.formula}'
    os.system(op1)

    
