import pandas as pd
import glob
import os
import sys
import argparse

parser = argparse.ArgumentParser(
    description=(
        "This script compares ground truth and predicted structure files. The files should "
        "have identical formulas in their names. For example, for a ground truth file named "
        "'SrTiO3.cif', the corresponding predicted file could be named 'TCSP_SrTiO3.cif'."
    )
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
    print('For structure: ' + p)
    print('------------------------------------------')
    op1 = f'python3 distance_mul.py --ground_truth {groundtruth_cif} --predicted {p} --formula {args.formula}'
    os.system(op1)

    
