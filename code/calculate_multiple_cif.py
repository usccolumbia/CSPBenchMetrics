import pandas as pd
import glob
import os
import sys

#df = pd.read_csv('formulas.csv')['full_formula']
# print(df)
#if os.path.exists('distance_test.csv'):
#    os.remove('distance_test.csv')

gt_cif = f"/mul_cif_folder/CSP-metrics/ground_truth_cif/BaNaBi.cif"
algo = 'pso'
formula = 'Ba3Na3Bi3'
path = f"/mul_cif_folder/predicted_cif/results_pso/structures/Ca4S4_*.cif"
files = glob.glob(path)
print(len(files))

for p in files:
    pred_path = p
    print(pred_path)
    op1 = f'python3 distance_v5_mul.py --cif {gt_cif} --predicted {pred_path} --formula {formula} --algo {algo}'
    os.system(op1)
    