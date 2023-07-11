# CSPMetrics: Crystal Structure Prediction Performance Metrics

This repository contains the code and datasets for the paper:
**
**link


by <a href="http://mleg.cse.sc.edu" target="_blank">Machine Learning and Evolution Laboratory</a>, University of South Carolina.

## How to install the CSPMetrics package

Follow the instructions in the requirement.txt file to install the corresponding packages.

Please notice that you need to install a version of pytorch that corresponds to cuda.

## CSPBenchMetrics
This is a new benchmark platform by using several common distance metrics, including m3gnet energy distance, minimal Wyckoff RMSE distance, minimal Wyckoff MAE distance, RMS distance, RMS anonymous distance, the Sinkhorn distance, the Chamfer distance, the Hausdorff distance, superpose RMSD distance, edit graph distance, XRD spectrum distance and fingerprint distance to standardize the comparison of material structure prediction models.  

## How to use the CSPMetrics for calculating the distance matrics

### Calculate the distance matrics between two structures

Go to the CSPBenchMetrics/code/ folder and then copy the command below. 

An example is to compute distance metrics between the ground truth formula and the predicted formula of SrTiO3. 
```
python3 distance_single.py --cif ../data/SrTiO3_gt.cif --predicted ../data/SrTiO3_pred.cif
```
Then you can check the result in "~/CSPBenchMetrics/results/distance_table.csv".

If you want to test your cases, you can upload the structure files and change the corresponding name.

### Calculate the distance metrics between a ground truth structure and a folder of predicted structures

Go to the CSPBenchMetrics/code/ folder and then copy the command below. 

An example is to compute distance metrics between the ground truth formula and the predicted formula of Ca4S4. 
```
python3 calculate_multiple_cif.py --formula Ca4S4
```
Then you can check the result in "~/CSPBenchMetrics/results/distance_table_formula_name.csv". The formula_name here is Ca4S4.

If you want to test your cases, you can upload the structure files and change the corresponding name.
