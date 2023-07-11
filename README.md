# CSPMetrics: Crystal Structure Prediction Performance Metrics

This repository contains the code and datasets for the paper [PDF]():

```
@article{PhysRevLett.120.145301,
  title = {Towards Quantitative Evaluation of Crystal Structure Prediction Performance},
  author = {Lai Wei, Qin Li, Sadman Omee, and Jianjun Hu},
  journal = {Arxiv},
  volume = {120},
  issue = {14},
  pages = {145301},
  numpages = {6},
  year = {2023},
  month = {Apr},
  publisher = {Arxiv},
  doi = {10.1103/PhysRevLett.120.145301},
  url = {https://link.aps.org/doi/10.1103/PhysRevLett.120.145301}
}
```


by <a href="http://mleg.cse.sc.edu" target="_blank">Machine Learning and Evolution Laboratory</a>, University of South Carolina.

## How to install the CSPMetrics package

Follow the instructions in the requirement.txt file to install the corresponding packages.

Please notice that you need to install a version of pytorch that corresponds to cuda.

## CSPMetrics
This is a package for benchmarking the performance of crystal structure prediction algorithms by calculating the quality of predicted structures compared to the ground truth structures using several common distance metrics, including m3gnet energy distance, minimal Wyckoff RMSE distance, minimal Wyckoff MAE distance, RMS distance, RMS anonymous distance, the Sinkhorn distance, the Chamfer distance, the Hausdorff distance, superpose RMSD distance, edit graph distance, XRD spectrum distance and fingerprint distance to standardize the comparison of material structure prediction models.  

## How to use the CSPMetrics for calculating the structure distance matrics

### Calculate the distance matrics between two structures

Go to the CSPBenchMetrics/code/ folder and then run the command below. 

The example here is to compute the distance metrics between the ground truth structure SrTiO3_gt.cif and the predicted structure SrTiO3_pred.cif. 
```
cd CSPBenchMetrics/code
python3 distance_single.py --cif ../data/SrTiO3_gt.cif --predicted ../data/SrTiO3_pred.cif
```
Then you can check the result in "~/CSPBenchMetrics/results/distance_table.csv".

If you want to test your own cases, you can upload the structure files and change the corresponding file names.

### Calculate the distance metrics between a ground truth structure and a folder of predicted structures

Here the example is to compute the distance metrics between the ground truth structure of Ca4S4 and the predicted structures of Ca4S4. 

First put the ground truth structure Ca4S4.cif into the data/ground_truth_structures folder with the file format as Ca4S4.cif

Next, put all the predicted structures into the data/predicted_structures/Ca4S4_0.cif, Ca4S4_1.cif, ...

Now, go to the CSPBenchMetrics/code/ folder and then run the command below. 

```
cd CSPBenchMetrics/code
python3 calculate_multiple_cif.py --formula Ca4S4
```
Then you can check the result in "~/CSPBenchMetrics/results/distance_table_formula_name.csv". The formula_name here is Ca4S4.

If you want to test your own cases, you can upload the structure files and change the corresponding name.
