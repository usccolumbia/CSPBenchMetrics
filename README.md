# Crystal Structure Prediction Performance Metrics
This repository contains the code and datasets for the paper:
**
**link



by <a href="http://mleg.cse.sc.edu" target="_blank">Machine Learning and Evolution Laboratory</a>, University of South Carolina.

### How to install the CSP metric package
Follow the requirement.txt file to install the corresponding packages.


### CSPBenchMetrics
Benchmark metrics for crystal structure prediction


### How to use the CSPMetrics 


### Calculate the distance matrics between two structures

An example is to compute distance metrics between the ground truth formula and the predicted formula.
```
python3 distance_single.py --cif {../data/SrTiO3_gt} --predicted {../data/SrTiO3_pred}
```

## Calculate the distance metrics between a ground truth structure and a folder of predicted structures

An example is to compute distance metrics between the ground truth formula and the predicted formula.
```
python3 calculate_multiple_cif.py --formula {Ca4S4}
```
