import point_cloud_utils as pcu
import numpy as np
from pymatgen.analysis.molecule_structure_comparator import CovalentRadius
from pymatgen.io.cif import CifParser
from sklearn.metrics import mean_squared_error, mean_absolute_error
import glob
import random
from scipy.stats import pearsonr 
import scipy.stats
import matplotlib.pyplot as plt
import numpy as np
from superpose3d import Superpose3D
from math import *
import sys
import argparse
import pandas as pd
import re
from sklearn.metrics import mean_squared_error, mean_absolute_error
from pymatgen.analysis.structure_matcher import StructureMatcher as SM
from cifparser import CifParserExpand
import glob
from m3gnet.models import M3GNet
from pymatgen.core import Structure
import torch
from scipy.stats import rankdata
from scipy import interpolate
from torch_geometric.utils import dense_to_sparse, degree, add_self_loops
from pymatgen.core.periodic_table import Element
from pymatgen.io.cif import CifParser
from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from cyged.graph_pkg_core.algorithm.graph_edit_distance import GED
from cyged.graph_pkg_core.edit_cost.edit_cost_vector import EditCostVector
from cyged.graph_pkg_core.graph.edge import Edge
from cyged.graph_pkg_core.graph.graph import Graph
from cyged.graph_pkg_core.graph.label.label_edge import LabelEdge
from cyged.graph_pkg_core.graph.label.label_node_vector import LabelNodeVector
from cyged.graph_pkg_core.graph.node import Node
import pickle
import time
import os
from matminer.featurizers.site import CrystalNNFingerprint  # matminer version = 0.6.2 (later version will give same calculation results).
from matminer.featurizers.structure import SiteStatsFingerprint
from scipy.spatial import distance_matrix
from pymatgen.analysis.molecule_structure_comparator import CovalentRadius
from cifparser import CifParserExpand
import itertools
from matminer.featurizers.structure.matrix import OrbitalFieldMatrix
import pymatgen.analysis.diffraction.xrd as xrd
from scipy.stats import gaussian_kde
m3gnet_e_form = M3GNet.load()

parser = argparse.ArgumentParser(
        description=(
            "rmsd,mae,rms and accuracy"
            "calculate rmsd,mae,rms and contact map accuracy"
        )
    )

parser.add_argument(
    "--cif",
    type=str,
    default="../data/SrTiO3_gt.cif",
    metavar="PATH",
    help="Path to target cif file",
)

parser.add_argument(
    "--predicted",
    default="../data/SrTiO3_pred.cif",
    metavar="PATH",
    help="Path to predicted cif file",
)
parser.add_argument(
    "--formula",
    type=str,
    default="Ca4S4",
    help="Formula name",
)

args = parser.parse_args(sys.argv[1:])

s1 = Structure.from_file(args.cif)
s2 = Structure.from_file(args.predicted)

ged = GED(EditCostVector(1., 1., 1., 1., 'euclidean'))

def threshold_sort(matrix, threshold, neighbors, reverse=False, adj=False):
    mask = matrix > threshold
    distance_matrix_trimmed = np.ma.array(matrix, mask=mask)
    if reverse == False:
        distance_matrix_trimmed = rankdata(
            distance_matrix_trimmed, method="ordinal", axis=1
        )
    elif reverse == True:
        distance_matrix_trimmed = rankdata(
            distance_matrix_trimmed * -1, method="ordinal", axis=1
        )
    distance_matrix_trimmed = np.nan_to_num(
        np.where(mask, np.nan, distance_matrix_trimmed)
    )
    distance_matrix_trimmed[distance_matrix_trimmed > neighbors + 1] = 0

    if adj == False:
        distance_matrix_trimmed = np.where(
            distance_matrix_trimmed == 0, distance_matrix_trimmed, matrix
        )
        return distance_matrix_trimmed
    elif adj == True:
        adj_list = np.zeros((matrix.shape[0], neighbors + 1))
        adj_attr = np.zeros((matrix.shape[0], neighbors + 1))
        for i in range(0, matrix.shape[0]):
            temp = np.where(distance_matrix_trimmed[i] != 0)[0]
            adj_list[i, :] = np.pad(
                temp,
                pad_width=(0, neighbors + 1 - len(temp)),
                mode="constant",
                constant_values=0,
            )
            adj_attr[i, :] = matrix[i, adj_list[i, :].astype(int)]
        distance_matrix_trimmed = np.where(
            distance_matrix_trimmed == 0, distance_matrix_trimmed, matrix
        )
        return distance_matrix_trimmed, adj_list, adj_attr

    ##Slightly edited version from pytorch geometric to create edge from gaussian basis
class GaussianSmearing(torch.nn.Module):
    def __init__(self, start=0.0, stop=5.0, resolution=50, width=0.05, **kwargs):
        super(GaussianSmearing, self).__init__()
        offset = torch.linspace(start, stop, resolution)
        # self.coeff = -0.5 / (offset[1] - offset[0]).item() ** 2
        self.coeff = -0.5 / ((stop - start) * width) ** 2
        self.register_buffer("offset", offset)

    def forward(self, dist):
        dist = dist.unsqueeze(-1) - self.offset.view(1, -1)
        return torch.exp(self.coeff * torch.pow(dist, 2))

def stripSymbol(elm):
  retstring = ""  
  for c in elm:
        if c.isalpha():
            retstring += c 
  return retstring

def read_cif_to_graph(struct):
    n = len(struct.sites)
    graph_source = Graph('gr_source', 'gr_source.gxl', n)
    for i, s in enumerate(struct.sites):
        #print(s.species_string)
        e=stripSymbol(s.species_string)
        E = Element(e)
        graph_source.add_node(Node(i, LabelNodeVector(np.array([float(E.Z)]))))
    m=struct.distance_matrix
    dis_cutoff = 2.5
    neig_cutoff = 12
    distance_matrix_trimmed = threshold_sort(
                m,
                dis_cutoff,
                neig_cutoff,
                adj=False,
            )

    distance_matrix_trimmed = torch.Tensor(distance_matrix_trimmed)
    out = dense_to_sparse(distance_matrix_trimmed)
    edge_index = out[0]
    edge_weight = out[1]
    for i, e in enumerate(edge_index[0]):
        row=edge_index[0][i] 
        column=edge_index[1][i]
        dis=m[row][column]
        graph_source.add_edge(Edge(edge_index[0][i], edge_index[1][i], LabelEdge(dis)))
        #break
    return graph_source


#graph edit distance
# https://github.com/CheshireCat12/graph-matching-core
def GraphEditDistance(s1,s2):
    g1 = read_cif_to_graph(s1)
    g2 = read_cif_to_graph(s2)
    edit_cost = ged.compute_edit_distance(g1, g2)
    return edit_cost

# Set a function for calculating the structure fingerprints.
ssf = SiteStatsFingerprint(CrystalNNFingerprint.from_preset('ops', distance_cutoffs=None, x_diff_weight=0),
                           stats=('mean', 'std_dev', 'minimum', 'maximum'))
def CrystalFingerDistance(s1,s2):
    f1 = np.array(ssf.featurize(s1))
    f2 = np.array(ssf.featurize(s2))
    dist = np.linalg.norm(f1 - f2)
    return dist

fingerPrint_dist = (CrystalFingerDistance(s1,s2))

#superpose point cloud distances with alignment by rotation and translation.
# https://github.com/jewettaij/superpose3d_cpp
def superpose3d_rmsd(s1, s2): # s1, and s2 are to pymatgen structures.
    gt = []
    pred = []
    for x in s1.sites:
        gt.append(x.coords)
    for y in s2.sites:
        pred.append(y.coords)

    a = np.array(gt)
    b = np.array(pred)
    X = a
    x = b

    result = Superpose3D(X,x, None, False, True)
    #print('weights: '+str([1.0,1.0,1.0,1.0]))
    #print(' (quaternion = '+str(result[1])+')\n')

    Xshifted = [ [X[i][0],X[i][1]+100, X[i][2]] for i in range(0,len(X))]
    xscaled  = [ [2*x[i][0],2*x[i][1],    2*x[i][2]] for i in range(0,len(x))]
    xscshift = [ [2*x[i][0],2*x[i][1]+200,2*x[i][2]] for i in range(0,len(x))]

    # Now try again using the translated, rescaled coordinates:

    # now test weights, rescale, and quaternions
    w = []
    w_rand = random.randint(1,3)
    for i in range(0, len(x)):
        w.append(w_rand)
    result = Superpose3D(X, xscshift, w, True)
    # Does the RMSD returned in result[0] match the RMSD calculated manually?
    R = np.array(result[1])              # rotation matrix
    T = np.array(result[2]).transpose()  # translation vector (3x1 matrix)
    c = result[3]                        # scalar
    if len(X) > 0:
        _x = np.array(xscshift).transpose()
        # _xprime = c*R*_x + T   <-- syntax is depreciated
        _xprime = c*np.matmul(R,_x) + np.outer(T, np.array([1]*len(X)))
        xprime = np.array(_xprime.transpose()) # convert to length 3 numpy array
    else:
        xprime = np.array([])

    RMSD = 0.0
    sum_w = 0.0
    for i in range(0, len(X)):
        RMSD += w[i]*((X[i][0] - xprime[i][0])**2 +
                      (X[i][1] - xprime[i][1])**2 +
                      (X[i][2] - xprime[i][2])**2)
        sum_w += w[i]

    if len(X) > 0:
        RMSD = sqrt(RMSD / sum_w)

    # assert(abs(RMSD - result[0]) < 1.0e-6)
    return RMSD



#energy distance
def en_dist(s1, s2):
    s_gt=s1
    s_pred=s2
    e_form_gt = m3gnet_e_form.predict_structure(s_gt)
    e_form_predict = m3gnet_e_form.predict_structure(s_pred)
    energy_distance = float(abs(e_form_gt - e_form_predict).numpy())
    return energy_distance

#calculate geometric distances..from point_cloud_3d package.
def point_cloud_distance(s1, s2):
    s_gt=s1
    s_pred=s2

    gt = []
    pred = []
    for x in s_gt.sites:
        gt.append(x.coords)
    for y in s_pred.sites:
        pred.append(y.coords)

    a = np.array(gt)
    b = np.array(pred)
    
    ##proximate Wasserstein (Sinkhorn) distance between two point clouds

    # M is a 100x100 array where each entry  (i, j) is the squared distance between point a[i, :] and b[j, :]
    M = pcu.pairwise_distances(a, b)

    # w_a and w_b are masses assigned to each point. In this case each point is weighted equally.
    w_a = np.ones(a.shape[0])
    w_b = np.ones(b.shape[0])

    # P is the transport matrix between a and b, eps is a regularization parameter, smaller epsilons lead to
    # better approximation of the true Wasserstein distance at the expense of slower convergence
    P = pcu.sinkhorn(w_a, w_b, M, eps=1e-3)

    # To get the distance as a number just compute the frobenius inner product <M, P>
    sinkhorn_dist = (M*P).sum()

    chamfer_dist = pcu.chamfer_distance(a, b)

    # Compute one-sided squared Hausdorff distances
    hausdorff_a_to_b = pcu.one_sided_hausdorff_distance(a, b)
    hausdorff_b_to_a = pcu.one_sided_hausdorff_distance(b, a)

    # Take a max of the one sided squared  distances to get the two sided Hausdorff distance
    hausdorff_dist = pcu.hausdorff_distance(a, b)

    # Find the index pairs of the two points with maximum shortest distancce
    hausdorff_b_to_a, idx_b, idx_a = pcu.one_sided_hausdorff_distance(b, a, return_index=True)
    assert np.abs(np.sum((a[idx_a] - b[idx_b])**2) - hausdorff_b_to_a**2) < 1e-5, "These values should be almost equal"

    # Find the index pairs of the two points with maximum shortest distancce
    hausdorff_dist, idx_b, idx_a = pcu.hausdorff_distance(b, a, return_index=True)
    assert np.abs(np.sum((a[idx_a] - b[idx_b])**2) - hausdorff_dist**2) < 1e-5, "These values should be almost equal"  
    return sinkhorn_dist, chamfer_dist, hausdorff_dist




def adj(structure):
    d2 = structure.distance_matrix
    for i, s1 in enumerate(structure.sites):
        for j, s2 in enumerate(structure.sites):
            if i == j:
                d2[i, j] = 0
            else:
                e1 = s1.as_dict()['species'][0]['element']
                e2 = s2.as_dict()['species'][0]['element']
                z='      '
                z1=z[:6-len(e2)]+e2
                z2=z[:6-len(e1)]+e1
                if e1 + z1 in thresholds:
                    max_len = float(thresholds[thresholds.rindex(e1 + z1) + len(e1 + z1) + 15:thresholds.rindex(
                        e1 + z1) + len(e1 + z1) + 22])
                    min_len = float(thresholds[thresholds.rindex(e1 + z1) + len(e1 + z1) + 4:thresholds.rindex(
                        e1 + z1) + len(e1 + z1) + 11])
                elif e2 + z2 in thresholds:
                    max_len = float(thresholds[
                                    thresholds.rindex(e2 + z2) + len(e2 + z2) + 15:thresholds.rindex(
                                        e2 + z2) + len(e2 + z2) + 22])
                    min_len = float(thresholds[
                                    thresholds.rindex(e2 + z2) + len(e2 + z2) + 4:thresholds.rindex(
                                        e2 + z2) + len(e2 + z2) + 11])
                else:
                    max_len = 0.00000
                    min_len = 0.00000
                if min_len <= d2[i, j] <= max_len:
                    d2[i, j] = 1
                else:
                    d2[i, j] = 0
    d2 = d2.astype(int)
    return d2


with open('style.ini', "r") as fileobj:
    thresholds = fileobj.read()


#distance matrix similarity score
def distance_matrix_score(gt_cif,predicted_cif):
    parser1 = CifParserExpand(gt_cif)
    structure1 = parser1.get_structures()[0]

    try:
        parser2 = CifParserExpand(predicted_cif)
        structure2 = parser2.get_structures()[0]
        d2 = adj(structure2)
    except ValueError:
        fitness = 9990
    else:
        d1 = adj(structure1)
        if d2.shape != d1.shape:
            fitness = 9999
        else:
            d2 = adj(structure2)
            n = np.sum(d2 * d1 == 1)
            n1 = np.sum(d1 == 1)
            n2 = np.sum(d2 == 1)

            try:
                fitness = 2 * n / (n1 + n2)
                if fitness <= 1:
                    fitness = 1 - fitness
            except ZeroDivisionError:
                fitness = 9995
    return fitness

def wyckoff_distance(targetcif,predictedcif):
    with open(targetcif, "r") as fileobj:
        text1 = fileobj.read()
        co = re.findall(r'\d+\.\d+', text1[text1.rindex('occupancy'):])
        cof1 = []
        for i in co:
            cof1.append(float(i))
        cof1 = np.array(cof1).reshape(int(len(cof1) / 3), 3)
        zs = text1[text1.rindex('_atom_site_occupancy') + len('_atom_site_occupancy') + 3:]
        mpl = zs.split('  ')
        mpl = np.array(mpl).reshape(int(len(mpl) / 7), 7)
        mulpl = []
        for u in range(len(mpl)):
            mulpl.append(mpl[u][2])

    with open(predictedcif, "r") as fileobj:
        textpr = fileobj.read()
        copr = re.findall(r'\d+\.\d+', textpr[textpr.rindex('occupancy'):])
        coprf1 = []
        for i in copr:
            coprf1.append(float(i))
        coprf1 = np.array(coprf1).reshape(int(len(coprf1) / 3), 3)


    parser1 = CifParserExpand(targetcif)
    structure1 = parser1.get_structures()[0]
 
    try:
        parser2 = CifParserExpand(predictedcif)
        structure2 = parser2.get_structures()[0]
        d2=adj(structure2)
    except ValueError:
        fitness = -9990
    else:
        d1 = adj(structure1)
        if d2.shape != d1.shape:
            fitness = -9999
        else:
            n = np.sum(d2 * d1 == 1)
            n1 = np.sum(d1 == 1)
            n2 = np.sum(d2 == 1)

            try:
                fitness = 2 * n / (n1 + n2)
            except ZeroDivisionError:
                fitness = -9995

    if len(cof1) == len(coprf1):
        cos2 = parser2.cos
        cos2t = []
        for i in range(len(cos2)):
            cos2t.append(cos2[-i - 1])
        tem = []
        for k in cos2t:
            tem += k

        te = np.array(tem).flatten()
        fz = []
        for w in mulpl:
            fz.append(te[:3 * int(w)])
            te = te[3 * int(w):]
        te = fz
        cos2f = []
        for i in range(len(te)):
            cos2f.append(te[i].reshape(int(len(te[i]) / 3), 3))

        co3 = []
        for i in range(len(cof1)):
            for j in range(len(cos2f[i])):
                co3.append(np.sqrt(mean_squared_error(cof1[i], cos2f[i][j])))

        tempe = []
        for w in mulpl:
            tempe.append(co3[:int(w)])
            co3 = co3[int(w):]
        co3 = tempe

        co4 = []
        for u in range(len(co3)):
            k = np.argmin(co3[u])
            co4.append(cos2f[u][k])

        wyckoff_rmse = np.sqrt(mean_squared_error(cof1.flatten(), np.array(co4).flatten()))
        wyckoff_mae = mean_absolute_error(cof1.flatten(), np.array(co4).flatten())
    else:
        wyckoff_rmse=5555
        wyckoff_mae=5555

    matrix_similarity = fitness
    return wyckoff_rmse,wyckoff_mae
    

wyckoff_rmse,wyckoff_mae = wyckoff_distance(args.cif,args.predicted)

#pymatgen RMS distances
def PMD(structure1, structure2):
    sm = SM(ltol=0.6, stol=0.6,angle_tol=20)
    try:
        rms = sm.get_rms_dist(structure1, structure2)[0]
    except TypeError:
        rms='None'
    try:
        rms_anonymous = sm.get_rms_anonymous(structure1, structure2)[0]
    except TypeError:
        rms_anonymous='None'
    return rms,rms_anonymous

rms,rms_anonymous = PMD(s1,s2)

#XRD and OFM features
ofm= OrbitalFieldMatrix()

def smooth(lt):
    X_data = list(itertools.chain.from_iterable([[x[0]] * int(x[1] * 10) for x in lt]))
    X_plot = np.linspace(0, 90, 901)
    kde = gaussian_kde(X_data, bw_method=.01)
    output = kde.evaluate(X_plot)
    return output

def convert_to_powder(cif_str, file):
    crystal_struct = Structure.from_file(file)
    powd = xrd.XRDCalculator()
    powd_patt = powd.get_pattern(crystal_struct)
    return [[x, y] for x, y in zip(powd_patt.x, powd_patt.intensity)]

def calc(gt_file, pred_file):
    s_gt = Structure.from_file(gt_file)
    s_pred = Structure.from_file(pred_file)
    gt_XRD_feature = smooth(convert_to_powder(open(gt_file, 'rb').read(), gt_file))
    pred_XRD_feature = smooth(convert_to_powder(open(pred_file, 'rb').read(), pred_file))
    gt_OFM_feature = ofm.featurize(s_gt)
    pred_OFM_feature = ofm.featurize(s_pred)
    # Compute Euclidean distance
    XRD_euclidean_distance = np.linalg.norm(gt_XRD_feature - pred_XRD_feature)
    OFM_euclidean_distance = np.linalg.norm(gt_OFM_feature - pred_OFM_feature)
    print("XRD Euclidean Distance:", XRD_euclidean_distance)
    print("OFM Euclidean Distance:", OFM_euclidean_distance)
    return XRD_euclidean_distance, OFM_euclidean_distance

XRD_dist, OFM_dist = calc(args.cif, args.predicted)
energy_dis = en_dist(s1, s2)

superpose_rmsd = superpose3d_rmsd(s1, s2)

a, b=  rms,rms_anonymous
x, y, z= point_cloud_distance(s1, s2)

edit_cost = GraphEditDistance(s1,s2)
result_dict = {'structure': args.predicted.split('/')[-1] , 'energy_distance':energy_dis, 
'wyckoff_rmse':wyckoff_rmse, 'wyckoff_mae':wyckoff_mae, 'rms_anonymous':rms_anonymous, 
'rms_dist': rms, 'sinkhorn_dist': x, 'chamfer_dist': y, 'hausdorff_dist': z, 
'superpose_rmsd': superpose_rmsd, 'edit_graph_distance':edit_cost, 'fingerPrint':fingerPrint_dist, 'XRD_dist': XRD_dist, 'OFM_dist': OFM_dist}
print(result_dict)
dct = {k:[v] for k,v in result_dict.items()}  # WORKAROUND
df = pd.DataFrame(dct)
df.drop_duplicates(inplace=True)

if not os.path.exists("../results"):
    os.mkdir("../results")
if not os.path.exists(f"../results/distance_table_{args.formula}.csv"):
    df.to_csv (f"../results/distance_table_{args.formula}.csv", index = False, header=True, float_format='%.4f')
else:
    df.to_csv (f"../results/distance_table_{args.formula}.csv", mode='a', index = False, header=False,float_format='%.4f')
