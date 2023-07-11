import sys
import argparse
import numpy as np
import pandas as pd
import re
from pymatgen.io.cif import CifParser
from sklearn.metrics import mean_squared_error, mean_absolute_error
from pymatgen.analysis.structure_matcher import StructureMatcher as SM
from cifparser import CifParserExpand

parser = argparse.ArgumentParser(
        description=(
            "rmsd,mae,rms and accuracy"
            "calculate rmsd,mae,rms and contact map accuracy"
        )
    )

parser.add_argument(
    "--cif",
    type=str,
    #default="2-14-mp-236.cif",
    default="dataset/SrTiO3_target.cif",
    metavar="PATH",
    help="Path to target cif file",
)

parser.add_argument(
    "--predicted",
    type=str,
    #default="2-14-mp-236_predicted.cif",
    default="dataset/SrTiO3_pred.cif",
    metavar="PATH",
    help="Path to predicted cif file",
)

args = parser.parse_args(sys.argv[1:])


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
                if e1 + z1 in yuzhi:
                    max_len = float(yuzhi[yuzhi.rindex(e1 + z1) + len(e1 + z1) + 15:yuzhi.rindex(
                        e1 + z1) + len(e1 + z1) + 22])
                    min_len = float(yuzhi[yuzhi.rindex(e1 + z1) + len(e1 + z1) + 4:yuzhi.rindex(
                        e1 + z1) + len(e1 + z1) + 11])
                elif e2 + z2 in yuzhi:
                    max_len = float(yuzhi[
                                    yuzhi.rindex(e2 + z2) + len(e2 + z2) + 15:yuzhi.rindex(
                                        e2 + z2) + len(e2 + z2) + 22])
                    min_len = float(yuzhi[
                                    yuzhi.rindex(e2 + z2) + len(e2 + z2) + 4:yuzhi.rindex(
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
    yuzhi = fileobj.read()

with open(args.cif, "r") as fileobj:
    text1 = fileobj.read()
    co = re.findall(r'\d+\.\d+', text1[text1.rindex('occupancy'):])
    cof1 = []
    for i in co:
        cof1.append(float(i))
    cof1=np.array(cof1).reshape(int(len(cof1)/3),3)
    zs=text1[text1.rindex('_atom_site_occupancy') + len('_atom_site_occupancy') + 3:]
    mpl=zs.split('  ')
    mpl=np.array(mpl).reshape(int(len(mpl) / 7), 7)
    mulpl=[]
    for u in range(len(mpl)):
        mulpl.append(mpl[u][2])

with open(args.predicted, "r") as fileobj:
    textpr = fileobj.read()
    copr = re.findall(r'\d+\.\d+', textpr[textpr.rindex('occupancy'):])
    coprf1 = []
    for i in copr:
        coprf1.append(float(i))
    coprf1=np.array(coprf1).reshape(int(len(coprf1)/3),3)

parser1 = CifParserExpand(args.cif)
structure1 = parser1.get_structures()[0]

try:
    parser2 = CifParserExpand(args.predicted)
    print (parser2)
    structure2 = parser2.get_structures()[0]
    d2=structure2.distance_matrix
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
            pass
        except ZeroDivisionError:
            fitness = 9995



sm = SM(ltol=0.6, stol=0.6,angle_tol=20)
try:
    rms = sm.get_rms_dist(structure1, structure2)[0]
except TypeError:
    rms='None'

print (cof1)

print('******')
      
print (coprf1)

print (len(cof1),len(coprf1))
      
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

    min_rmse = np.sqrt(mean_squared_error(cof1.flatten(), np.array(co4).flatten()))
    min_mae = mean_absolute_error(cof1.flatten(), np.array(co4).flatten())
else:
    min_rmse=5555
    min_mae=5555


fp = open(args.predicted.replace('.cif','_rmsd'), "w")
fp.truncate()
with open(args.predicted.replace('.cif','_rmsd'), "w") as fw:
    fw.write(args.predicted+'\n')
    fw.write('rmsd:'+str(min_rmse)+'\n')
    fw.write('mae:'+str(min_mae)+'\n')
    fw.write('rms:' + str(rms) + '\n')



