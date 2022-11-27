from pylada.vasp import Extract
from pylada.crystal import read
from glob import iglob
import pymatgen.io.vasp as pyv
import numpy as np
import re
import pickle

def read_procar(folder, metal):
    
    A = read.poscar(folder+"/POSCAR")
    for i,a in enumerate(A):
        if a.type == metal:
            metal_id = i + 1
            break;

    data = []
    try: 
        with open(folder+"/PROCAR","r") as f:
            for line in f:
                if "band " in line:
                    data.append([None, None, None])
                    tot = 0
                    continue
                match = re.match("^\s*%d.*\s([\-0-9.]+)\s*$"%metal_id, line)
                if match and data[-1][0] is None:
                    data[-1][0] = float(match.group(1))
                match = re.match("^\s*tot.*\s([\-0-9.]+)\s*$", line)
                if match:
                    tot += 1
                    if tot == 1:
                        data[-1][1] = float(match.group(1))
                    if tot == 4:
                        data[-1][2] = float(match.group(1))
    except FileNotFoundError:
        return None, None

    data = np.array(data)

    s = data[:,2]

    c = data[:,0]/data[:,1]

    return s, c
            

def bandgap(obj):
    return obj.eigenvalues[0,int(obj.nelect)] - obj.eigenvalues[0,int(obj.nelect)-1]

plt_data = []
for mol in iglob("ene_all/doi-10-1002-anie-200604733-entry-picjuh-molecule-0/barrier/0/00*"):

    # Non Radiative Lifetime
    s, mc = read_procar(mol, "Ir")

    soc = Extract(mol)

    print(mol, soc.fermi_energy, soc.total_energy, soc.total_energies[0])

    plt_data.append([soc.eigenvalues[:,0,:], mc, soc.fermi_energy])

pickle.dump(plt_data,open("images_weird.dat","wb"))
