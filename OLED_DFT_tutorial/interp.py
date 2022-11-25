from pylada.crystal import read, write, Structure, supercell
import numpy.linalg as la
import numpy as np
import os
import sys

def interp(A,B,n_steps):

    iAcell = la.inv(A.cell)
    iBcell = la.inv(B.cell)
    
    T = B.cell.dot(la.inv(A.cell))
    
    disps = []
    for i,a in enumerate(A):
        disps.append(iBcell.dot(B[i].pos) - iAcell.dot(a.pos))
    
    for i in range(n_steps):
        curMat = (T-np.eye(3))*i/(n_steps-1) + np.eye(3)
        curStruc = Structure(curMat.dot(A.cell))
        for j,d in enumerate(disps):
            curDisp = d*i/(n_steps-1)
            curPos = curStruc.cell.dot(iAcell.dot(A[j].pos) + curDisp)
            curStruc.add_atom(*curPos,A[j].type)
            curStruc.name = "%d"%(i/(n_steps-1))
        yield curStruc

if __name__ == "__main__":
    A = read.poscar("vtst_images_new/00/POSCAR")
    B = read.poscar("vtst_images_new/10/POSCAR")
    
    outdir = "vtst_images_new"
    
    for i, curStruc in enumerate(interp(A,B,11)):
        os.makedirs(outdir+"/%02d"%i, exist_ok=True)
        write.poscar(curStruc, vasp5=True, file=outdir+"/%02d/POSCAR"%i)
