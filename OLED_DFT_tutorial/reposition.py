from pylada.crystal import read, write, Structure
import numpy.linalg as la
from copy import deepcopy
import sys

def is_same(A,B, tol = 1.e-5):
    C = reposition(A,B)
    for a in A:
        for c in C:
            if a.type != c.type or any(abs(a.pos - c.pos)>tol):
                break
        else:
            continue
        break
    else:
        return True
    return False

def reposition(A,B):

    A = deepcopy(A)

    iAcell = la.inv(A.cell)
    iBcell = la.inv(B.cell)
    
    for i, a in enumerate(A):
        distmin = None
        for j in [0,-1]:
            for k in [0,-1]:
                for l in [0,-1]:
                    dist = la.norm(a.pos + A.cell.dot([j,k,l]) - B[i].pos)
                    if distmin is None or dist < distmin:
                        distmin = dist
                        shift = [j,k,l]
    
        if shift != [0,0,0]:
            print("Shifting:", i, shift)
            
        a.pos = a.pos + A.cell.dot(shift)

    return A

if __name__ == "__main__":

    A = read.poscar(sys.argv[1])
    B = read.poscar(sys.argv[2])
    
    write.poscar(reposition(A,B), sys.argv[3])
    



