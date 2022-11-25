from pylada.crystal import read, write, Structure
import numpy as np
import numpy.linalg as la
import sys

extra_length = 1.5
shell = 2.5
tol = 1e-5

def list_bonds(A,idx, bond_length):
    list_atoms = []
    for i,a in enumerate(A):
        if i!=idx and la.norm(a.pos - A[idx].pos) < bond_length:
            list_atoms.append(i)
    return list_atoms

def to_move(A,idx,i,bond_length):
    for j in list_bonds(A,i,bond_length):
        if j not in idx:
            idx.append(j)
            idx = to_move(A,idx,j,bond_length)

    return idx

def stretch(A, extra_length, IrN=2.25, shell = 2, tol = 1e-5, with_bond=False):
    n_bonds = 0
    bond = np.zeros(3)
    bonds = []
    for i, a in enumerate(A):
        if a.type == "Ir":
            for j in list_bonds(A,i,IrN):
                if True: # A[j].type == "N":
                    for b in bonds:
                        if abs(la.norm(b) - la.norm(A[j].pos - a.pos)) < tol:
                            print("NO!", la.norm(b), la.norm(A[j].pos - a.pos))
                            break
                    else:
                        n_bonds += 1
                        bond = A[j].pos - a.pos
                        print("Bond length:", la.norm(bond))
                        bonds.append(bond)
                        new_bond = bond * (la.norm(bond) + extra_length)/la.norm(bond)
                        shift = new_bond - bond
                        atoms_to_move = to_move(A,[i,j],j,shell)
                        print("Ir-N", atoms_to_move, len(atoms_to_move))
                        B = Structure(A.cell)
                        for k,b in enumerate(A):
                            if k in atoms_to_move[1:]:
                                B.add_atom(*(b.pos + shift), b.type)
                            else:
                                B.add_atom(*b.pos, b.type)
                        if with_bond:
                            yield B, i, j
                        else:
                            yield B

if __name__ == "__main__":
    A = read.poscar(sys.argv[1])
    
    for i,B in enumerate(stretch(A, 1.5)):
        write.poscar(B, sys.argv[1]+"_stretch_%d"%n_bonds)
