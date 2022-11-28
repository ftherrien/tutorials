from pylada.crystal import read, write, Structure
import numpy as np
import sys

A = read.poscar(sys.argv[1])

for i,a in enumerate(A):
    if a.type in ["Ir","Pt","Pd","Au","Rh"]:
        pos = a.pos
        break
else:
    print("Not in list")

allpos = []
for i,a in enumerate(A):
    a.pos = a.pos - pos
    allpos.append(a.pos)
    
allpos = np.array(allpos)

print(np.max(allpos,axis=0) - np.min(allpos,axis=0))

A.cell = np.diag(np.max(allpos,axis=0) - np.min(allpos,axis=0) + 10)

write.poscar(A,sys.argv[1])
