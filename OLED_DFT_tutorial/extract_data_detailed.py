from pylada.vasp import Extract
from pylada.crystal import read
from glob import iglob
import pymatgen.io.vasp as pyv
import numpy as np
import mcu
import re
import os
import pickle

output_name = "detailed_results.txt"

result_folder = "ene_all"

def get_soc(folder):

    A = read.poscar(folder+"/POSCAR")
    for i,a in enumerate(A):
        if a.type == "Ir":
            metal_id = i + 1
            break;

    with open(folder + "/OUTCAR") as f:
        s = f.read()

    return float(re.findall("Ion:\s+%d\s+E_soc:\s+([0-9\.\-]+)"%(metal_id), s)[0])

def get_barrier_alt(mol):
    ene_init = 0
    char_init = 0
    barriers = []
    ks = []
    for bond in sorted(iglob(mol + "/barrier/[0-9]/")):
        shape = []
        energies = []

        if ene_init != 0:
            shape.append(char_init)
            energies.append(0)

        for step in sorted(iglob(bond + "/*/")):
            prevcalc = Extract(step)
    
            s, mc = read_procar(step, "Ir")

            if ene_init == 0:
                ene_init = prevcalc.total_energy
            
            eig = np.array([j for i in prevcalc.eigenvalues[:,0,:] for j in i])
            occ = np.array([j for i in prevcalc.occupations[:,0,:] for j in i])

            o_eig = eig[occ > 0.3]
            o_mc = mc[occ > 0.3]

            char_LUMO = o_mc[np.argmax(o_eig)]

            shape.append(char_LUMO)
            energies.append(prevcalc.total_energy - ene_init)

            if char_LUMO  > 0.3:
                
                if len(shape) > 2:
                    # print((np.mean(np.array(shape)[1:-1]-np.array(shape)[:-2]))/(shape[-1]-shape[-2]), shape[-1]/shape[-2])
                    MC = (np.mean(np.array(shape)[1:-1]-np.array(shape)[:-2]))/(shape[-1]-shape[-2]) < 1/2 and shape[-1]/shape[-2] > 2
                elif len(shape) == 2:
                    MC = shape[-1]/shape[-2] > 2 
                else:
                    MC = True
                    
                if MC:
                    barriers.append(prevcalc.total_energy - ene_init)
                    break
        else:
            barriers.append(1000)

        char_init = shape[0]
        # print(shape)

        #print(np.mean(2*np.array(energies[1:min(4, len(energies))])/(np.arange(1,min(4, len(energies)))*1/6)**2), 2*np.array(energies[1:min(4, len(energies))])/(np.arange(1,min(4, len(energies)))*1/6)**2)

        ks.append(np.mean(2*np.array(energies[1:min(4, len(energies))])/(np.arange(1,min(4, len(energies)))*1/6)**2))

    if len(barriers) == 0:
        print("No MC found for %s"%(mol))
        return 0

    barrier = np.min(barriers)

    bond_type = "N"
    try:
        old_barrier, bond_type, rerelax, allB = get_barrier(mol)

        if abs(barrier - old_barrier) > 1e-12:
            print("Barrier changed for %s from %f to %f"%(mol, old_barrier, barrier))
        else:
            print("same")
    except:
        pass

    print("barrieres:", np.argmin(barriers), np.argmin(ks))
    print(ks)
    print(barriers)

    return barrier, np.min(ks), bond_type, rerelax, allB

def get_barrier(mol):
    with open(mol + "/barrier/BARRIER") as f:
        s = f.read()

    bond_types = re.findall("[0-9]\s+[0-9\.\-]+\s+([A-Z])\s*\n", s)

    print("BOND TYPES:", bond_types)

    return float(re.findall("Min barrier:\s*([0-9\.\-]+)", s)[0]), re.findall("Min barrier:\s*[0-9\.\-]+\s*(\w)", s)[0], re.search("Relaxing a new barrier!", s) is not None and "B" in bond_types, len(set(bond_types)) == 1 and bond_types[0] == "B"
    # return np.mean([float(a) for a in re.findall("[0-9]+\s+([0-9\.\-]+)\s+eV", s)])

def get_rel_barrier(mol):
    with open(mol + "/barrier/BARRIER") as f:
        s = f.read()

    return float(re.findall("New Min barrier:\s*([0-9\.\-]+)", s)[0])
    # return np.mean([float(a) for a in re.findall("[0-9]+\s+([0-9\.\-]+)\s+eV", s)])

def get_1b_barrier(mol):
    with open(mol + "/barrier/BARRIER") as f:
        s = f.read()

    return float(re.findall("New SB barrier:\s*([0-9\.\-]+)", s)[0])
    # return np.mean([float(a) for a in re.findall("[0-9]+\s+([0-9\.\-]+)\s+eV", s)])

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

all_times = dict()

with open(output_name,"w") as f:

    print("Entry", "1240/Ediff", 
              "k_r", "Ea (not realxed)", "Ea (relaxed)", "Ea (fully relaxed)", "[WEAK BOND]",
              "[Not N bond]",
              "[rerelaxed]",
              "[All carbenes]",
              file=f)

    for mol in iglob(result_folder + "/*"):

        sp_rel_nup = Extract(mol+"/SP_B3LYP_nup")
        nup2 = Extract(mol+"/Nupdown2_B3LYP")
        # bar = Extract(mol+"/barrier_B3LYP")
        optics_name = "/Optics_B3LYP"
        opt = Extract(mol + optics_name)

        if not nup2.success:
            print(mol, "Nupdown2 failed")
            continue
        if not sp_rel_nup.success:
            print(mol, "SP_B3LYP_nup failed")
            continue
            
        Ediff = float(nup2.total_energy-sp_rel_nup.total_energy)

        if Ediff < 2.35:
            print(mol, 1240/Ediff, file=f)
            continue

        if not opt.success:
            print(mol, 1240/Ediff, "NF", "NF", "NF", "NF", file=f)
            continue

        if not os.path.isfile(mol + "/barrier/BARRIER"):
            print(mol, 1240/Ediff, "NF", "NF", "NF", "NF", file=f)
            continue
        # if not bar.success:
        #     print(mol, "barrier_B3LYP failed")
        #     continue

        # Skip Open shell
        if opt.nelect%2 != 0:
            print(mol, "is open shell")
            continue

        # try:
        #     Ea = get_barrier_alt(mol)
        # except IndexError:
        #     print(mol, "Could not find info in BARRIER file")
        try:
            Ea, k, bond_type, rerelax, allB = get_barrier_alt(mol)
            print("Found old barrier (alt)")
        except:
            print(mol, 1240/Ediff, "NF", "NF", "NF", "NF", file=f)
            print(mol, "could not find old barrier")
            continue

        try:
            nEa = get_rel_barrier(mol)
        except IndexError:
            print(mol, "could not find NEW barrier")
            nEa = -1

        try:
            sbEa = get_1b_barrier(mol)
        except IndexError:
            print(mol, "could not find SB barrier")
            sbEa = -1

        # Eab3 = bar.total_energy - nup2.total_energy

        try:
            Esoc = get_soc(mol+optics_name)
        except IndexError:
            print(mol, "could not find SOC element")
            continue

        # Radiative Lifetime
        # s = np.fromregex(mol + "/Optics_b3lyp/PROCAR", "\s*([^\s]*)\n\s*\n\s*band", [("num", float)])
        # s = s.view((float, len(s.dtype.names)))

        s, mc = read_procar(mol + optics_name, "Ir")

        f_osc = mcu.read_WAVEDER(mol + optics_name + "/WAVEDER")

        idx_order = np.argsort(opt.eigenvalues[0,:])

        M = np.empty(3)
        E = np.empty(3)
        K = np.empty(3)

        if int(opt.nelect)%2 == 0:
        
            idx = np.array(range(opt.nbands))
            
            # print("Order is the same?", all(idx_order == idx))
            
            idx_up = idx[s[idx_order]>0]
            
            idx_down = idx[s[idx_order]<=0]
            
            # Down -> Down*
            i = idx_order[idx_down[idx_down < opt.nelect][-1]]
            j = idx_order[idx_down[idx_down >= opt.nelect][0]]
            dldh = f_osc[0][0,0,i,j,:]
            dldh_s = f_osc[0][0,0,j,i,:].conjugate()
            E_dldh = opt.eigenvalues[0,j] - opt.eigenvalues[0,i]
            # print(i,j,s[i],s[j])
            
            # Up -> Up*
            i = idx_order[idx_up[idx_up < opt.nelect][-1]]
            j = idx_order[idx_up[idx_up >= opt.nelect][0]]
            uluh = f_osc[0][0,0,i,j,:]
            uluh_s = f_osc[0][0,0,j,i,:].conjugate()
            E_uluh = opt.eigenvalues[0,j] - opt.eigenvalues[0,i]
            # print(i,j,s[i],s[j])
            
            # Down -> Up*
            i = idx_order[idx_down[idx_down < opt.nelect][-1]]
            j = idx_order[idx_up[idx_up >= opt.nelect][0]]
            # print(i,j,s[i],s[j])
            dluh = f_osc[0][0,0,i,j,:]
            dluh_s = f_osc[0][0,0,j,i,:].conjugate()
            E_dluh = opt.eigenvalues[0,j] - opt.eigenvalues[0,i]
            
            # Up -> Down*
            i = idx_order[idx_up[idx_up < opt.nelect][-1]]
            j = idx_order[idx_down[idx_down >= opt.nelect][0]]
            # print(i, j,s[i],s[j])
            uldh = f_osc[0][0,0,i,j,:]
            uldh_s = f_osc[0][0,0,j,i,:].conjugate()
            E_uldh = opt.eigenvalues[0,j] - opt.eigenvalues[0,i]            
            
            M[0] =  np.sum(np.abs(uldh)**2)
            E[0] =  E_uldh
            
            M[1] =  np.sum(np.abs(dluh)**2)
            E[1] =  E_dluh
            
            M[2] =  1/2 * np.sum(np.abs(dldh - uluh)**2)
            E[2] =  0.5 * (E_dldh + E_uluh)

        else:
            print(mol, "Odd number of electrons")
            for k in range(3):
                i = idx_order[int(opt.nelect) - 1]
                j = idx_order[int(opt.nelect) + 1 + k]
                M[k] = np.sum(np.abs(f_osc[0][0,0,i,j,:])**2)
                E[k] = opt.eigenvalues[0,j] - opt.eigenvalues[0,i]
            
        C = 3.79634e6
        Kb = 8.617e-5
        T = 300

        idmin = np.argmin(E)

        if Ediff < 1e-2:
            print(mol, "triplet ground state")
            Ediff = E[idmin]

        for i in range(3):
            # print(np.exp(-E[i]/(Kb*T))/np.sum(np.exp(-E/(Kb*T))))
            K[i] = C * (Ediff + E[i] - E[idmin])**3 * M[i]

        k_r = np.sum(K * np.exp(-E/(Kb*T)))/np.sum(np.exp(-E/(Kb*T)))

        k_nr = np.exp(-Ea/(Kb*T))

        timing = pickle.load(open(mol + "/timing.pkl","rb"))

        for key in timing:
            if key in all_times:
                all_times[key].append(timing[key])
            else:
                all_times[key] = [timing[key]]

        #print(Ediff, E, M, K, k_r, k_nr, "<======")

        print(mol, "->OK<-")
        #print(mol, Ediff, k_r, Esoc, Ea, sbEa, k, file=f)
        k_r = k_r/15
        knr_sb = 1e7*np.exp(-sbEa/3/0.02585) + np.exp(-Ediff)*1e5
        knr_n = 1e7*np.exp(-nEa/3/0.02585) + np.exp(-Ediff)*1e5
        knr = 6.251e12*np.exp(-Ea/2.4/0.02585) + np.exp(-Ediff)*1e5

        print(mol, 1240/Ediff, 
              k_r, Ea, nEa if nEa != -1 else "NF", 
              sbEa if sbEa != -1 else "NF",
              "" if k > 5 else "WEAK BOND",
              "" if bond_type is "N" else "Not N bond",
              "rerelaxed" if rerelax else "",
              "All carbenes" if allB else "",
              file=f)
    
pickle.dump(all_times, open(output_name + ".times","wb"))
