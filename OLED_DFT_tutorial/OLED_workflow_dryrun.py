######################################
class CustomChain(object):

    # Defining the folder tree structure
    def __init__(self,vaspobj=None):
        self.names=['Nupdown2', 'Nupdown2_B3LYP', 'SP_B3LYP_nup', "Optics_B3LYP", "barrier"]
        self.vasp=vaspobj

    # Defining the Extraction object
    def Extract(self, jobdir):
        from pylada.vasp import MassExtract
        from pylada.vasp.extract import Extract
        extract = MassExtract(jobdir)
        success={}
        for name in self.names:
            success[name] = Extract(jobdir+'/'+name).success
        success = all(success.values())
        extract.success=success
        return extract
        
    # Creating the workflow
    def __call__(self, structure, outdir=None, **kwargs ):

        from copy import deepcopy
        from os import getcwd, makedirs
        from os.path import join, isdir, isfile
        from pylada.misc import RelativePath
        from pylada.error import ExternalRunFailed
        from pylada.vasp.extract import Extract
        from pylada.vasp import Vasp
        from pylada.vasp.relax import Relax
        from stretch import stretch
        from interp import interp
        from extract_images import read_procar
        from pylada.crystal import read, write
        from reposition import reposition, is_same
        import numpy as np
        import time
        import pickle
        #from constraints import stretchcombo
        from ase.constraints import FixBondLengths, FixBondLength
        from ase.io import read as aseread
        import ase.calculators.vasp as vasp_calculator
        from ase.optimize import BFGS
        import shutil
        from glob import iglob

        from rdkit import Chem
        
        def is_carbene(idx, outdir):
            try:
                entry = outdir.split("entry-")[1].split("-")[0].lower()
            except IndexError:
                entry = outdir.split("/")[-1].lower()
            
            mol = Chem.MolFromMolFile("./dry_run_mols/"+entry+".mol", sanitize = False, removeHs = False)
            
            # atom = mol.GetAtomWithIdx(idx)

            atoms = {}
            for atom in mol.GetAtoms():
                if atom.GetSymbol() in atoms:
                    atoms[atom.GetSymbol()].append(atom)
                else:
                    atoms[atom.GetSymbol()] = [atom]
        
            mol_list = sum(atoms.values(), [])

            atom = mol_list[idx]

            print("Carbene:", idx, atom.GetSymbol(), sum([b.GetBondTypeAsDouble() for b in atom.GetBonds()]))

            return sum([b.GetBondTypeAsDouble() for b in atom.GetBonds()]) == 3 and atom.GetSymbol() == "C"


        # make this function stateless.
        structure_ = structure.copy()
        outdir = getcwd() if outdir is None else RelativePath(outdir).path

        try:
            timing = pickle.load(open(outdir + "/timing.pkl","rb"))
        except FileNotFoundError:
            timing = dict()

        ############ Part 1: Nupdown = 2 ############### 
        name  = self.names[0]

        nupdown = Relax(copy=deepcopy(self.vasp))
        nupdown.relaxation = "ionic"
        nupdown.nsw        = 200
        nupdown.add_keyword("NUPDOWN", 2)
        nupdown.ispin = 2
        nupdown.lwave = True
        
        params = deepcopy(kwargs)
        fulldir = join(outdir, name)
        
        t = time.time()
        ## if this calculation has not been done, run it
        outputnup = nupdown(structure_, outdir=fulldir, convergence=1e-3, **params)
        if not outputnup.success: 
            raise ExternalRunFailed("NUPDOWN2 calculation did not succeed.")

        if self.names[0] not in timing:
            timing[self.names[0]] = time.time() - t

        pickle.dump(timing, open(outdir + "/timing.pkl","wb"))

        ############ Part 2: Nupdown2 + B3LYP ###############  
        name  = self.names[1]

        nupdown_b3lyp = Vasp(copy=deepcopy(nupdown))
        nupdown_b3lyp.relaxation = "static"
        nupdown_b3lyp.nsw        = 0
        nupdown_b3lyp.algo = "Damped"
        nupdown_b3lyp.add_keyword("LHFCALC", True)
        nupdown_b3lyp.add_keyword("GGA", "B3")    
        nupdown_b3lyp.add_keyword("AEXX", 0.2)   
        nupdown_b3lyp.add_keyword("AGGAX", 0.72)  
        nupdown_b3lyp.add_keyword("AGGAC", 0.81)  
        nupdown_b3lyp.add_keyword("ALDAC", 0.19)  
        
        params = deepcopy(kwargs)
        fulldir = join(outdir, name)
        
        t = time.time()
        ## if this calculation has not been done, run it
        output = nupdown_b3lyp(structure_, outdir=fulldir, restart=outputnup, **params)
        if not output.success: 
            raise ExternalRunFailed("NUPDOWN_B3LYP2 calculation did not succeed.")

        high = output.total_energy

        if self.names[1] not in timing:
            timing[self.names[1]] = time.time() - t

        pickle.dump(timing, open(outdir + "/timing.pkl","wb"))

        ############ Part 3: SP + B3LYP from NUP ###############  
        name  = self.names[2]
        
        sp_b3lyp_nup = Vasp(copy=deepcopy(self.vasp))
        sp_b3lyp_nup.relaxation = "static"
        sp_b3lyp_nup.nsw        = 0 
        sp_b3lyp_nup.ispin = 2
        sp_b3lyp_nup.algo = "Damped"
        sp_b3lyp_nup.add_keyword("LHFCALC", True)
        sp_b3lyp_nup.add_keyword("GGA", "B3")    
        sp_b3lyp_nup.add_keyword("AEXX", 0.2)   
        sp_b3lyp_nup.add_keyword("AGGAX", 0.72)  
        sp_b3lyp_nup.add_keyword("AGGAC", 0.81)  
        sp_b3lyp_nup.add_keyword("ALDAC", 0.19)  
        
        params = deepcopy(kwargs)
        fulldir = join(outdir, name)
        
        t = time.time()
        ## if this calculation has not been done, run it
        output = sp_b3lyp_nup(structure_, outdir=fulldir, restart=outputnup, **params)
        if not output.success: 
            raise ExternalRunFailed("SP_B3LYP_NUP calculation did not succeed.")

        low = output.total_energy

        if self.names[2] not in timing:
            timing[self.names[2]] = time.time() - t

        pickle.dump(timing, open(outdir + "/timing.pkl","wb"))


        ################# Continue if blue #####################
        if float(high - low) < 2.35:
            print("Not Blue! - Stopping", float(high - low))
            return self.Extract(fulldir)

        ############ Part 4: Spin orbit coupling ###############  
        name  = self.names[3]

        soc = Vasp(copy=deepcopy(self.vasp))
        sp_b3lyp_nup.relaxation = "static"
        sp_b3lyp_nup.nsw        = 0 
        sp_b3lyp_nup.ispin = 2
        soc.add_keyword("LSORBIT", True)
        soc.nbands = output.nbands*2
        soc.loptics=True
        
        params = deepcopy(kwargs)
        fulldir = join(outdir, name)
        
        t = time.time()
        output = soc(output.structure, outdir=fulldir, **params)
        if not output.success: 
            raise ExternalRunFailed("SOC calculation did not succeed.")

        if self.names[3] not in timing:
            timing[self.names[3]] = time.time() - t

        pickle.dump(timing, open(outdir + "/timing.pkl","wb"))

        ############ Part 5: Barrier calculation ###############

        t = time.time()

        max_steps = 7
        min_char = 0.3

        A = outputnup.structure

        A = reposition(A, structure)

        makedirs(outdir + "/" + self.names[4], exist_ok = True)

        with open(outdir + "/" + self.names[4] + "/BARRIER", "w") as f:

            barriers = []
            calcs = []
            bond_indices = []
            ene_init = 0
            bond_types = []
            for i, (A_stretched, bond_i, bond_j) in enumerate(stretch(A, 1, IrN = 2.5, shell = 1.9, with_bond = True)):

                attachment = bond_i if A[bond_j].type == "Ir" else bond_j

                print(is_carbene(attachment, outdir))

                type_bond = "B" if is_carbene(attachment, outdir) else A[attachment].type
                
                makedirs(outdir + "/" + self.names[4] + "/%d"%i, exist_ok = True)
                shape = []
                energies = []
            
                if ene_init != 0:
                    shape.append(char_init)
                    energies.append(0)
            
                for j, A_step in enumerate(interp(A, A_stretched, max_steps)):
                    
                    if j == 0 and ene_init != 0:
                        continue

                    if j == 4:
                        if type_bond == "C":
                            print("Stopped!", file=f)
                            break

                    curcalc = Vasp(copy=deepcopy(self.vasp))
                    curcalc.add_keyword("NUPDOWN", 2)
                    curcalc.ispin = 2
                    params = deepcopy(kwargs)
                    fulldir = outdir + "/" + self.names[4] + "/%d"%i + "/%02d"%j

                    if isdir(outdir + "/barrier_disordered"):
                        if j == 0:
                            try:
                                shutil.copytree(outdir + "/barrier_disordered/0/00", fulldir)
                            except FileExistsError:
                                pass
                        if j == 1:
                            already_done = []
                            for precalc_barrier in sorted(iglob(outdir + "/barrier_disordered/[0-9]*")):
                                pidx = int(precalc_barrier.split("/")[-1])
                                if is_same(A_step, read.poscar(precalc_barrier + "/01/POSCAR")):
                                    already_done.append(pidx)
                            if len(already_done) > 0:
                                if i in already_done:
                                    pidx = i
                                else:
                                    pidx = already_done[0]
                                
                                try:
                                    shutil.copytree(outdir + "/barrier_disordered/" + pidx, outdir + "/" + self.names[4] + "/%d"%i)
                                except FileExistsError:
                                    pass
        
                    prevcalc = curcalc(A_step, outdir=fulldir, **params)
                    if not prevcalc.success: 
                        raise ExternalRunFailed("Barrier calc failed at %d, %d"%(i,j))

                    if ene_init == 0:
                        ene_init = float(prevcalc.total_energy)
            
                    s, mc = read_procar(fulldir, "Ir")
            
                    eig = np.array([j for i in prevcalc.eigenvalues[:,0,:] for j in i])
                    occ = np.array([j for i in prevcalc.occupations[:,0,:] for j in i])
            
                    o_eig = eig[occ > 0.3]
                    o_mc = mc[occ > 0.3]
            
                    char_LUMO = o_mc[np.argmax(o_eig)]
            
                    shape.append(char_LUMO)
                    energies.append(float(prevcalc.total_energy) - ene_init)
                    k = np.mean(2*np.array(energies[1:min(4, len(energies))])/(np.arange(1,min(4, len(energies)))*1/6)**2)

                    if char_LUMO > 0.3:
                        
                        if len(shape) > 2:
                            MC = (np.mean(np.array(shape)[1:-1]-np.array(shape)[:-2]))/(shape[-1]-shape[-2]) < 1/2
                        elif len(shape) == 2:
                            MC = shape[-1]/shape[-2] > 2 
                        else:
                            MC = True
                            
                        if MC:
                            barriers.append(float(prevcalc.total_energy) - ene_init)
                            calcs.append(prevcalc)
                            bond_indices.append([bond_i, bond_j])
                            bond_types.append(type_bond)
                            if isdir(outdir + "/barrier_disordered") and len(already_done) > 0 and pidx == i:
                                same_staring_pos = []
                                for precalc_arpess in sorted(iglob(outdir + "/barrier_disordered/ARPESS/*/")):
                                    if is_same(prevcalc.structure, read.poscar(precalc_arpess + "/POSCAR")):
                                        same_staring_pos.append(precalc_arpess)
                                if len(same_staring_pos) == 1:
                                    print("COPYING FROM DISORDERED!!")
                                    shutil.copytree(same_staring_pos[0], outdir + "/" + self.names[4] + "/ARPESS/%d"%i)
                                else:
                                    raise RuntimeError("Multiple identical barrier starting points!")
                            break
    
                char_init = shape[0]
                print("Character:", shape, file=f)
                print("Energies:", energies, file=f)
                print("Bond strength:", k, file=f)
                print("Bond type:", type_bond, file=f)

                if k < 3:
                    print("Very weak bond!", file=f)
                    return self.Extract(fulldir)

                if j==0:
                    print("Excited state is MC!", file=f)
                    break

            if len(barriers) == 0:
                print("No MC found", file=f)
            else:
                print("Barriers (eV)", file=f)
                for i,b in enumerate(barriers):
                    print(i, b, bond_types[i], file=f)
                    
                print(file=f)
                minidx = np.argmin(barriers)
                print("Min barrier:", barriers[minidx], bond_types[minidx], file=f)

            if self.names[4] not in timing:
                timing[self.names[4]] = time.time() - t
            
            pickle.dump(timing, open(outdir + "/timing.pkl","wb"))
            

           ############ Part 5.3: Relax Barrier 1 bond only ###############

            print("More RELAXATION", bond_indices)

            t = time.time()
            
            new_barriers = []
            old_calcdir = calcdir
            calcdir = outdir + "/" + self.names[4] + "/SINGLE_BOND/%d"%i
            makedirs(calcdir, exist_ok = True)
            shutil.copyfile(old_calcdir+ "/POSCAR_final", calcdir + "/POSCAR_init")
            A = read.poscar(calcdir + "/POSCAR_init")

            if not isfile(calcdir + "/POSCAR_final"):

                if not isdir(calcdir + "/asedir"):
                    atoms = aseread(calcdir + "/POSCAR_init")
                else:
                    atoms = aseread(calcdir + "/asedir/POSCAR")
                
                liste=[[bond_i,bond_j,1.0]]
                
                magmom = np.zeros(len(atoms))
                magmom[0] = 0.6
                
                initial_values = [np.linalg.norm(A[bi].pos - A[bj].pos) for bi, bj in bond_indices]
                
                #cc=stretchcombo(initial_value, liste, atoms)
                
                cc=FixBondLength(bond_i, bond_j)
                
                print("INITIAL VALUES:", initial_values)
                # print(cc.return_adjusted_prositions())
                
                # atoms.set_positions(cc.return_adjusted_prositions())
                
                atoms.set_constraint([cc])

                print("AVANT", file=f)
                calc = vasp_calculator.Vasp(command='srun vasp_std',
                                            xc = 'PBE',
                                            # ivdw=11,
                                            kpts  = (1,1,1),
                                            npar = 5,
                                            istart = False,
                                            gamma = True,
                                            ismear = 0,
                                            nelm = 250,
                                            nupdown = 2,
                                            algo = 'Fast',
                                            sigma = 0.2,
                                            ibrion = -1,
                                            ediffg = 0.001,
                                            ediff = 1e-5,
                                            prec = 'Accurate',
                                            lcharg = False,
                                            lwave = False,
                                            lorbit = 10,
                                            nsw = 0,
                                            lreal = False,
                                            ispin = 2,
                                            lmaxmix = 4,
                                            magmom = magmom,
                                            label="image",
                                            directory=calcdir + "/asedir")

                print("APRES", file=f)

                atoms.set_calculator(calc)
                
                dyn = BFGS(atoms, restart=calcdir + "/asedir/restart.pkl")
                dyn.run(fmax=0.05)
            
                new_barriers.append(float(atoms.get_potential_energy())-ene_init)
            
                atoms.write(calcdir + "/POSCAR_final")
            
                print("ENERGY:", atoms.get_potential_energy(), file = f)
            
            else:
                asecalc = Extract(calcdir + "/asedir")
                new_barriers.append(float(asecalc.total_energy)-ene_init)

            B = read.poscar(calcdir + "/POSCAR_final")

            final_values = [np.linalg.norm(B[bi].pos - B[bj].pos) for bi, bj in bond_indices]

            print("SINGLE BOND", barriers, file=f)
            if len(barriers) != 0:
                #with open(outdir + "/" + self.names[4] + "/BARRIER", "a") as f:
                print("SB Barriers (eV)", file=f)
                for i,b in enumerate(new_barriers):
                    print(i, b, file=f)
                    
                print(file=f)
                new_minidx = np.argmin(new_barriers)
                print("New SB barrier:", new_barriers[new_minidx], file=f)
                
                if minidx != new_minidx:
                    print("Not the same barrier!", file=f)
            
            print("CHANGE IN LENGTH:", (np.array(initial_values) - np.array(final_values))/np.array(initial_values), file=f)
            
            if "SINGLE_BOND" not in timing:
                timing["SINGLE_BOND"] = time.time() - t
            
            pickle.dump(timing, open(outdir + "/timing.pkl","wb"))


        return self.Extract(fulldir)
