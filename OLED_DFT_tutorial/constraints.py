from math import sqrt
from ase.utils.geometry import find_mic

import numpy as np

__all__ = ['FixCartesian', 'FixBondLength', 'FixedMode', 'FixConstraintSingle',
           'FixAtoms', 'UnitCellFilter', 'FixScaled', 'StrainFilter',
           'FixedPlane', 'Filter', 'FixConstraint', 'FixedLine',
           'FixBondLengths', 'FixInternals', 'Hookean','stretchcombo']


def dict2constraint(dct):
    if dct['name'] not in __all__:
        raise ValueError
    return globals()[dct['name']](**dct['kwargs'])

            
def slice2enlist(s, n):
    """Convert a slice object into a list of (new, old) tuples."""

    print(s,n)

    if isinstance(s, (list, tuple)):
        return enumerate(s)
    return enumerate(range(*s.indices(n)))


class FixConstraint:
    """Base class for classes that fix one or more atoms in some way."""

    def index_shuffle(self, atoms, ind):
        """Change the indices.

        When the ordering of the atoms in the Atoms object changes,
        this method can be called to shuffle the indices of the
        constraints.

        ind -- List or tuple of indices.

        """
     #   raise NotImplementedError
        """Shuffle the indices of the two atoms in this constraint"""
        print('##',self,ind)
        newa = [-1, -1, -1]  # Signal error
        for new, old in slice2enlist(ind, len(atoms)):
            for i, a in enumerate(self.indices):
                if old == a:
                    newa[i] = new
        if newa[0] == -1 or newa[1] == -1 or newa[2] == -1:
            raise IndexError('Constraint not part of slice')
        self.indices = newa

    def repeat(self, m, n):
        """ basic method to multiply by m, needs to know the length
        of the underlying atoms object for the assignment of
        multiplied constraints to work.
        """
        msg = ("Repeat is not compatible with your atoms' constraints."
               ' Use atoms.set_constraint() before calling repeat to '
               'remove your constraints.')
        raise NotImplementedError(msg)

    def adjust_momenta(self, atoms, momenta):
        """Adjusts momenta in identical manner to forces."""
        self.adjust_forces(atoms, momenta)

    def copy(self):
        return dict2constraint(self.todict().copy())
        

class FixConstraintSingle(FixConstraint):
    """Base class for classes that fix a single atom."""

    def index_shuffle(self, atoms, ind):
        """The atom index must be stored as self.a."""
        newa = -1   # Signal error
        for new, old in slice2enlist(ind, len(atoms)):
            if old == self.a:
                newa = new
                break
        if newa == -1:
            raise IndexError('Constraint not part of slice')
        self.a = newa


class FixAtoms(FixConstraint):
    """Constraint object for fixing some chosen atoms."""
    def __init__(self, indices=None, mask=None):
        """Constrain chosen atoms.

        Parameters
        ----------
        indices : list of int
           Indices for those atoms that should be constrained.
        mask : list of bool
           One boolean per atom indicating if the atom should be
           constrained or not.

        Examples
        --------
        Fix all Copper atoms:

        >>> mask = [s == 'Cu' for s in atoms.get_chemical_symbols()]
        >>> c = FixAtoms(mask=mask)
        >>> atoms.set_constraint(c)

        Fix all atoms with z-coordinate less than 1.0 Angstrom:

        >>> c = FixAtoms(mask=atoms.positions[:, 2] < 1.0)
        >>> atoms.set_constraint(c)
        """

        if indices is None and mask is None:
            raise ValueError('Use "indices" or "mask".')
        if indices is not None and mask is not None:
            raise ValueError('Use only one of "indices" and "mask".')

        if mask is not None:
            indices = np.arange(len(mask))[np.asarray(mask, bool)]
        else:
            # Check for duplicates:
            srt = np.sort(indices)
           # print 'srt',srt
           # print 'indices',indices
            assert (sorted(srt) == sorted(indices))
            #assert (srt == indices).all()
            for i in range(len(indices) - 1):
                if srt[i] == srt[i + 1]:
                    raise ValueError(
                        'FixAtoms: The indices array contained duplicates. '
                        'Perhaps you wanted to specify a mask instead, but '
                        'forgot the mask= keyword.')
        self.index = np.asarray(indices, int)

        if self.index.ndim != 1:
            raise ValueError('Wrong argument to FixAtoms class!')

    def adjust_positions(self, atoms, new):
        new[self.index] = atoms.positions[self.index]

    def adjust_forces(self, atoms, forces):
        forces[self.index] = 0.0

    def index_shuffle(self, atoms, ind):
        # See docstring of superclass
        index = []
        for new, old in slice2enlist(ind, len(atoms)):
            if old in self.index:
                index.append(new)
        if len(index) == 0:
            raise IndexError('All indices in FixAtoms not part of slice')
        self.index = np.asarray(index, int)

    def __repr__(self):
        return 'FixAtoms(indices=%s)' % ints2string(self.index)

    def todict(self):
        return {'name': 'FixAtoms',
                'kwargs': {'indices': self.index}}

    def repeat(self, m, n):
        i0 = 0
        natoms = 0
        if isinstance(m, int):
            m = (m, m, m)
        index_new = []
        for m2 in range(m[2]):
            for m1 in range(m[1]):
                for m0 in range(m[0]):
                    i1 = i0 + n
                    index_new += [i + natoms for i in self.index]
                    i0 = i1
                    natoms += n
        self.index = np.asarray(index_new, int)
        return self

    def delete_atom(self, ind):
        """ Removes atom number ind from the index array, if present.
        Required for removing atoms with existing FixAtoms constraints.
        """
        if ind in self.index:
            i = list(self.index).index(ind)
            self.index = np.delete(self.index, i)
        for i in range(len(self.index)):
            if self.index[i] >= ind:
                self.index[i] -= 1


def ints2string(x, threshold=None):
    """Convert ndarray of ints to string."""
    if threshold is None or len(x) <= threshold:
        return str(x.tolist())
    return str(x[:threshold].tolist())[:-1] + ', ...]'


class FixBondLengths(FixConstraint):
    def __init__(self, pairs, iterations=10):
        self.constraints = [FixBondLength(a1, a2) for a1, a2 in pairs]
        self.iterations = iterations

    def adjust_positions(self, atoms, new):
        for i in range(self.iterations):
            for constraint in self.constraints:
                constraint.adjust_positions(atoms, new)

    def adjust_forces(self, atoms, forces):
        for i in range(self.iterations):
            for constraint in self.constraints:
                constraint.adjust_forces(atoms, forces)

    def todict(self):
        return {'name': 'FixBondLengths',
                'kwargs': {'pairs': [constraint.indices
                                     for constraint in self.constraints],
                           'iterations': self.iterations}}


class FixBondLength(FixConstraint):
    """Constraint object for fixing a bond length."""
    def __init__(self, a1, a2):
        """Fix distance between atoms with indices a1 and a2. If mic is
        True, follows the minimum image convention to keep constant the
        shortest distance between a1 and a2 in any periodic direction.
        atoms only needs to be supplied if mic=True.
        """
        self.indices = [a1, a2]
        self.constraint_force = None

    def adjust_positions(self, atoms, new):
        p1, p2 = atoms.positions[self.indices]
        d, p = find_mic(np.array([p2 - p1]), atoms._cell, atoms._pbc)
        q1, q2 = new[self.indices]
        d, q = find_mic(np.array([q2 - q1]), atoms._cell, atoms._pbc)
        d *= 0.5 * (p - q) / q
        new[self.indices] = (q1 - d[0], q2 + d[0])

    def adjust_forces(self, atoms, forces):
        d = np.subtract.reduce(atoms.positions[self.indices])
        d, p = find_mic(np.array([d]), atoms._cell, atoms._pbc)
        d = d[0]
        d *= 0.5 * np.dot(np.subtract.reduce(forces[self.indices]), d) / p**2
        self.constraint_force = d
        forces[self.indices] += (-d, d)

    def index_shuffle(self, atoms, ind):
        """Shuffle the indices of the two atoms in this constraint"""
        newa = [-1, -1]  # Signal error
        for new, old in slice2enlist(ind, len(atoms)):
            for i, a in enumerate(self.indices):
                if old == a:
                    newa[i] = new
        if newa[0] == -1 or newa[1] == -1:
            raise IndexError('Constraint not part of slice')
        self.indices = newa

    def get_constraint_force(self):
        """Return the (scalar) force required to maintain the constraint"""
        return self.constraint_force

    def __repr__(self):
        return 'FixBondLength(%d, %d)' % tuple(self.indices)

    def todict(self):
        return {'name': 'FixBondLength',
                'kwargs': {'a1': self.indices[0], 'a2': self.indices[1]}}

        
class FixedMode(FixConstraint):
    """Constrain atoms to move along directions orthogonal to
    a given mode only."""

    def __init__(self, mode):
        self.mode = (np.asarray(mode) / np.sqrt((mode**2).sum())).reshape(-1)

    def adjust_positions(self, atoms, newpositions):
        newpositions = newpositions.ravel()
        oldpositions = atoms.positions.ravel()
        step = newpositions - oldpositions
        newpositions -= self.mode * np.dot(step, self.mode)

    def adjust_forces(self, atoms, forces):
        forces = forces.ravel()
        forces -= self.mode * np.dot(forces, self.mode)

    def index_shuffle(self, atoms, ind):
        eps = 1e-12
        mode = self.mode.reshape(-1, 3)
        excluded = np.ones(len(mode), dtype=bool)
        excluded[ind] = False
        if (abs(mode[excluded]) > eps).any():
            raise IndexError('All nonzero parts of mode not in slice')
        self.mode = mode[ind].ravel()

    def todict(self):
        return {'name': 'FixedMode',
                'kwargs': {'mode': self.mode}}

    def __repr__(self):
        return 'FixedMode(%s)' % self.mode.tolist()


class FixedPlane(FixConstraintSingle):
    """Constrain an atom index *a* to move in a given plane only.

    The plane is defined by its normal vector *direction*."""

    def __init__(self, a, direction):
        self.a = a
        self.dir = np.asarray(direction) / sqrt(np.dot(direction, direction))

    def adjust_positions(self, atoms, newpositions):
        step = newpositions[self.a] - atoms.positions[self.a]
        newpositions[self.a] -= self.dir * np.dot(step, self.dir)

    def adjust_forces(self, atoms, forces):
        forces[self.a] -= self.dir * np.dot(forces[self.a], self.dir)

    def todict(self):
        return {'name': 'FixedPlane',
                'kwargs': {'a': self.a, 'direction': self.dir}}
        
    def __repr__(self):
        return 'FixedPlane(%d, %s)' % (self.a, self.dir.tolist())


class FixedLine(FixConstraintSingle):
    """Constrain an atom index *a* to move on a given line only.

    The line is defined by its vector *direction*."""

    def __init__(self, a, direction):
        self.a = a
        self.dir = np.asarray(direction) / sqrt(np.dot(direction, direction))

    def adjust_positions(self, atoms, newpositions):
        step = newpositions[self.a] - atoms.positions[self.a]
        x = np.dot(step, self.dir)
        newpositions[self.a] = atoms.positions[self.a] + x * self.dir

    def adjust_forces(self, atoms, forces):
        forces[self.a] = self.dir * np.dot(forces[self.a], self.dir)

    def __repr__(self):
        return 'FixedLine(%d, %s)' % (self.a, self.dir.tolist())
        
    def todict(self):
        return {'name': 'FixedLine',
                'kwargs': {'a': self.a, 'direction': self.dir}}


class FixCartesian(FixConstraintSingle):
    'Fix an atom index *a* in the directions of the cartesian coordinates.'
    def __init__(self, a, mask=(1, 1, 1)):
        self.a = a
        self.mask = ~np.asarray(mask, bool)

    def adjust_positions(self, atoms, new):
        step = new[self.a] - atoms.positions[self.a]
        step *= self.mask
        new[self.a] = atoms.positions[self.a] + step

    def adjust_forces(self, atoms, forces):
        forces[self.a] *= self.mask

    def __repr__(self):
        return 'FixCartesian(a={0}, mask={1})'.format(self.a,
                                                      list(~self.mask))

    def todict(self):
        return {'name': 'FixCartesian',
                'kwargs': {'a': self.a, 'mask': ~self.mask}}

        
class FixScaled(FixConstraintSingle):
    'Fix an atom index *a* in the directions of the unit vectors.'
    def __init__(self, cell, a, mask=(1, 1, 1)):
        self.cell = np.asarray(cell)
        self.a = a
        self.mask = np.array(mask)

    def adjust_positions(self, atoms, new):
        scaled_old = np.linalg.solve(self.cell.T, atoms.positions.T).T
        scaled_new = np.linalg.solve(self.cell.T, new.T).T
        for n in range(3):
            if self.mask[n]:
                scaled_new[self.a, n] = scaled_old[self.a, n]
        new[self.a] = np.dot(scaled_new, self.cell)[self.a]

    def adjust_forces(self, atoms, forces):
        scaled_forces = np.linalg.solve(self.cell.T, forces.T).T
        scaled_forces[self.a] *= -(self.mask - 1)
        forces[self.a] = np.dot(scaled_forces, self.cell)[self.a]

    def todict(self):
        return {'name': 'FixScaled',
                'kwargs': {'a': self.a,
                           'cell': self.cell,
                           'mask': self.mask}}

    def __repr__(self):
        return 'FixScaled(%s, %d, %s)' % (repr(self.cell),
                                          self.a,
                                          repr(self.mask))


# TODO: Better interface might be to use dictionaries in place of very
# nested lists/tuples
class FixInternals(FixConstraint):
    """Constraint object for fixing multiple internal coordinates.

    Allows fixing bonds, angles, and dihedrals."""
    def __init__(self, bonds=None, angles=None, dihedrals=None,
                 epsilon=1.e-7):
        self.bonds = bonds or []
        self.angles = angles or []
        self.dihedrals = dihedrals or []
        
        # Initialize these at run-time:
        self.n = 0
        self.constraints = []
        self.epsilon = epsilon
        
        self.initialized = False
    
    def initialize(self, atoms):
        if self.initialized:
            return
        masses = atoms.get_masses()
        self.n = len(self.bonds) + len(self.angles) + len(self.dihedrals)
        self.constraints = []
        for bond in self.bonds:
            masses_bond = masses.take(bond[1])
            self.constraints.append(self.FixBondLengthAlt(bond[0], bond[1],
                                                          masses_bond))
        for angle in self.angles:
            masses_angle = masses.take(angle[1])
            self.constraints.append(self.FixAngle(angle[0], angle[1],
                                                  masses_angle))
        for dihedral in self.dihedrals:
            masses_dihedral = masses.take(dihedral[1])
            self.constraints.append(self.FixDihedral(dihedral[0],
                                                     dihedral[1],
                                                     masses_dihedral))
        self.initialized = True

    def todict(self):
        return {'name': 'FixInternals',
                'kwargs': {'bonds': self.bonds,
                           'angles': self.angles,
                           'dihedrals': self.dihedrals,
                           'epsilon': self.epsilon}}
        
    def adjust_positions(self, atoms, new):
        self.initialize(atoms)
        for constraint in self.constraints:
            constraint.set_h_vectors(atoms.positions)
        for j in range(50):
            maxerr = 0.0
            for constraint in self.constraints:
                constraint.adjust_positions(atoms.positions, new)
                maxerr = max(abs(constraint.sigma), maxerr)
            if maxerr < self.epsilon:
                return
        raise ValueError('Shake did not converge.')

    def adjust_forces(self, atoms, forces):
        """Project out translations and rotations and all other constraints"""
        self.initialize(atoms)
        positions = atoms.positions
        N = len(forces)
        list2_constraints = list(np.zeros((6, N, 3)))
        tx, ty, tz, rx, ry, rz = list2_constraints

        list_constraints = [r.ravel() for r in list2_constraints]
        
        tx[:, 0] = 1.0
        ty[:, 1] = 1.0
        tz[:, 2] = 1.0
        ff = forces.ravel()
        
        # Calculate the center of mass
        center = positions.sum(axis=0) / N

        rx[:, 1] = -(positions[:, 2] - center[2])
        rx[:, 2] = positions[:, 1] - center[1]
        ry[:, 0] = positions[:, 2] - center[2]
        ry[:, 2] = -(positions[:, 0] - center[0])
        rz[:, 0] = -(positions[:, 1] - center[1])
        rz[:, 1] = positions[:, 0] - center[0]
        
        # Normalizing transl., rotat. constraints
        for r in list2_constraints:
            r /= np.linalg.norm(r.ravel())
        
        # Add all angle, etc. constraint vectors
        for constraint in self.constraints:
            constraint.adjust_forces(positions, forces)
            list_constraints.insert(0, constraint.h)
        # QR DECOMPOSITION - GRAM SCHMIDT

        list_constraints = [r.ravel() for r in list_constraints]
        aa = np.column_stack(list_constraints)
        (aa, bb) = np.linalg.qr(aa)
        # Projection
        hh = []
        for i, constraint in enumerate(self.constraints):
            hh.append(aa[:, i] * np.row_stack(aa[:, i]))

        txx = aa[:, self.n] * np.row_stack(aa[:, self.n])
        tyy = aa[:, self.n + 1] * np.row_stack(aa[:, self.n + 1])
        tzz = aa[:, self.n + 2] * np.row_stack(aa[:, self.n + 2])
        rxx = aa[:, self.n + 3] * np.row_stack(aa[:, self.n + 3])
        ryy = aa[:, self.n + 4] * np.row_stack(aa[:, self.n + 4])
        rzz = aa[:, self.n + 5] * np.row_stack(aa[:, self.n + 5])
        T = txx + tyy + tzz + rxx + ryy + rzz
        for vec in hh:
            T += vec
        ff = np.dot(T, np.row_stack(ff))
        forces[:, :] -= np.dot(T, np.row_stack(ff)).reshape(-1, 3)

    def __repr__(self):
        constraints = repr(self.constraints)
        return 'FixInternals(_copy_init=%s, epsilon=%s)' % (constraints,
                                                            repr(self.epsilon))
    
    def __str__(self):
        return '\n'.join([repr(c) for c in self.constraints])

    # Classes for internal use in FixInternals
    class FixBondLengthAlt:
        """Constraint subobject for fixing bond length within FixInternals."""
        def __init__(self, bond, indices, masses, maxstep=0.01):
            """Fix distance between atoms with indices a1, a2."""
            self.indices = indices
            self.bond = bond
            self.h1 = None
            self.h2 = None
            self.masses = masses
            self.h = []
            self.sigma = 1.

        def set_h_vectors(self, pos):
            dist1 = pos[self.indices[0]] - pos[self.indices[1]]
            self.h1 = 2 * dist1
            self.h2 = -self.h1

        def adjust_positions(self, old, new):
            h1 = self.h1 / self.masses[0]
            h2 = self.h2 / self.masses[1]
            dist1 = new[self.indices[0]] - new[self.indices[1]]
            dist = np.dot(dist1, dist1)
            self.sigma = dist - self.bond**2
            lamda = -self.sigma / (2 * np.dot(dist1, (h1 - h2)))
            new[self.indices[0]] += lamda * h1
            new[self.indices[1]] += lamda * h2

        def adjust_forces(self, positions, forces):
            self.h1 = 2 * (positions[self.indices[0]] -
                           positions[self.indices[1]])
            self.h2 = -self.h1
            self.h = np.zeros([len(forces) * 3])
            self.h[(self.indices[0]) * 3] = self.h1[0]
            self.h[(self.indices[0]) * 3 + 1] = self.h1[1]
            self.h[(self.indices[0]) * 3 + 2] = self.h1[2]
            self.h[(self.indices[1]) * 3] = self.h2[0]
            self.h[(self.indices[1]) * 3 + 1] = self.h2[1]
            self.h[(self.indices[1]) * 3 + 2] = self.h2[2]
            self.h /= np.linalg.norm(self.h)

        def __repr__(self):
            return 'FixBondLengthAlt(%s, %d, %d)' % \
                (repr(self.bond), self.indices[0], self.indices[1])

    class FixAngle:
        """Constraint object for fixing an angle within
        FixInternals."""
        def __init__(self, angle, indices, masses):
            """Fix atom movement to construct a constant angle."""
            self.indices = indices
            self.a1m, self.a2m, self.a3m = masses
            self.angle = np.cos(angle)
            self.h1 = self.h2 = self.h3 = None
            self.h = []
            self.sigma = 1.

        def set_h_vectors(self, pos):
            r21 = pos[self.indices[0]] - pos[self.indices[1]]
            r21_len = np.linalg.norm(r21)
            e21 = r21 / r21_len
            r23 = pos[self.indices[2]] - pos[self.indices[1]]
            r23_len = np.linalg.norm(r23)
            e23 = r23 / r23_len
            angle = np.dot(e21, e23)
            self.h1 = -2 * angle * ((angle * e21 - e23) / (r21_len))
            self.h3 = -2 * angle * ((angle * e23 - e21) / (r23_len))
            self.h2 = -(self.h1 + self.h3)

        def adjust_positions(self, oldpositions, newpositions):
            r21 = newpositions[self.indices[0]] - newpositions[self.indices[1]]
            r21_len = np.linalg.norm(r21)
            e21 = r21 / r21_len
            r23 = newpositions[self.indices[2]] - newpositions[self.indices[1]]
            r23_len = np.linalg.norm(r23)
            e23 = r23 / r23_len
            angle = np.dot(e21, e23)
            self.sigma = (angle - self.angle) * (angle + self.angle)
            h1 = self.h1 / self.a1m
            h3 = self.h3 / self.a3m
            h2 = self.h2 / self.a2m
            h21 = h1 - h2
            h23 = h3 - h2
            # Calculating new positions
            deriv = (((np.dot(r21, h23) + np.dot(r23, h21))
                      / (r21_len * r23_len))
                     - (np.dot(r21, h21) / (r21_len * r21_len)
                        + np.dot(r23, h23) / (r23_len * r23_len)) * angle)
            deriv *= 2 * angle
            lamda = -self.sigma / deriv
            newpositions[self.indices[0]] += lamda * h1
            newpositions[self.indices[1]] += lamda * h2
            newpositions[self.indices[2]] += lamda * h3

        def adjust_forces(self, positions, forces):
            r21 = positions[self.indices[0]] - positions[self.indices[1]]
            r21_len = np.linalg.norm(r21)
            e21 = r21 / r21_len
            r23 = positions[self.indices[2]] - positions[self.indices[1]]
            r23_len = np.linalg.norm(r23)
            e23 = r23 / r23_len
            angle = np.dot(e21, e23)
            self.h1 = -2 * angle * (angle * e21 - e23) / r21_len
            self.h3 = -2 * angle * (angle * e23 - e21) / r23_len
            self.h2 = -(self.h1 + self.h3)
            self.h = np.zeros([len(positions) * 3])
            self.h[(self.indices[0]) * 3] = self.h1[0]
            self.h[(self.indices[0]) * 3 + 1] = self.h1[1]
            self.h[(self.indices[0]) * 3 + 2] = self.h1[2]
            self.h[(self.indices[1]) * 3] = self.h2[0]
            self.h[(self.indices[1]) * 3 + 1] = self.h2[1]
            self.h[(self.indices[1]) * 3 + 2] = self.h2[2]
            self.h[(self.indices[2]) * 3] = self.h3[0]
            self.h[(self.indices[2]) * 3 + 1] = self.h3[1]
            self.h[(self.indices[2]) * 3 + 2] = self.h3[2]
            self.h /= np.linalg.norm(self.h)

        def __repr__(self):
            return 'FixAngle(%s, %f)' % (tuple(self.indices),
                                         np.arccos(self.angle))
 
    class FixDihedral:
        """Constraint object for fixing an dihedral using
        the shake algorithm. This one allows also other constraints."""
        def __init__(self, angle, indices, masses):
            """Fix atom movement to construct a constant dihedral angle."""
            self.indices = indices
            self.a1m, self.a2m, self.a3m, self.a4m = masses
            self.angle = np.cos(angle)
            self.h1 = self.h2 = self.h3 = self.h4 = None
            self.h = []
            self.sigma = 1.

        def set_h_vectors(self, pos):
            r12 = pos[self.indices[1]] - pos[self.indices[0]]
            r23 = pos[self.indices[2]] - pos[self.indices[1]]
            r23_len = np.linalg.norm(r23)
            e23 = r23 / r23_len
            r34 = pos[self.indices[3]] - pos[self.indices[2]]
            a = -r12 - np.dot(-r12, e23) * e23
            a_len = np.linalg.norm(a)
            ea = a / a_len
            b = r34 - np.dot(r34, e23) * e23
            b_len = np.linalg.norm(b)
            eb = b / b_len
            angle = np.dot(ea, eb).clip(-1.0, 1.0)
            self.h1 = (eb - angle * ea) / a_len
            self.h4 = (ea - angle * eb) / b_len
            self.h2 = self.h1 * (np.dot(-r12, e23) / r23_len - 1)
            self.h2 += np.dot(r34, e23) / r23_len * self.h4
            self.h3 = -self.h4 * (np.dot(r34, e23) / r23_len + 1)
            self.h3 += np.dot(r12, e23) / r23_len * self.h1

        def adjust_positions(self, oldpositions, newpositions):
            r12 = newpositions[self.indices[1]] - newpositions[self.indices[0]]
            r23 = newpositions[self.indices[2]] - newpositions[self.indices[1]]
            r34 = newpositions[self.indices[3]] - newpositions[self.indices[2]]
            n1 = np.cross(r12, r23)
            n1_len = np.linalg.norm(n1)
            n1e = n1 / n1_len
            n2 = np.cross(r23, r34)
            n2_len = np.linalg.norm(n2)
            n2e = n2 / n2_len
            angle = np.dot(n1e, n2e).clip(-1.0, 1.0)
            self.sigma = (angle - self.angle) * (angle + self.angle)
            h1 = self.h1 / self.a1m
            h2 = self.h2 / self.a2m
            h3 = self.h3 / self.a3m
            h4 = self.h4 / self.a4m
            h12 = h2 - h1
            h23 = h3 - h2
            h34 = h4 - h3
            deriv = ((np.dot(n1, np.cross(r34, h23) + np.cross(h34, r23))
                      + np.dot(n2, np.cross(r23, h12) + np.cross(h23, r12)))
                     / (n1_len * n2_len))
            deriv -= (((np.dot(n1, np.cross(r23, h12) + np.cross(h23, r12))
                        / n1_len**2)
                       + (np.dot(n2, np.cross(r34, h23) + np.cross(h34, r23))
                          / n2_len**2)) * angle)
            deriv *= -2 * angle
            lamda = -self.sigma / deriv
            newpositions[self.indices[0]] += lamda * h1
            newpositions[self.indices[1]] += lamda * h2
            newpositions[self.indices[2]] += lamda * h3
            newpositions[self.indices[3]] += lamda * h4

        def adjust_forces(self, positions, forces):
            r12 = positions[self.indices[1]] - positions[self.indices[0]]
            r23 = positions[self.indices[2]] - positions[self.indices[1]]
            r23_len = np.linalg.norm(r23)
            e23 = r23 / r23_len
            r34 = positions[self.indices[3]] - positions[self.indices[2]]
            a = -r12 - np.dot(-r12, e23) * e23
            a_len = np.linalg.norm(a)
            ea = a / a_len
            b = r34 - np.dot(r34, e23) * e23
            b_len = np.linalg.norm(b)
            eb = b / b_len
            angle = np.dot(ea, eb).clip(-1.0, 1.0)
            self.h1 = (eb - angle * ea) / a_len
            self.h4 = (ea - angle * eb) / b_len
            self.h2 = self.h1 * (np.dot(-r12, e23) / r23_len - 1)
            self.h2 += np.dot(r34, e23) / r23_len * self.h4
            self.h3 = -self.h4 * (np.dot(r34, e23) / r23_len + 1)
            self.h3 -= np.dot(-r12, e23) / r23_len * self.h1

            self.h = np.zeros([len(positions) * 3])
            self.h[(self.indices[0]) * 3] = self.h1[0]
            self.h[(self.indices[0]) * 3 + 1] = self.h1[1]
            self.h[(self.indices[0]) * 3 + 2] = self.h1[2]
            self.h[(self.indices[1]) * 3] = self.h2[0]
            self.h[(self.indices[1]) * 3 + 1] = self.h2[1]
            self.h[(self.indices[1]) * 3 + 2] = self.h2[2]
            self.h[(self.indices[2]) * 3] = self.h3[0]
            self.h[(self.indices[2]) * 3 + 1] = self.h3[1]
            self.h[(self.indices[2]) * 3 + 2] = self.h3[2]
            self.h[(self.indices[3]) * 3] = self.h4[0]
            self.h[(self.indices[3]) * 3 + 1] = self.h4[1]
            self.h[(self.indices[3]) * 3 + 2] = self.h4[2]
            self.h /= np.linalg.norm(self.h)
        
        def __repr__(self):
            return 'FixDihedral(%s, %f)' % (tuple(self.indices), self.angle)


class Hookean(FixConstraint):
    """Applies a Hookean restorative force between a pair of atoms, an atom
    and a point, or an atom and a plane."""

    def __init__(self, a1, a2, k, rt=None):
        """Forces two atoms to stay close together by applying no force if
        they are below a threshold length, rt, and applying a Hookean
        restorative force when the distance between them exceeds rt. Can
        also be used to tether an atom to a fixed point in space or to a
        distance above a plane.

        a1 : int
           Index of atom 1
        a2 : one of three options
           1) index of atom 2
           2) a fixed point in cartesian space to which to tether a1
           3) a plane given as (A, B, C, D) in A x + B y + C z + D = 0.
        k : float
           Hooke's law (spring) constant to apply when distance
           exceeds threshold_length. Units of eV A^-2.
        rt : float
           The threshold length below which there is no force. The
           length is 1) between two atoms, 2) between atom and point.
           This argument is not supplied in case 3. Units of A.

        If a plane is specified, the Hooke's law force is applied if the atom
        is on the normal side of the plane. For instance, the plane with
        (A, B, C, D) = (0, 0, 1, -7) defines a plane in the xy plane with a z
        intercept of +7 and a normal vector pointing in the +z direction.
        If the atom has z > 7, then a downward force would be applied of
        k * (atom.z - 7). The same plane with the normal vector pointing in
        the -z direction would be given by (A, B, C, D) = (0, 0, -1, 7).
        """

        if isinstance(a2, int):
            self._type = 'two atoms'
            self.indices = [a1, a2]
        elif len(a2) == 3:
            self._type = 'point'
            self.index = a1
            self.origin = np.array(a2)
        elif len(a2) == 4:
            self._type = 'plane'
            self.index = a1
            self.plane = a2
        else:
            raise RuntimeError('Unknown type for a2')
        self.threshold = rt
        self.spring = k

    def todict(self):
        dct = {'name': 'Hookean'}
        dct['kwargs'] = {'rt': self.threshold,
                         'k': self.spring}
        if self._type == 'two atoms':
            dct['kwargs']['a1'] = self.indices[0]
            dct['kwargs']['a2'] = self.indices[1]
        elif self._type == 'point':
            dct['kwargs']['a1'] = self.index
            dct['kwargs']['a2'] = self.origin
        elif self._type == 'plane':
            dct['kwargs']['a1'] = self.index
            dct['kwargs']['a2'] = self.plane
        else:
            raise NotImplementedError('Bad type: %s' % self._type)
        return dct

    def adjust_positions(self, atoms, newpositions):
        pass

    def adjust_momenta(self, atoms, momenta):
        pass

    def adjust_forces(self, atoms, forces):
        positions = atoms.positions
        if self._type == 'plane':
            A, B, C, D = self.plane
            x, y, z = positions[self.index]
            d = ((A * x + B * y + C * z + D) /
                 np.sqrt(A**2 + B**2 + C**2))
            if d < 0:
                return
            magnitude = self.spring * d
            direction = - np.array((A, B, C)) / np.linalg.norm((A, B, C))
            forces[self.index] += direction * magnitude
            return
        if self._type == 'two atoms':
            p1, p2 = positions[self.indices]
        elif self._type == 'point':
            p1 = positions[self.index]
            p2 = self.origin
        displace = p2 - p1
        bondlength = np.linalg.norm(displace)
        if bondlength > self.threshold:
            magnitude = self.spring * (bondlength - self.threshold)
            direction = displace / np.linalg.norm(displace)
            if self._type == 'two atoms':
                forces[self.indices[0]] += direction * magnitude
                forces[self.indices[1]] -= direction * magnitude
            else:
                forces[self.index] += direction * magnitude

    def adjust_potential_energy(self, atoms):
        """Returns the difference to the potential energy due to an active
        constraint. (That is, the quantity returned is to be added to the
        potential energy.)"""
        positions = atoms.positions
        if self._type == 'plane':
            A, B, C, D = self.plane
            x, y, z = positions[self.index]
            d = ((A * x + B * y + C * z + D) /
                 np.sqrt(A**2 + B**2 + C**2))
            if d > 0:
                return 0.5 * self.spring * d**2
            else:
                return 0.
        if self._type == 'two atoms':
            p1, p2 = positions[self.indices]
        elif self._type == 'point':
            p1 = positions[self.index]
            p2 = self.origin
        displace = p2 - p1
        bondlength = np.linalg.norm(displace)
        if bondlength > self.threshold:
            return 0.5 * self.spring * (bondlength - self.threshold)**2
        else:
            return 0.

    def index_shuffle(self, atoms, ind):
        # See docstring of superclass
        if self._type == 'two atoms':
            self.indices = [ind.index(self.indices[0]),
                            ind.index(self.indices[1])]
        elif self._type == 'point':
            self.index = ind.index(self.index)
        elif self._type == 'plane':
            self.index = ind.index(self.index)

    def __repr__(self):
        if self._type == 'two atoms':
            return 'Hookean(%d, %d)' % tuple(self.indices)
        elif self._type == 'point':
            return 'Hookean(%d) to cartesian' % self.index
        else:
            return 'Hookean(%d) to plane' % self.index

class stretchcombo:
    """Constrain a linear combination of bond-length . """
    """bondlenght[i] =[[ atom_i_i, atom_i_j, weight-factor_i]]"""
    def __init__(self, a, bondlist, atoms,invert=False):
        xyz = atoms.get_positions()
        self.invert=invert
        self.dim=len(xyz)
        self.xyz=np.copy(xyz)
        self.cdim=len(bondlist)
        self.thr=1.0e-12
        self.bondlist=list(bondlist)
        self.projected_force = 0.0
        self.f_thresh  = 0.0
        self.full_force = 0.0
        self.update_constraint(a,xyz)

    def index_shuffle(self, atoms, ind):
        """Shuffle the indices of the two atoms in this constraint"""

        for new, old in slice2enlist(ind.tolist(), len(atoms)):
            for i, a in enumerate(self.bondlist):
                for j in range(2):
                    if old == a[j]:
                        print("Reposition:", old, new)
                        self.bondlist[i][j] == new

    def constraint_to_dimermethod(self):
        displacement_vector = [[0.0,0.0,0.0] for ii in range(self.dim)]
        maskding=[0 for i in range(self.dim)]
        for schieber in self.bondlist:
          diff_vec = schieber[2]* ( self.xyz[schieber[1]] - self.xyz[schieber[0]] )
          ii=0
          for at in [schieber[0],schieber[1]]:
            displacement_vector[at] = (-1)**ii * diff_vec
            maskding[at] = 1
            ii+=1
        displacement_vector /= np.linalg.norm(displacement_vector)
        return maskding,displacement_vector


    def update_constraint(self,a,xyz):
        self.a = a
        self.val0=0.0
        self.update(xyz)
        try:
           self.val0=float(a)
           print('set constrained coordinate to',self.val0,'from current value',self.val)
           oldxyz=np.copy(xyz)
           self.adjust_positions(oldxyz,xyz)
        except:
           self.val0=self.val
           self.a = self.val0
           print('fix constrained coordinate to initial value',self.val0)

    def reset_a(self, a,atoms):
        self.update_constraint(a,atoms.get_positions())
        return self.xyz

    def update(self, xyz):
        self.val=0.0
        self.xyz = xyz
        self.dir=np.zeros([self.dim,3])
        self.bond=np.zeros(self.cdim)
        self.involved_atoms=[]
        ii=-1
        for i,j,fac in self.bondlist:
          ii+=1
          if i not in self.involved_atoms:
            self.involved_atoms += [i]
          if j not in self.involved_atoms:
            self.involved_atoms += [j]
          self.bond[ii]=np.linalg.norm(xyz[i]-xyz[j])
          self.val+=fac*self.bond[ii]
          self.dir[i] += fac*(xyz[i]-xyz[j]) / self.bond[ii]
          self.dir[j] -= fac*(xyz[i]-xyz[j]) / self.bond[ii]
        self.der = np.linalg.norm(self.dir)
        self.maxdir=max([np.linalg.norm(ddd) for ddd in self.dir])
        self.dir = self.dir / self.der

    def return_adjusted_positions(self):
        xyz=self.xyz
        oldxyz=np.copy(xyz)
        self.adjust_positions(oldxyz,xyz)
        return self.xyz

    def adjust_positions(self, atoms       , newpositions):
        if self.invert:
          return
        try:
          oldpositions = atoms.get_positions()
        except:
          oldpositions = atoms
        self.update(oldpositions)
        step=np.zeros([self.dim,3])
        for i in range(0,self.dim):
          step[i] = newpositions[i] - oldpositions[i]
        proj=0.0
        for i in range(0,self.dim):
          proj += np.dot(step[i], self.dir[i])
        for i in range(0,self.dim):
          newpositions[i] = newpositions[i] - self.dir[i] * proj
        self.update(newpositions)

        i_geo=0
        imax_geo=100
        while (abs(self.val-self.val0)>self.thr): 
          i_geo+=1
          if i_geo > imax_geo:
            print('failed to prepare appropriate structure')
            sys.exit(1)
          correction_direction = self.dir.copy()
          tonorm = np.linalg.norm(correction_direction) 

          for i in range(0,self.dim):
            newpositions[i] = newpositions[i] - (self.val-self.val0) * correction_direction[i] / ( tonorm * self.der )
          self.update(newpositions)

        if i_geo != 0:
          print('... adjusted constraint to new val = ',self.val,'required ',i_geo,'iterations')
          # is called in practise twice --> apears in output.
          # reason: when, in scan_constraints, atoms.set_positions() is called, this calls again
          # adjust position. no harm..

          
        if(abs(self.val-self.val0)>self.thr): 
          print('primitive first-come-first-served correction from',self.val,'to',self.val0)

          for i in range(0,self.dim):
            newpositions[i] = newpositions[i] - self.dir[i] *(self.val-self.val0) * self.der
            self.update(newpositions)
          print('... done: new val',self.val)


    def adjust_forces(self,atoms,  forces):
        positions = atoms.get_positions()
        self.update(positions)
        proj=0.0
        self.forces_weighted_dir = []
        for i in range(0,len(forces)):
          proj += np.dot(forces[i], self.dir[i])
          self.forces_weighted_dir += [ self.dir[i] * proj *  np.dot(forces[i], self.dir[i]) ]
        self.projected_force = proj / self.der
        self.f_thresh  = self.projected_force * self.der**2 / (self.maxdir * float(len(self.involved_atoms)))
        self.full_force = max([np.linalg.norm(ff) for ff in forces])
        self.projected_forces=[]
        if self.invert:
          proj *= 2.0
        for i in range(0,len(forces)):
          forces[i] = forces[i] - self.dir[i] * proj
          self.projected_forces += [ self.dir[i] * proj ]

class Filter:
    """Subset filter class."""
    def __init__(self, atoms, indices=None, mask=None):
        """Filter atoms.

        This filter can be used to hide degrees of freedom in an Atoms
        object.

        Parameters
        ----------
        indices : list of int
           Indices for those atoms that should remain visible.
        mask : list of bool
           One boolean per atom indicating if the atom should remain
           visible or not.

        If a Trajectory tries to save this object, it will instead
        save the underlying Atoms object.  To prevent this, delete
        the atoms_for_saving attribute.
        """

        self.atoms = atoms
        self.constraints = []
        # Make self.info a reference to the underlying atoms' info dictionary.
        self.info = self.atoms.info

        if indices is None and mask is None:
            raise ValueError('Use "indices" or "mask".')
        if indices is not None and mask is not None:
            raise ValueError('Use only one of "indices" and "mask".')

        if mask is not None:
            self.index = np.asarray(mask, bool)
            self.n = self.index.sum()
        else:
            self.index = np.asarray(indices, int)
            self.n = len(self.index)

        # Present the real atoms object to Trajectory and friends
        self.atoms_for_saving = self.atoms

    def get_cell(self):
        """Returns the computational cell.

        The computational cell is the same as for the original system.
        """
        return self.atoms.get_cell()

    def get_pbc(self):
        """Returns the periodic boundary conditions.

        The boundary conditions are the same as for the original system.
        """
        return self.atoms.get_pbc()

    def get_positions(self):
        'Return the positions of the visible atoms.'
        return self.atoms.get_positions()[self.index]

    def set_positions(self, positions):
        'Set the positions of the visible atoms.'
        pos = self.atoms.get_positions()
        pos[self.index] = positions
        self.atoms.set_positions(pos)

    positions = property(get_positions, set_positions,
                         doc='Positions of the atoms')

    def get_momenta(self):
        'Return the momenta of the visible atoms.'
        return self.atoms.get_momenta()[self.index]

    def set_momenta(self, momenta, **kwargs):
        'Set the momenta of the visible atoms.'
        mom = self.atoms.get_momenta()
        mom[self.index] = momenta
        self.atoms.set_momenta(mom, **kwargs)

    def get_atomic_numbers(self):
        'Return the atomic numbers of the visible atoms.'
        return self.atoms.get_atomic_numbers()[self.index]

    def set_atomic_numbers(self, atomic_numbers):
        'Set the atomic numbers of the visible atoms.'
        z = self.atoms.get_atomic_numbers()
        z[self.index] = atomic_numbers
        self.atoms.set_atomic_numbers(z)

    def get_tags(self):
        'Return the tags of the visible atoms.'
        return self.atoms.get_tags()[self.index]

    def set_tags(self, tags):
        'Set the tags of the visible atoms.'
        tg = self.atoms.get_tags()
        tg[self.index] = tags
        self.atoms.set_tags(tg)

    def get_forces(self, *args, **kwargs):
        return self.atoms.get_forces(*args, **kwargs)[self.index]

    def get_stress(self):
        return self.atoms.get_stress()

    def get_stresses(self):
        return self.atoms.get_stresses()[self.index]

    def get_masses(self):
        return self.atoms.get_masses()[self.index]

    def get_potential_energy(self, **kwargs):
        """Calculate potential energy.

        Returns the potential energy of the full system.
        """
        return self.atoms.get_potential_energy(**kwargs)

    def get_chemical_symbols(self):
        return self.atoms.get_chemical_symbols()

    def get_initial_magnetic_moments(self):
        return self.atoms.get_initial_magnetic_moments()

    def get_calculator(self):
        """Returns the calculator.

        WARNING: The calculator is unaware of this filter, and sees a
        different number of atoms.
        """
        return self.atoms.get_calculator()
    
    def get_celldisp(self):
        return self.atoms.get_celldisp()

    def has(self, name):
        'Check for existence of array.'
        return self.atoms.has(name)

    def __len__(self):
        'Return the number of movable atoms.'
        return self.n

    def __getitem__(self, i):
        'Return an atom.'
        return self.atoms[self.index[i]]
    

class StrainFilter(Filter):
    """Modify the supercell while keeping the scaled positions fixed.

    Presents the strain of the supercell as the generalized positions,
    and the global stress tensor (times the volume) as the generalized
    force.

    This filter can be used to relax the unit cell until the stress is
    zero.  If MDMin is used for this, the timestep (dt) to be used
    depends on the system size. 0.01/x where x is a typical dimension
    seems like a good choice.

    The stress and strain are presented as 6-vectors, the order of the
    components follow the standard engingeering practice: xx, yy, zz,
    yz, xz, xy.

    """
    def __init__(self, atoms, mask=None):
        """Create a filter applying a homogeneous strain to a list of atoms.

        The first argument, atoms, is the atoms object.

        The optional second argument, mask, is a list of six booleans,
        indicating which of the six independent components of the
        strain that are allowed to become non-zero.  It defaults to
        [1,1,1,1,1,1].

        """

        self.atoms = atoms
        self.strain = np.zeros(6)

        if mask is None:
            self.mask = np.ones(6)
        else:
            self.mask = np.array(mask)

        self.index = np.arange(len(atoms))
        self.n = self.index.sum()

        self.origcell = atoms.get_cell()

    def get_positions(self):
        return self.strain.reshape((2, 3)).copy()

    def set_positions(self, new):
        new = new.ravel() * self.mask
        eps = np.array([[1.0 + new[0], 0.5 * new[5], 0.5 * new[4]],
                        [0.5 * new[5], 1.0 + new[1], 0.5 * new[3]],
                        [0.5 * new[4], 0.5 * new[3], 1.0 + new[2]]])

        self.atoms.set_cell(np.dot(self.origcell, eps), scale_atoms=True)
        self.strain[:] = new

    def get_forces(self):
        stress = self.atoms.get_stress()
        return -self.atoms.get_volume() * (stress * self.mask).reshape((2, 3))

    def get_potential_energy(self):
        return self.atoms.get_potential_energy()

    def has(self, x):
        return self.atoms.has(x)

    def __len__(self):
        return 2


class UnitCellFilter(Filter):
    """Modify the supercell and the atom positions. """
    def __init__(self, atoms, mask=None):
        """Create a filter that returns the atomic forces and unit cell
        stresses together, so they can simultaneously be minimized.

        The first argument, atoms, is the atoms object. The optional second
        argument, mask, is a list of booleans, indicating which of the six
        independent components of the strain are relaxed.

        - True = relax to zero
        - False = fixed, ignore this component

        You can still use constraints on the atoms, e.g. FixAtoms, to control
        the relaxation of the atoms.

        >>> # this should be equivalent to the StrainFilter
        >>> atoms = Atoms(...)
        >>> atoms.set_constraint(FixAtoms(mask=[True for atom in atoms]))
        >>> ucf = UnitCellFilter(atoms)

        You should not attach this UnitCellFilter object to a
        trajectory. Instead, create a trajectory for the atoms, and
        attach it to an optimizer like this:

        >>> atoms = Atoms(...)
        >>> ucf = UnitCellFilter(atoms)
        >>> qn = QuasiNewton(ucf)
        >>> traj = Trajectory('TiO2.traj', 'w', atoms)
        >>> qn.attach(traj)
        >>> qn.run(fmax=0.05)

        Helpful conversion table:

        - 0.05 eV/A^3   = 8 GPA
        - 0.003 eV/A^3  = 0.48 GPa
        - 0.0006 eV/A^3 = 0.096 GPa
        - 0.0003 eV/A^3 = 0.048 GPa
        - 0.0001 eV/A^3 = 0.02 GPa
        """

        Filter.__init__(self, atoms, indices=list(range(len(atoms))))

        self.atoms = atoms
        self.strain = np.zeros(6)

        if mask is None:
            self.mask = np.ones(6)
        else:
            self.mask = np.array(mask)

        self.origcell = atoms.get_cell()
        self.copy = self.atoms.copy
        self.arrays = self.atoms.arrays

    def get_positions(self):
        '''
        this returns an array with shape (natoms + 2,3).

        the first natoms rows are the positions of the atoms, the last
        two rows are the strains associated with the unit cell
        '''

        atom_positions = self.atoms.get_positions()
        strains = self.strain.reshape((2, 3))

        natoms = len(self.atoms)
        all_pos = np.zeros((natoms + 2, 3), np.float)
        all_pos[0:natoms, :] = atom_positions
        all_pos[natoms:, :] = strains

        return all_pos

    def set_positions(self, new):
        '''
        new is an array with shape (natoms+2,3).

        the first natoms rows are the positions of the atoms, the last
        two rows are the strains used to change the cell shape.

        The atom positions are set first, then the unit cell is
        changed keeping the atoms in their scaled positions.
        '''

        natoms = len(self.atoms)

        atom_positions = new[0:natoms, :]
        self.atoms.set_positions(atom_positions)

        new = new[natoms:, :]  # this is only the strains
        new = new.ravel() * self.mask
        eps = np.array([[1.0 + new[0], 0.5 * new[5], 0.5 * new[4]],
                        [0.5 * new[5], 1.0 + new[1], 0.5 * new[3]],
                        [0.5 * new[4], 0.5 * new[3], 1.0 + new[2]]])

        self.atoms.set_cell(np.dot(self.origcell, eps), scale_atoms=True)
        self.strain[:] = new

    def get_forces(self, apply_constraint=False):
        '''
        returns an array with shape (natoms+2,3) of the atomic forces
        and unit cell stresses.

        the first natoms rows are the forces on the atoms, the last
        two rows are the stresses on the unit cell, which have been
        reshaped to look like "atomic forces". i.e.,

        f[-2] = -vol*[sxx,syy,szz]*mask[0:3]
        f[-1] = -vol*[syz, sxz, sxy]*mask[3:]

        apply_constraint is an argument expected by ase
        '''

        stress = self.atoms.get_stress()
        atom_forces = self.atoms.get_forces()

        natoms = len(self.atoms)
        all_forces = np.zeros((natoms + 2, 3), np.float)
        all_forces[0:natoms, :] = atom_forces

        vol = self.atoms.get_volume()
        stress_forces = -vol * (stress * self.mask).reshape((2, 3))
        all_forces[natoms:, :] = stress_forces
        return all_forces

    def get_potential_energy(self):
        return self.atoms.get_potential_energy()

    def has(self, x):
        return self.atoms.has(x)

    def __len__(self):
        return (2 + len(self.atoms))

