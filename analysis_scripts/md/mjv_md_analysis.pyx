from __future__ import division
import MDAnalysis
import numpy as np
cimport numpy as np


class Atom:
    def __init__(self, x, y, z, mass):
        self.x = x
        self.y = y
        self.z = z
        self.mass = mass

    @classmethod
    def from_mdanalysis(cls, positions, mass):
        return cls(positions[0], positions[1], positions[2], mass)

    def translate(self, arr):
        self.x += arr[0]
        self.y += arr[1]
        self.z += arr[2]

    def display(self):
        print "[%.8f, %.8f, %.8f], %.8f" % (self.x, self.y, self.z, self.mass)


class AtomGroup:
    def __init__(self, atoms, box=None):
        self.atoms = atoms
        self.box = box if box is not None else [0, 0, 0]

    @classmethod
    def create_from_mdanalysis(cls, atomGroup, numAtoms):
        tmp = cls([], atomGroup.dimensions[:3])
        for i in range(0, len(atomGroup.atoms)):
            tmp.add_atom(Atom.from_mdanalysis(atomGroup.atoms.positions[i], atomGroup.atoms.masses[i]))
        return tmp

    def center_of_mass(self):
        totalMass = sum(atom.mass for atom in self.atoms)
        com = np.array([0.0,0.0,0.0])
        for atom in self.atoms:
            com[0] += atom.mass * atom.x
            com[1] += atom.mass * atom.y
            com[2] += atom.mass * atom.z
        return com/totalMass

    def add_atom(self, atom):
        self.atoms.append(atom)


# def wrapMD(pro1, pro2, system):
#     if (len(pro1.atoms) != len(pro2.atoms)):
#         print "Atom count mismatch, exiting."
#         exit(1)
#
#     numAtoms = len(pro1.atoms)
#     x, y, z = dimensions[:3]
#
#     # Determine if any pairwise distance is not using the nearest image
#     wrapX, wrapY, wrapZ = False, False, False
#     for i in range(0, numAtoms):
#         for j in range(0, numAtoms):
#             if (np.abs(pro1.atoms.positions[i][0] - pro2.atoms.positions[j][0]) > x):
#                 wrapX = True
#             if (np.abs(pro1.atoms.positions[i][1] - pro2.atoms.positions[j][1]) > y):
#                 wrapY = True
#             if (np.abs(pro1.atoms.positions[i][2] - pro2.atoms.positions[j][2]) > z):
#                 wrapZ = True
#
#             if (wrapX == wrapY == wrapZ == True):
#                 break
#         else:
#             break
#
#     # Wrap coordinates
#     if wrapX or wrapY or wrapZ:
#         for i in range(0, numAtoms):
#             if wrapX:
#                 (pro1.atoms.positions[i][0] += x/2) if (pro1.atoms.positions[i][0] < x/2)
#                 (pro2.atoms.positions[i][0] += x/2) if (pro2.atoms.positions[i][0] < x/2)
#             if wrapY:
#                 (pro1.atoms.positions[i][1] += y/2) if (pro1.atoms.positions[i][1] < y/2)
#                 (pro2.atoms.positions[i][1] += y/2) if (pro2.atoms.positions[i][1] < y/2)
#             if wrapZ:
#                 (pro1.atoms.positions[i][2] += z/2) if (pro1.atoms.positions[i][2] < z/2)
#                 (pro2.atoms.positions[i][2] += z/2) if (pro2.atoms.positions[i][2] < z/2)
#
#     return
#
#

cpdef wrap(pro1, pro2):
    # if (len(pro1.atoms) != len(pro2.atoms)):
    #     print "Atom count mismatch, exiting."
    #     exit(1)
    if ((pro1.box != pro2.box).any() or (pro1.box == 0).any() or (pro2.box == 0).any()):
        print "Box dimension issue, exiting."
        exit(1)

    cpdef int numAtoms = len(pro1.atoms)
    # cpdef double x, y, z = pro1.box[0], pro1.box[1], pro1.box[2]
    cpdef float x = pro1.box[0]
    cpdef float y = pro1.box[1]
    cpdef float z = pro1.box[2]

    # Determine if any pairwise distance is not using the nearest image
    wrapX, wrapY, wrapZ = False, False, False
    for i in range(0, numAtoms):
        for j in range(0, numAtoms):
            if (np.abs(pro1.atoms[i].x - pro2.atoms[j].x) > x/2):
                wrapX = True
            if (np.abs(pro1.atoms[i].y - pro2.atoms[j].y) > y/2):
                wrapY = True
            if (np.abs(pro1.atoms[i].z - pro2.atoms[j].z) > z/2):
                wrapZ = True

            if (wrapX == wrapY == wrapZ == True):
                i = j = numAtoms

    print wrapX, wrapY, wrapZ

    # Wrap coordinates (inefficient right now)
    if wrapX or wrapY or wrapZ:
        for i in range(0, numAtoms):
            if wrapX:
                if (pro1.atoms[i].x < x/2):
                    pro1.atoms[i].x += x
                if (pro2.atoms[i].x < x/2):
                    pro2.atoms[i].x += x
            if wrapY:
                if (pro1.atoms[i].y < y/2):
                    pro1.atoms[i].y += y
                if (pro2.atoms[i].y < y/2):
                    pro2.atoms[i].y += y
            if wrapZ:
                if (pro1.atoms[i].z < z/2):
                    pro1.atoms[i].z += z
                if (pro2.atoms[i].z < z/2):
                    pro2.atoms[i].z += z

    return


cpdef wrap_opt(pro1, pro2):
    # if (len(pro1.atoms) != len(pro2.atoms)):
    #     print "Atom count mismatch, exiting."
    #     exit(1)
    if ((pro1.box != pro2.box).any() or (pro1.box == 0).any() or (pro2.box == 0).any()):
        print "Box dimension issue, exiting."
        exit(1)

    cpdef int numAtoms = len(pro1.atoms)
    # cpdef double x, y, z = pro1.box[0], pro1.box[1], pro1.box[2]
    cpdef float x = pro1.box[0]
    cpdef float y = pro1.box[1]
    cpdef float z = pro1.box[2]

    cdef np.ndarray[float, ndim=1] pro1X = np.array([atom.x for atom in pro1.atoms])
    cdef np.ndarray[float, ndim=1] pro1Y = np.array([atom.y for atom in pro1.atoms])
    cdef np.ndarray[float, ndim=1] pro1Z = np.array([atom.z for atom in pro1.atoms])

    cdef np.ndarray[float, ndim=1] pro2X = np.array([atom.x for atom in pro2.atoms])
    cdef np.ndarray[float, ndim=1] pro2Y = np.array([atom.y for atom in pro2.atoms])
    cdef np.ndarray[float, ndim=1] pro2Z = np.array([atom.z for atom in pro2.atoms])

    # Determine if any pairwise distance is not using the nearest image
    wrapX, wrapY, wrapZ = False, False, False
    for i in range(0, numAtoms):
        for j in range(0, numAtoms):
            if (np.abs(pro1X[i] - pro2X[i]) > x/2):
                wrapX = True
            if (np.abs(pro1Y[i] - pro2Y[j]) > x/2):
                wrapY = True
            if (np.abs(pro1Z[i] - pro2Z[j]) > x/2):
                wrapZ = True

            if (wrapX == wrapY == wrapZ == True):
                i = j = numAtoms

    print wrapX, wrapY, wrapZ

    # Wrap coordinates (inefficient right now)
    if wrapX or wrapY or wrapZ:
        for i in range(0, numAtoms):
            if wrapX:
                if (pro1.atoms[i].x < x/2):
                    pro1.atoms[i].x += x
                if (pro2.atoms[i].x < x/2):
                    pro2.atoms[i].x += x
            if wrapY:
                if (pro1.atoms[i].y < y/2):
                    pro1.atoms[i].y += y
                if (pro2.atoms[i].y < y/2):
                    pro2.atoms[i].y += y
            if wrapZ:
                if (pro1.atoms[i].z < z/2):
                    pro1.atoms[i].z += z
                if (pro2.atoms[i].z < z/2):
                    pro2.atoms[i].z += z

    return


def distance(atom1, atom2):
    vec = np.array([atom1.x - atom2.x, atom1.y - atom2.y, atom1.z - atom2.z])
    return np.linalg.norm(np.abs(vec))


def distance_pbc(atom1, atom2):
    atom1Copy, atom2Copy = atom1, atom2
    if (atom1Copy.box != atom2Copy.box):
        print "Box mismatch, exiting."
        exit(1)

    wrapX, wrapY, wrapZ = False, False, False
    x, y, z = atom1Copy.box

    if (np.abs(atom1Copy.x - atom2Copy.x) > x/2):
        wrapX = True
    if (np.abs(atom1Copy.y - atom2Copy.y) > y/2):
        wrapY = True
    if (np.abs(atom1Copy.z - atom2Copy.z) > z/2):
        wrapZ = True

    if wrapX:
        if (atom1Copy.x < x/2):
            atom1Copy.x += x
        if (atom2Copy.x < x/2):
            atom2Copy.x += x
    if wrapY:
        if (atom1Copy.y < y/2):
            atom1Copy.y += y
        if (atom2Copy.y < y/2):
            atom2Copy.y += y
    if wrapZ:
        if (atom1Copy.z < z/2):
            atom1Copy.z += z
        if (atom2Copy.z < z/2):
            atom2Copy.z += z

    return distance(atom1Copy, atom2Copy)
