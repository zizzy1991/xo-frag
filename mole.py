

import math
from collections import defaultdict
from copy import copy

# some chem info, bond length, etc. copy from GUO?
ElementsTable = ['', 
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 
    'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 
    'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 
    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'X']
AtomRad=(0.0E0,
        0.643E0, 0.643E0,2.457E0, 1.909E0, 1.587E0, 1.436E0, 1.209E0,
        1.096E0, 1.020E0,0.945E0, 2.986E0, 2.646E0, 2.400E0, 2.192E0,
        2.060E0, 1.890E0,1.795E0, 1.701E0, 3.836E0, 3.288E0, 2.721E0,
        2.494E0, 2.305E0,2.230E0, 2.211E0, 2.211E0, 2.192E0, 2.173E0,
        2.211E0, 2.362E0,2.381E0, 2.305E0, 2.268E0, 2.192E0, 2.154E0,
        2.116E0, 4.082E0,3.609E0, 3.061E0, 2.740E0, 2.532E0, 2.457E0,
        2.400E0, 2.362E0,2.362E0, 2.419E0, 2.532E0, 2.797E0, 2.721E0,
        2.665E0, 2.646E0,2.570E0, 2.513E0, 2.476E0, 4.441E0, 3.742E0,
        3.194E0, 3.118E0,3.118E0, 3.099E0, 3.080E0, 3.061E0, 3.496E0,
        3.042E0, 3.005E0,3.005E0, 2.986E0, 2.967E0, 2.948E0, 2.948E0,
        2.948E0, 2.721E0,2.532E0, 2.457E0, 2.419E0, 2.381E0, 2.400E0,
        2.457E0, 2.532E0,2.816E0, 2.797E0, 2.778E0, 2.759E0, 2.759E0,
        2.740E0
     )
BondRanges = (
#      H-H    He-H    He-He   Li-H    Li-He   Li-Li   Be-H
      0.730E0,0.000E0,0.000E0,1.636E0,0.000E0,2.812E0,1.334E0 ,
#      Be-He   Be-Li   Be-Be    B-H    B-He    B-Li    B-Be
      0.000E0,2.469E0,2.123E0,1.188E0,0.000E0,2.234E0,1.903E0 ,
#       B-B     C-H    C-He    C-Li    C-Be     C-B     C-C
      1.679E0,1.084E0,0.000E0,2.001E0,1.699E0,1.574E0,1.528E0 ,
#       N-H    N-He    N-Li    N-Be     N-B     N-C     N-N
      1.002E0,0.000E0,1.750E0,1.503E0,1.389E0,1.453E0,1.413E0 ,
#       O-H    O-He    O-Li    O-Be     O-B     O-C     O-N
      0.947E0,0.000E0,1.592E0,1.377E0,1.344E0,1.400E0,1.404E0 ,
#       O-O     F-H    F-He    F-Li    F-Be     F-B     F-C
      1.396E0,0.911E0,0.000E0,1.555E0,1.366E0,1.313E0,1.365E0 ,
#       F-N     F-O     F-F    Ne-H    Ne-He   Ne-Li   Ne-Be
      1.386E0,1.376E0,1.345E0,0.000E0,0.000E0,0.000E0,0.000E0 ,
#      Ne-B    Ne-C    Ne-N    Ne-O    Ne-F    Ne-Ne   Na-H
      0.000E0,0.000E0,0.000E0,0.000E0,0.000E0,0.000E0,1.914E0 ,
#      Na-He   Na-Li   Na-Be   Na-B    Na-C    Na-N    Na-O
      0.000E0,2.999E0,2.742E0,2.537E0,2.324E0,2.080E0,1.921E0 ,
#      Na-F    Na-Ne   Na-Na   Mg-H    Mg-He   Mg-Li   Mg-Be
      1.885E0,0.000E0,3.189E0,1.718E0,0.000E0,2.847E0,2.531E0 ,
#      Mg-B    Mg-C    Mg-N    Mg-O    Mg-F    Mg-Ne   Mg-Na
      2.319E0,2.106E0,1.894E0,1.756E0,1.730E0,0.000E0,3.087E0 ,
#      Mg-Mg   Al-H    Al-He   Al-Li   Al-Be   Al-B    Al-C
      2.912E0,1.584E0,0.000E0,2.693E0,2.373E0,2.151E0,1.972E0 ,
#      Al-N    Al-O    Al-F    Al-Ne   Al-Na   Al-Mg   Al-Al
#     1.771E0,1.697E0,1.640E0,0.000E0,2.963E0,2.771E0,2.613E0 ,
      1.771E0,1.797E0,1.640E0,0.000E0,2.963E0,2.771E0,2.613E0 ,
#      Si-H    Si-He   Si-Li   Si-Be   Si-B    Si-C    Si-N
      1.475E0,0.000E0,2.524E0,2.217E0,2.040E0,1.888E0,1.724E0 ,
#      Si-O_1.647 Si-F    Si-Ne   Si-Na   Si-Mg   Si-Al   Si-Si
      1.647E0,1.594E0,0.000E0,2.815E0,2.617E0,2.478E0,2.352E0 ,
#       P-H    P-He    P-Li    P-Be     P-B     P-C     P-N
      1.403E0,0.000E0,2.375E0,2.075E0,1.902E0,1.860E0,1.706E0 ,
#       P-O     P-F    P-Ne    P-Na    P-Mg    P-Al    P-Si
      1.650E0,1.599E0,0.000E0,2.683E0,2.479E0,2.341E0,2.266E0 ,
#       P-P     S-H    S-He    S-Li    S-Be     S-B     S-C
      2.214E0,1.326E0,0.000E0,2.191E0,1.917E0,1.791E0,1.818E0 ,
#       S-N     S-O     S-F    S-Ne    S-Na    S-Mg    S-Al
      1.695E0,1.654E0,1.612E0,0.000E0,2.515E0,2.317E0,2.196E0 ,
#      S-Si     S-P     S-S    Cl-H    Cl-He   Cl-Li   Cl-Be
      2.151E0,2.127E0,2.063E0,1.266E0,0.000E0,2.072E0,1.810E0 ,
#      Cl-B    Cl-C    Cl-N    Cl-O    Cl-F    Cl-Ne   Cl-Na
      1.754E0,1.785E0,1.732E0,1.670E0,1.613E0,0.000E0,2.397E0 ,
#      Cl-Mg   Cl-Al   Cl-Si   Cl-P    Cl-S    Cl-Cl
      2.211E0,2.111E0,2.068E0,2.072E0,2.034E0,1.990E0
)

AtomBondingRange = defaultdict(lambda : range(9))
AtomBondingRange["H"] = (1,)
AtomBondingRange["O"] = (2, 6)
AtomBondingRange["F"] = (1,)
AtomBondingRange["Si"] = (4,)
AtomBondingRange["Al"] = (3, )
AtomBondingRange["C"] = (2, 3, 4)
AtomBondingRange["P"] = (3, 5)

ele1, ele2 = 1, 1 
singleBond = {}
for i in ElementsTable:
	singleBond[i] = {}
for i in BondRanges:
	if i > 0.0000001:
		singleBond[ElementsTable[ele1]][ElementsTable[ele2]] = i
		singleBond[ElementsTable[ele2]][ElementsTable[ele1]] = i
	if ele2 >= ele1:
		ele1 = ele1 + 1
		ele2 = 1
	else: ele2 = ele2 + 1


class Vector(object):
	def __init__(self, x, y, z):
		self.x, self.y, self.z = (x, y, z)

	def __str__(self):
		return "%18.8f%18.8f%18.8f" % (self.x, self.y, self.z)

	@property
	def xyz(self):
		return self.x, self.y, self.z
	@xyz.setter
	def xyz(self, another):
		self.x, self.y, self.z = another

	@property
	def length(self):
		return math.sqrt(self.x**2 + self.y**2 + self.z**2)
	
	def distance(self, another):
		return math.sqrt((self.x - another.x)**2 + (self.y - another.y)**2 + (self.z - another.z)**2)

	def dot(self, another):
		return self.x * another.x + self.y * another.y + self.z * another.z

	def cross(self, another):
		xx = self.y * another.z - self.z * another.y
		yy = self.z * another.x - self.x * another.z
		zz = self.x * another.y - self.y * another.x
		result = Vector(xx, yy, zz)
		return result

	def __add__(self, another):
		xx = self.x + another.x
		yy = self.y + another.y
		zz = self.z + another.z
		return Vector(xx, yy, zz)

	def __sub__(self, another):
		xx = self.x - another.x
		yy = self.y - another.y
		zz = self.z - another.z
		return Vector(xx, yy, zz)

	
class Coord(Vector):
	def __init__(self, x, y, z):
		self.x, self.y, self.z = (x,y,z)
	def __str__(self):
		return "%18.8f%18.8f%18.8f" % (self.x, self.y, self.z)


class Atom(Coord):
	def __init__(self, x = 0, y = 0, z = 0, tp = "H", idx = 0):
		''' new Atom(x,y,z,atomType,idx)
		'''
		Coord.__init__(self, x, y, z)
		if tp not in ElementsTable:
			self.atomType = ElementsTable[tp]
		else: self.atomType = tp
		self.index = idx

	def __str__(self):
		return "%s %18.8f%18.8f%18.8f" % (self.atomType, self.x, self.y, self.z)

class Bond(object):
	def __init__(self, atom1, atom2, bondType = 0):
		if atom1.index < atom2.index:
			self.atom1, self.atom2 = atom1, atom2
		else:
			self.atom1, self.atom2 = atom2, atom1
		if bondType > 0: 
			self._order = bondType
		else: self._order = Bond.calc_order(atom1, atom2)
			
	@property
	def length(self):
		return self.atom1.distance(self.atom2)

	@staticmethod
	def calc_order(atom1, atom2):
		""" calculate bond order according to distance"""
		bo = 0
		if atom2.atomType in singleBond[atom1.atomType]:
			b = singleBond[atom1.atomType][atom2.atomType]
		else: b = 2
		# TODO is that reasonable?!!  and support c=c c=o etc
		#if atom1.distance(atom2) < (b * 1.1)**2:
		if atom1.distance(atom2) < (b * 1.21):
			bo = 1
		return bo
	
	@property
	def order(self):
		return self._order

	def __str__(self):
		return "%3d %3d %3f" % (self.atom1.index + 1, self.atom2.index + 1, self.order)
	

class Molecule(object):
	""" contains atoms and bonds
	"""
	def __init__(self, name = "initial molecule"):
		self.atoms = []
		self.atomsMap = []
		self.bonds = []
		self.charge = 0
		self.multiplicity = 1
		self.title = name
		self.bondMatrix = []
		self.isFinished = False

	def addAtom(self, atom):
		if atom not in self.atoms:
			self.atoms.append(atom)
			self.atomsMap.append(str(atom))
			self.bondMatrix.append([])
			for x in self.bondMatrix: x.append(0)
			self.bondMatrix[-1] = [0] * len(self.atoms)

	def addBond(self, bond):
		if bond.atom1 in self.atoms and bond.atom2 in self.atoms:
			if bond not in self.bonds:
				self.bonds.append(bond)
				self.bondMatrix[bond.atom1.index][bond.atom2.index] = bond.order
				self.bondMatrix[bond.atom2.index][bond.atom1.index] = bond.order
				
	def __len__(self):
		''' return len(self.atoms), amount of atoms
		'''
		return len(self.atoms)

	def __iter__(self):
		''' return iter(self.atoms)
		'''
		return iter(self.atoms)

	def _rebond(self):
		''' re-build all bonds from coords data
			note: use it only when no bond data loaded
		'''
		self.bonds = []
		for i in range(len(self.atoms)):
			for j in range(i + 1, len(self.atoms)):
				order = Bond.calc_order(self.atoms[i], self.atoms[j])
				if order > 0:
					self.addBond(Bond(self.atoms[i], self.atoms[j], order))
	
	@property
	def center(self):
		xyzlist = [(a.x, a.y, a.z) for a in self.atoms]
		sx, sy, sz = map(sum, zip(*xyzlist))
		n = len(self.atoms)
		return Coord(sx*1.0/n, sy*1.0/n, sz*1.0/n)
	
	@property
	def size(self):
		c = self.center
		mx = max([c.distance(a) for a in self.atoms])
		return mx

	def __str__(self):
		# generate a longstring containing cartisian coordinates atom
		# types and charge, etc. in a gausian-job-file-like way
		# this can be modified in sub-class
		
		result = "%s\n%2d%3d\n" % (self.title, self.charge, self.multiplicity)
		for i in self.atoms:
			result = result + str(i) + "\n"
		result += "\n"
		for i in self.bonds:
			result = result + str(i) + "\n"
		result += "\n"
		return result


	def hasBond(self, aBond):
		return self.bondMatrix[aBond.atom1.index][aBond.atom2.index] > 0

	def hasAtom(self, anAtom):
		return str(anAtom) in self.atomsMap	
	
	def finish(self):
		# generate bondmap atommap bonds
		if len(self.bonds) == 0: self._rebond()
		#self.atomsMap = [str(x) for x in self.atoms]
		#self.bondsMap = [str(x) for x in self.bonds]
		self.isFinished = True

class Group(Atom):
	""" build unfragmentable group
		Currently support non X-X only, e.g. X=X, X#X, X-OH, -COOH, -CO-NH
		only single bond between atoms with same type can be cut
	"""
		# TODO RINGS!!!
	def __init__(self, parentMole, origIndex, curIndex):
		''' origIndex : index of the starting atom
			curIndex : index of this group
		'''
		self.molecule = parentMole
		a = self.molecule.atoms[origIndex]
		self.atoms = [a]
		self.index0 = origIndex
		Atom.__init__(self, a.x, a.y, a.z, a.atomType, curIndex)


	def xout(self):
		c = [a.index + 1 for a in self.atoms]
		c.sort()
		s = str(c)
		s = s.replace(' ','')
		return s[1:-1]
		

	def __str__(self):
		s = "<Group:"
		for i in self.atoms:
			s = s + i.atomType
		return s + " no.%d>" % self.atoms[0].index

	def __len__(self):
		return len(self.atoms)

	def __iter__(self):
		return iter(self.atoms)
	
	@property
	def center(self):
		a = self.atoms[0]
		return Coord(a.x*1.0/n, a.y*1.0/n, a.z*1.0/n)
	
	@property
	def size(self):
		c = self.atoms[0]
		mx = max([c.distance(a) for a in self.atoms])
		return mx

	def expand(self):
		q = [self.atoms[0]]
		while q != []:
			# check connected atoms
			for x in self.molecule.atoms:
				if (x not in self.atoms) and (self.molecule.bondMatrix[x.index][q[0].index] > 1 or (x.atomType != q[0].atomType and self.molecule.bondMatrix[x.index][q[0].index] == 1)):
					self.atoms.append(x)
					q.append(x)
			# TODO check rings!!

			q = q[1:]
	
class Groups(object):
	""" split a molecule into several unfragmentable groups and build a
		linkage matrix of these groups
		steps:
			1. pick an atom
			2. generate a group, exclude these atom from the molecule
			3. repeat 12
	"""
	def __init__(self, parentMole, title = '', singleAtomGroup = False):
		self.molecule = parentMole
		self.bondMatrix = []
		self.groups = []
		self.excludedAtoms = []
		self.searchAtom = [-1] * len(self.molecule)
		self.title = title
		self.singleAtomGroup = singleAtomGroup

	def groupize(self):
		if self.singleAtomGroup:
			self.groups = copy(self.molecule.atoms)
			self.groups = deepcopy(self.molecule.bondMatrix)
			return self.groups
		k = 0
		for i in self.molecule:
			if i not in self.excludedAtoms:
				self.groups.append(Group(self.molecule, i.index, k))
				self.groups[-1].expand()
				for j in self.groups[-1]:
					self.excludedAtoms.append(j)	
					self.searchAtom[j.index] = k
				k = k + 1
		self.bondMatrix = []
		for i in range(k):
			self.bondMatrix.append([0] * k)

		for i in self.molecule.bonds:
			x = i.atom1.index
			y = i.atom2.index
			if self.searchAtom[x] != self.searchAtom[y]:
				# print x, y
				self.bondMatrix[self.searchAtom[x]][self.searchAtom[y]] = 1
				self.bondMatrix[self.searchAtom[y]][self.searchAtom[x]] = 1
		return self.groups

	def xout(self):
		s = ''
		for i in self.groups:
			s = s + i.xout() + '\n'
		
		s = s + str(self.searchAtom) + '\n'
		#return s
		for i, x in enumerate(self.bondMatrix):
			for j, y in enumerate(x):
				if self.bondMatrix[i][j] > 0:
					s = s + "%d %d\n" % (i, j)
		return s
	
	def __iter__(self):
	    return iter(self.groups)

	def __str__(self):
		return "< Groups: with %d groups>" % len(self.groups)
		
	def __len__(self):
	    return len(self.groups)


