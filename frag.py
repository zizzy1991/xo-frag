

from mole import *
import math 
from copy import copy
from copy import deepcopy


#CHratio = singleBond["C"]["H"] / singleBond["C"]["C"]

#define a stardard H2
H2 = Molecule("H2")
H2.addAtom(Atom(0,0,0,"H",-1))
H2.addAtom(Atom(0,0,singleBond['H']['H'],'H',-2))


def addHydrogen(c1, c2):
	'''This function determine an added-Hydrogen atoms near c2 position, 
	replacing c2 atom, after breaking c1-c2 bond
	'''
	# c1, c2 should be 2 carbon atom
	# id -> c1 * 100000 + c2
	rat = 1. * singleBond[c1.atomType]["H"] / single[c1.atomType][c2.atomType]
	delx = c2.x - c1.x
	dely = c2.y - c1.y
	delz = c2.z - c1.z
	f = Atom(c1.x + delx * rat, c1.y + dely * rat, c1.z + delz * rat, "H", c1.index * 100000 + c2.index)
	return f


class Fragment(object):
	'''Fragment() class: support merge(|) and intersect(&) operation,
	use check() 
	to automatedly finish this Fragment, incl generate Hydrogen that missing from 
	the fragment or to satisfy valences
	'''
	def __init__(self, parentMole, name = ""):
		''' two params, 1. parent molecule 2. fragment name
		'''
		# i think Fragment do not need bonds=[]
		self.atoms = []  #this denotes skeletonAtoms, no longer all atoms
		# H could also be skeleton atoms
		self.hydrogen = []
		if name == "":
			self.title = parentMole.title
		else:
			self.title = name
		self.molecule = parentMole
		self.isFinished = False 
		self.flag = 1
		self.charge = 0
		self.multiplicity = 1


	def addAtom(self, atm):
		''' add an atom to this fragment,
				non-H atoms -> self.atoms
				H atoms     -> self.hydrogen
		'''
		if atm not in self.atoms and atm.atomType != 1:
			self.atoms.append(atm)
		elif atm not in self.hydrogen and atm.atomType == 1:
			self.atoms.hydrogen(atm)
		return self


	def __str__(self):
		''' output coord string of all atoms
		'''
		s = "%2d %s\n\n%d %d\n" % (self.flag, self.title, self.charge, self.multiplicity)
		for i in self.atoms:
			s = s + ("%s\n" % str(i))
		for i in self.hydrogen:
			s = s + ("%s\n" % str(i))
		return s


	def xout(self):
		sortedAtoms = copy(self.atoms + self.hydrogen)
		sortedAtoms.sort(lambda a, b : cmp(a.index, b.index))
		xo = [a.index + 1 for a in sortedAtoms if a.index < 100000]		
		return str(xo).replace(' ','')

	def __len__(self):
		return len(self.atoms) + len(self.hydrogen) * 0 
	
	def __iter__(self):
		return iter(self.atoms)

	def __or__(self, another):
		''' merge two fragments
		'''
		s = Fragment(self.molecule, "%so%s" % (self.title, another.title))
		for i in self.atoms:
			s.addAtom(i)
		for i in another.atoms:
			if i not in s.atoms:
				s.addAtom(i)
		return s

	def __and__(self, another):
		''' calc intersection of two fragments
		'''
		s = Fragment(self.molecule, "%sa%s" % (self.title, another.title))
		s.flag = -self.flag * another.flag
		for i in self.atoms:
			if i in another.atoms:
				s.addAtom(i)
		return s

	def generateHydrogen(self):
		''' generate Hydrogen to satisfy the valences of all carbon atoms
		'''
		# check all skeletonAtoms, c1 in skeleton, c2 in all atoms. had c1-c2,
		# add hydrogen,   use a BFS way
		for i in self.atoms:
			for x in self.molecule.atoms:
				if x not in self.atoms and x not in self.hydrogen and self.molecule.hasBond(Bond(i, x)):
					if x.atomType == 'H':
						self.hydrogen.append(x)
					elif x.atomType == 'C':
						self.hydrogen.append(addHydrogen(i, x))
						
	def check(self):
		''' in order, not to cut function groups as c=c n-c, c-o, c=o etc.
		that have double bonds or atoms other than CH.
		'''
		# only single-bond between same type atoms can be cut
		i = 0
		rear = len(self.atoms)
		while i < rear:
			for x in self.molecule.atoms:
				if x not in self.atoms and x not in self.hydrogen:
					if (x.atomType != self.atoms[i].atomType) and self.molecule.bondMatrix[self.atoms[i].index][x.index] >= 0.9999:
						# i - x, i or x is connected to a non-H, non-C atom
						self.addAtom(x)
						rear = rear + 1
					elif self.molecule.bondMatrix[self.atoms[i].index][x.index] > 1.0001:
						# i - x, bond order > 1
						self.addAtom(x)
						rear = rear + 1
			i = i + 1
		#self.generateHydrogen()
		self.isFinished = True

	def _addByBond(self, bnd = 0):
		'''
		'''
		front = 0
		b = [0] * len(self)
		while front < len(self):
			for x in self.molecule:
				if x not in self.atoms and x not in self.hydrogen and b[front] < bnd:
					if self.molecule.bondMatrix[x.index][self.atoms[front].index] > 0:
						self.addAtom(x)
						b.append(b[front] + 1)
			front = front + 1
		return self.check()


	def _addByDist(self, dist = 0.001):
		front = 0
		while front < len(self):
			for x in self.molecule:
				if x not in self.atoms and x not in self.hydrogen:
					if self.atoms[front].distance(x) < dist:
						self.addAtom(x)
			front = front + 1
		return self.check()

class Fragmentor(object):
	''' not a real Fragmentor, actually a father class, support common function
	for a detailed fragmentation method. as output, etc.
	'''
	# support formatted output, generate all fragments
	def __init__(self, mole, R = 3, maxsize = 15):
		self.molecule = mole
		self.title = "Fragmentation job:" + mole.title 
		self.frags = []
		self.maxSize = maxsize
		self.R = R

	def output(self):
		result = []
		for i in self.frags:
			result.append(str(i))
		return result

	def xout(self):
		''' give an xo schemefile-like output
			return:
			frag1_title:[frag1.atoms]
			frag2_title:[frag2.atoms]
		'''
		xot = ''
		for i in self.frags:
			xot = xot + 'f' + i.title + ' : ' + i.xout() + '\n'
		return xot


	def __str__(self):
		return "Fragmentor: should use output(self)"

	def _fixH(self):
		#fix H2, R1-H + R2-H - R1-R2 = H-H
		total = 0
		for x in self.frags:
			total += len(x) * x.flag
		if total != len(self.molecule):
			H2_f = Fragment(H2, "H2 Fragment")
			H2_f.atoms = copy(H2.atoms)
			H2_f.isFinished = True
			H2_f.flag = -(total - len(self.molecule)) / 2
			self.frags.append(H2_f)

	def _sortfrags(self):
		for i in self.frags:
			i.sort()

	def _clean(self):
		''' delete fragments embedded by another
		'''
		# self._sortfrags()
		changed = True
		while changed:
			changed = False
			for i, x in enumerate(self.frags):
				for j, y in enumerate(self.frags):
					xisiny = True
					if j != i:
						for s in x:
							xisiny = xisiny and (s in y)
						if xisiny:
							break
					xisiny = False
				if xisiny:
					self.frags.remove(x)
					changed = True
					break


class MTAfrag(Fragment):
	def __init__(self, parentMole, title = ''):
		super.__init__(parentMole, title)

	def isMergable(self, Rgood = 3, maxsize = 15):
		return self.size < Rgood and len(self) < maxsize
	
class MTA(Fragmentor):
	''' MTA process:
		1. for each atom, generate an fragment
		2. combine small fragments
	'''
	isFinished = False
	def __init__(self, mole, R = 3, maxsize = 15):
		Fragmentor.__init__(self, mole, R, maxsize)
		

	def _generateFragment(self, idx):
		tmpf = Fragment(self.molecule, "%d" % (idx + 1))
		point = self.molecule.atoms[idx]
		tmpf.addAtom(point)
		for x in self.molecule.atoms:
			if x is not point and x.distance(point) < self.R:
				tmpf.addAtom(x)
		self.frags.append(tmpf)

	def _autoCombine(self):
		''' Combine small fragments
		'''
		# 1. find the smallest frag (with fewest atoms)
		# 2. combine it with another frag, which has most common atoms with it
		# 3. repeat 1.2 until size limitation
		l = 0
		while l != len(self.frags):
			l = len(self.frags)
			for tmpf in self.frags:
				mx = 0
				k = 0
				for i, x in enumerate(self.frags):
					c = x & tmpf
					if len(c) > mx and x is not tmpf:
						mx = len(c)
						k = x
				if k != 0 and (len(tmpf | k) <= self.maxSize or mx == len(tmpf) or mx == len(k)):
					self.frags.remove(tmpf)
					self.frags.remove(k)
					self.frags.append(tmpf | k)
					self.frags[-1].check()
					break


	def Run(self):
		self.frags = []
		for i, atm in enumerate(self.molecule.atoms):
			if atm.atomType != 'H':
				self._generateFragment(i)
				self.frags[-1].check()
		self._autoCombine()
		self._clean()

	def __str__(self):
		return "MTA: should use output() or xout()"


class iFragment(Fragment):
	def __init__(self, parentMole, name, forceLst, Rmax):
		Fragment.__init__(self, parentMole, name)
		self.Rmax = Rmax
		self.forceList = forceLst
	
	def check(self):
		''' in order, not to cut function groups as c=c n-c, c-o, c=o etc.
		that have double bonds or atoms other than CH.
		'''
		# TODO  size limitation
		
		front = 0
		while front < len(self.atoms):
			for x in self.molecule.atoms:
				order = self.molecule.bondMatrix[self.atoms[front].index][x.index]
				if order > 0 and x not in self.atoms and x not in self.hydrogen and self.atoms[0].distance(x) < self.Rmax:
					if x.atomType != self.atoms[front].atomType or order > 1.001:
						self.atoms.append(x)
					elif self.forceList[self.atoms[front].index] > self.forceList[x.index]:
						self.atoms.append(x)
			front = front + 1
		#Fragment.check(self)

class iFrag(Fragmentor):
	''' an iFrag Fragmentor:
	procedures:
	1. pick atom with largest force(-error)
    2. generate fragment by placing a sphere of Rmin at the center of atom step1 picked
    3. expand the fragment according to the forcelist,   in a bfs way
       acceptable: safe bond(c-c si-o si-si etc.),  and have smaller force(or smaller enough)
       until: reach unsafe bond, strong force, too many atoms
    4. remove the atom step1 picked
    5. repeat 1234 until whole molecule is covered and all atoms with force > RMS / Fratio are picked
       default: Fradio = 1
    6. merge fragments (in order to reduce amounts) (maybe this step is no need)
    7. calc intersection
    8. check H-bug
	'''
	def __init__(self, mole, f = 3.5, ra = 3.):
		''' should have 4 params, 
		mole: initial molecule
		f   : the times in the complexity expression
		ra  : radius that the fragment need to expand from its core
		'''
		Fragmentor.__init__(self, mole, 0, 0)
		self.Rmin = ra /  (f - 1) * singleBond['C']['C'] + 0.05
		# Rmin, smallest radius that keeps force data stable
		self.Rmax = self.Rmin + ra * singleBond['C']['C'] + 0.05
		self.forceList = []
		# self.forceList list of tuples (id, (x,y,z)^2)
		self.frags = []
		# sum([len(x) for x in self.frags]) / len(self.molecule)
		# update when fragmentation changed. 
		self.coverage = 0.
		# set of atoms that excluded from atom set which we do not pick
		# unstable atom					
		self.excludedAtoms = set([])
		# set of atoms that has been part of fragment
		# this should equal the whole molecule when fragmentation finished
		self.calcedAtoms = set([])

	def loadForceMatrix(self, forceMat):
		''' load force(-error) matrix from param: forceMat,
		    forceMat should be generated from fio.forceMatrixFromGaussOutput
			or fio.forceMatrixFromGaussOutput(filename)
			or fio.forceErrorMatrixFromForceMatrix(fmat1, fmat2)
		'''
		self.excludedAtoms = set([])
		self.forceList = [x[1]**2 + x[2]**2 + x[3]**2 for x in forceMat]
		#for i in self.forceList: print i

	def _getUnstableAtom(self):
		''' return ID of atom that has strongest force(-error)
		'''
		mxForce = -1
		mxI = -1
		for i, x in enumerate(self.forceList):
			if x > mxForce and self.molecule.atoms[i] not in self.excludedAtoms and self.molecule.atoms[i].atomType != 'H':
				mxI, mxForce = i, x
		return mxI

	def _generateFragment(self, idx):
		''' return a fragment generated aroung the No.idx atom
			follow step of 2 and 3
			2. put a sphere at the center of No.idx atom
			3. expand the fragment, expand the size and check whether the fragment
			is appropriate
			4. exclude atoms at distance < r = a / (f - 1)
		'''
		k = idx
		if k == -1: return 0
		# all atoms coverred
		# print k
		self.frags.append(iFragment(self.molecule, '%d' % k, self.forceList, self.Rmax))
		for i in self.molecule.atoms:
			# include atoms in this fragment
			if i.distance(self.molecule.atoms[k]) < self.Rmax:
				self.frags[-1].addAtom(i)
			# exclude atoms from the core set
			if i.distance(self.molecule.atoms[k]) < self.Rmin:
				self.excludedAtoms.add(i)
		self.frags[-1].check()
		return 1

	def Run(self):
		''' repeat:
		    1. pick unstable atom
			2. generete fragment
			until: all atoms are in self.excludedAtoms
		'''
		self.excludedAtoms = set([])
		self.frags = []
		p = self._getUnstableAtom();
		while p >= 0:
			self._generateFragment(p)
			p = self._getUnstableAtom()

	def calcIntersection(self):
		''' calc intersections of all frags,
			doFragmentation(self) should be done before this
			algorithm:
			  1. empty the intersection 
			  2. for all fragments do this.addIntersection
		'''
		# use self.intersection = [...] to store intersection 
		# meanwhile calc coverage
		# NOTE intersection are excluded when we check coverage
		# this can be copied from former-MTA-fragmentor
		pass

	def _addIntersection(self, idOfFragment):
		''' when new block of fragmentation is added, calc intersection
			algorithm: scan over all fragments, incl. intersection already calced
		'''
		# new flag = -fragment.flag * this.flag
	    # NOTE in fact generating a XO-job, pointing out intersection apparently
		# explicitly is no need
		pass

	def __str__(self):
		if self.isFinished: return Fragmentor.__str__(self)
		return "iFrag %s, not finished yet" % self.title

			

class FragForOpt(Fragmentor):
	'''Fragmentation Method ONLY for opt
	'''
	def __init__(self, mole, f = 3.5, a = 8):
		self.molecule = mole
		self.rs = math.ceil(0.5 * a/(f - 1)) * (singleBond['C']['C'] + 0.01)
		# r is minimal radius
		self.rl = self.rs + 0.5 * a * singleBond['C']['C']
		self.frags = []
		self.atombelong = [-1] * len(self.molecule)

	def _generateFragment(self):
		while True:
			mx = 0
			addf = None
			center = 0
			for i in self.molecule.atoms:	
				if self.atombelong[i.index] < 0 and i.atomType != 'H':
					tmpf = Fragment(self.molecule, '%d' % (i.index + 1))
					for j in self.molecule.atoms:
						if j.distance(i) < self.rl:
							# this atom is included in this fragment
							tmpf.addAtom(j)

					tmpf.check()
					if len(tmpf) > mx :
						mx = len(tmpf)
						addf = tmpf
						center = i.index
			if mx != 0:
				self.frags.append(addf)
				self.atombelong[center] = len(self.frags) - 1
				for j in self.molecule.atoms:
					if j.distance(self.molecule.atoms[center]) < self.rs:
						# this atom belongs to the core of the fragment
						self.atombelong[j.index] = len(self.frags) - 1
			else: break
	
	def Run(self):
		self._generateFragment()
		self._clean()

	def __str__(self):
		st = ''
		for i in self.frags:	
			st = st + (i.title + ' :' + i.xout() + '\n')
		return st	
	

class Frouper(object):
	''' fragmentation method based on dividing groups
	'''
	def __init__(self, mole):
		self.groups = Groups(mole)
		self.groups.groupize()
		self.molecule = mole
		self.frags = []

	def __str__(self):
		return "< Frouper based on molecule: %s >" % self.molecule.title

	def xout(self):
		s = ''
		for i in self.frags:
			i.sort(lambda x, y : cmp(x.index0, y.index0))
			l = 'f'
			for j in i:
				l = l + 'o%s' % (j.index0 + 1)
			l = l + ' : ['
			for j in i[:-1]:
				l = l + j.xout() + ','
			l = l + i[-1].xout() + ']'
			s = s + l + '\n'
		return s

	def _addByBond(self, grpIdx, maxdist = 0):
		tmpf = [self.groups.groups[grpIdx]]
		dist = [0]
		front = 0
		while front < len(tmpf):
			for x in self.groups.groups:
				if x not in tmpf and dist[front] < maxdist :
					if self.groups.bondMatrix[x.index][tmpf[front].index] > 0:
						tmpf.append(x)
						dist.append(dist[front] + 1)
			front = front + 1
		self.frags.append(tmpf)

	def Run(self):
		pass

	def _sortfrags(self):
		for i in self.frags:
			i.sort()

	def _clean(self):
		''' delete fragments embedded by another
		'''
		self._sortfrags()
		changed = True
		while changed:
			changed = False
			for i, x in enumerate(self.frags):
				for j, y in enumerate(self.frags):
					xisiny = True
					if j != i:
						for s in x:
							xisiny = xisiny and (s in y)
						if xisiny:
							break
					xisiny = False
				if xisiny:
					self.frags.remove(x)
					changed = True
					break

class SFM(Frouper):
	''' recursive?
	'''
	def __init__(self, mole, lvl = 0):
		Frouper.__init__(self, mole)
		self.level = lvl
		self.frags = []
		self.isRing = [False] * len(self.groups)


	def _path(self, start, groupset):
		''' a BFS Shortest Path algorithm, from a single point : start
			allowed pass-way points are included in groupset
			DETERMINE rings, which are not allowed to be cut-off during the fragmentation
		'''
		dist = [[]] * (len(self.groups))
		front = 0
		dq = [start]
		dist[start] = [start]
		while front < len(dq):
			for j, b in enumerate(self.groups.bondMatrix[dq[front]]):
				if b > 0 and j in groupset and j != dq[front]:
					if dist[j] == []:
						dist[j] = dist[dq[front]] + [j]
						dq.append(j)
					elif dq[front] not in dist[j]:
						# in this case, has a ring.
						# ring member dq[front] - j - others...
						self.isRing[j] = True
						self.isRing[dq[front]] = True
						rb = copy(dist[j])
						rb.reverse()
						# find start point of this ring
						#  - k - ... - dq[front] - j - ... - k -
						for k in rb:
							if k in dist[dq[front]]:
								break
						for x in dist[j][dist[j].index(k):]:
							self.isRing[x] = True
						for x in dist[dq[front]][dist[dq[front]].index(k):]:
							self.isRing[x] = True
					else:
						print "This line never runs"

			front = front + 1
		return dist

	def _correctFrag(self, f, st):
		''' after breaking deleting an atom, figure out the rest part of the fragment
		'''
		q = [st]
		front = 0
		while front < len(q):
			for i, b in enumerate(self.groups.bondMatrix[q[front]]):
				if b > 0 and i in f and i not in q:
					q.append(i)
			front = front + 1
		return q

	def _blocks(self):
		''' in case a molecule has multiple parts. use a linking-blocks to determine
			how many parts which will be fragmented individually next
		'''
		bks = []
		flag = [False] * len(self.groups)
		for i, g in enumerate(self.groups.groups):
			if not flag[i]:
				blk = [i]
				d = self._path(i, range(len(self.groups)))
				flag[i] = True
				for j, g2 in enumerate(d):
					if g2 != [] and not flag[j]:
						blk.append(j)
						flag[j] = True
				blk.sort()
				bks.append(blk)
		return bks

	def Run(self):
		q = self._blocks()
		while q != []:
			f = True
			for i in q[0]:	
				d = self._path(i, q[0])
				# print d
				for j, dj in enumerate(d):
					if len(dj) == self.level + 2:
						frag1 = copy(q[0])
						frag1.remove(dj[-1])
						frag2 = copy(q[0])
						frag2.remove(dj[dj.index(i)])
						frag1 = self._correctFrag(frag1, i)
						frag2 = self._correctFrag(frag2, j)
						frag1.sort()
						frag2.sort()
						q.append(frag1)
						q.append(frag2)
						f = False
						break
				if not f: break		
			if f:
				# add q to self.frags[]
				self.frags.append(q[0])
			q = q[1:]
		self._clean()

	def xout(self):
		s = ''
		#for i in self.frags:
		#	a = ''
		#	for j in i:
		#		a = a + self.groups.groups[j].xout() + ','
		#	print a
		# print '============================================================'
		for i in self.frags:
			# i.sort(lambda x, y : cmp(x.index0, y.index0))
			i.sort()	
			l = 'f'
			l = l + '%s' % (self.groups.groups[i[0]].index0 + 1)
			for j in i[1:]:
				l = l + 'o%s' % (self.groups.groups[j].index0 + 1)
			l = l + ' : ['
			for j in i[:-1]:
				l = l + self.groups.groups[j].xout() + ','
			l = l + self.groups.groups[i[-1]].xout() + ']'
			s = s + l + '\n'
		return s

class CFM(Frouper):
	''' a CFM fragmentor
	'''
	def __init__(self, mole, r = 0):
		Frouper.__init__(self, mole)
		self.rd = r
	
	def Run(self):
		for i in self.groups:
			self._addByBond(i.index, self.rd)
		self._clean()
	
	def __str__(self):
		return "< CFM fragmentor for molecule %s >" % self.molecule.title


