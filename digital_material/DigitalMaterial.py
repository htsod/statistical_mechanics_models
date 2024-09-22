################################
# Digital Material 
#
# ListOfAtoms
# Initializers
# Transformers
# Observers
# NeighborLocators
# BoundaryConditions
# Potentials
# Movers
# MDSystem
#
################################

#########################################
# Import the library(s)
#########################################
# Windows
# import matplotlib.pyplot as plt
# import matplotlib
# matplotlib.use("QTAgg")
# plt.plot([1,2],[3,4])
# plt.show()


# Linux
import vpython as vi
# kludge required to get visual window up and running before import np
c = vi.cone(radius=1.0e-10)
c.visible = 1


# vi.scene.caption = """Right button drag or Ctrl-drag to rotate "camera" to view scene.
# To zoom, drag with middle button or Alt/Option depressed, or use scroll wheel.
#     On a two-button mouse, middle is left + right.
# Shift-drag to pan left/right and up/down.
# Touch screen: pinch/extend to zoom, swipe or two-finger rotate."""




import numpy as np
RandomArray=np.random
#import psyco
#psyco.full()

#########################################
# ListOfAtoms
#########################################

class ListOfAtoms:
    #
	def __init__(self, mass=1.0, radius=0.5, color=vi.color.green, positions=None, velocities=None):
		"""
		Sets up blank list of atoms
		"""
		self.positions = positions
		self.velocities = velocities
		self.mass = mass
		self.radius = radius
		self.color = color


    #
	def KineticEnergy(self):
		"""
		Sums m v^2 / 2:
		"""
		#return 0.5*sum(np.sum(self.mass*self.velocities*self.velocities))
		return 0.5*np.sum(self.mass*self.velocities*self.velocities)



#########################################
# Initializers
#########################################


def RandomListOfAtoms(L, nAtoms=1000, temperature=1.0, 
			mass=1.0, radius=0.5, color=vi.color.green):
    """
    Put atoms of radius r in box (0,L)^dim
    """
    positions = RandomArray.uniform(radius, L-radius, (nAtoms,dim))
    velocities = np.zeros((nAtoms,dim))
    atoms = ListOfAtoms(mass, radius, color, positions, velocities)
    ThermalizingTransformer(atoms, temperature)
    return atoms


def RandomNonoverlappingListOfAtoms(L, neighborLocator, minDist=1.0, 
				nAtoms=10, temperature=1.0, maxTriesQuiet=1000,
				mass=1.0, radius=0.5, color=vi.color.green):
	"""
	Put atoms of given radius in box (0,L)^dim, at random except without
	overlaps less than minDist.
	"""
	pos = RandomArray.uniform(radius, L-radius, (nAtoms,dim))
	vel = np.zeros((nAtoms,dim))
	atoms = ListOfAtoms(mass, radius, color, pos, vel)
	# Enforce no overlaps
	# neighborLocator.HalfNeighbors returns (n1, n2, r, dr) 
	# with n1 > n2: remove n1 atoms
	overlappingAtoms = neighborLocator.HalfNeighbors(atoms, minDist)[0]
	nPairOverlap = len(overlappingAtoms)
	tries = 0
	while nPairOverlap > 0:
		tries += 1
		if tries%maxTriesQuiet == 0:
			print("After ", tries, "attempts, still have ", nPairOverlap, \
			"overlapping pairs of atoms: may need fewer atoms or larger box")
		posNew = [pos for n, pos in enumerate(atoms.positions) \
				if n not in overlappingAtoms]
		nOverlap = nAtoms - len(posNew)
		newAtomPositions = RandomArray.uniform(radius, L-radius, (nOverlap,dim))
		posNew.extend(newAtomPositions)
		atoms.positions=np.array(posNew)
		overlappingAtoms = neighborLocator.HalfNeighbors(atoms, minDist)[0]
		nPairOverlap = len(overlappingAtoms)
	ThermalizingTransformer(atoms, temperature)
	return atoms


def SquareListOfAtoms(L, nAtoms=25, temperature=1.0, 
					  mass=1.0, radius=0.5, color=vi.color.green):
	"""
	Atoms in a square lattice
	"""
	positions = np.zeros((nAtoms,dim), float)
	velocities = np.zeros((nAtoms,dim), float)
	GridSize = int(np.ceil(float(nAtoms)**(1.0/dim)))
	atno = 0
	corner = L/2.0 - GridSize*radius + radius
	# Avoid atom at boundary of box
	if dim==2:
		zeroPosition = np.array([corner, corner])
	if dim==3:
		zeroPosition = np.array([corner, corner, corner])

	for i in range(GridSize):
		for j in range(GridSize):
			if (atno<nAtoms):
				positions[atno] = np.array([i,j])*2.*radius+zeroPosition
				atno+=1
	atoms = ListOfAtoms(mass, radius, color, positions, velocities)
	ThermalizingTransformer(atoms, temperature)
	return atoms


def ClusterListOfAtoms(inside, R, latticeVectors, origin=None, temperature=1.0, 
		center=None, maxIndex=None, 
		mass=1.0, radius=0.5, color=vi.color.green):
	"""
    Fills a volume defined by 'inside(dr, R)' of characteristic volume R with
    atoms with the given lattice vectors. Shifts lattice so one atom is at 
    'origin'
    #
    Creates a lattice of atoms determined by latticeVectors, 
	origin + m[1] lV[1] + m[2] lV[2] {+ m [3] lV[3]}
    filling a sphere of radius R centered at 'center'.
    Assumes lattice isn't very 'skewed', so max index for m is 
    less than 3*R/r, where r is minimum length of three
    lattice vectors.
    """
	if R<0:
		positions = np.array([])
		velocities = np.array([])
		return
	if maxIndex==None:
		rMin = np.sqrt(min( \
	    	[np.dot(latticeVectors[i], latticeVectors[i]) for i in range(dim)]))
	maxIndex = int(3*R/rMin)
	if type(center) == None:
		center = np.zeros(dim, float)
	if type(origin) == None:
		origin = np.zeros(dim, float)
	
	m = np.zeros(dim)
	atomList = []
	mCenter = [int(np.dot(center-origin, latticeVectors[i])) for i in range(dim)]

	for m[0] in range(mCenter[0]-maxIndex, mCenter[0]+maxIndex):
		for m[1] in range(mCenter[1]-maxIndex, mCenter[1]+maxIndex):
			if (dim==2):
				pos = m[0]*latticeVectors[0]+m[1]*latticeVectors[1]
				if inside(pos-center, R):
					atomList.append(pos)
			if (dim==3):
				for m[2] in range(mCenter[2]-maxIndex, mCenter[2]+maxIndex):
					pos = m[0]*latticeVectors[0]+m[1]*latticeVectors[1]+m[2]*latticeVectors[2]
					if inside(pos-center, R):
						atomList.append(pos)

	positions = np.array(atomList)
	velocities = np.zeros(np.shape(atomList), float)
	atoms = ListOfAtoms(mass, radius, color, positions, velocities)
	ThermalizingTransformer(atoms, temperature)
	return atoms


def SphericalClusterListOfAtoms(R, latticeVectors, origin=None, 
		temperature=1.0, center=None, maxIndex=None, 
		mass=1.0, radius=0.5, color=vi.color.green):
    """
    Fills a circle (d=2) or sphere (d=3) centered at 'center' with atoms 
    with the given lattice vectors. Shifts lattice so one atom is at 'origin'
    """
    #
    def inside(dr, R):
        return np.dot(dr, dr) < R*R
    #
    return ClusterListOfAtoms(inside, R, latticeVectors, origin,
    			temperature, center, maxIndex, mass, radius, color)


def CubicalClusterListOfAtoms(L, latticeVectors, origin=None, 
		temperature=1.0, maxIndex=None, 
		mass=1.0, radius=0.5, color=vi.color.green):
	
	"""
    Fills a cube [0,L]x[0,L]x[0,L] with atoms 
    with the given lattice vectors. Shifts lattice so one atom is at 'origin'
    """

	def inside(dr, R):
		#min = -L/2 + radius
		#max = L/2-radius
		min = -L/2 
		max = L/2
		Inside = True
		
		for d in range(dim):
			if dr[d] < min:
				Inside = False
			if dr[d] > max:
				Inside = False
		return Inside
    #
	center = L/2.*np.ones(dim, float)
	return ClusterListOfAtoms(inside, L, latticeVectors, origin,
    			temperature, center, maxIndex, mass, radius, color)


def SphericalClusterListOfAtoms_1(R, latticeVectors, origin=None, 
		temperature=1.0, center=None, maxIndex=None, 
		mass=1.0, radius=0.5, color=vi.color.green):
	"""
    Fills a circle (d=2) or sphere (d=3) centered at 'center' with atoms 
    with the given lattice vectors. Shifts lattice so one atom is at 'origin'
    #
    Creates a lattice of atoms determined by latticeVectors, 
	origin + m[1] lV[1] + m[2] lV[2] {+ m [3] lV[3]}
    filling a sphere of radius R centered at 'center'.
    Assumes lattice isn't very 'skewed', so max index for m is 
    less than 3*R/r, where r is minimum length of three
    lattice vectors.
    """
	if R<0:
		positions = np.array([])
		velocities = np.array([])
		return
	if maxIndex==None:
		rMin = np.sqrt(min([np.dot(latticeVectors[i], latticeVectors[i]) for i in range(dim)]))
		maxIndex = int(3*R/rMin)
	if type(center) == None:
		center = np.zeros(dim, float)
	if type(origin) == None:
		origin = np.zeros(dim, float)

	m = np.zeros(dim)
	atomList = []
	mCenter = [int(np.dot(center-origin, latticeVectors[i]))
		for i in range(dim)]
	
	for m[0] in range(mCenter[0]-maxIndex, mCenter[0]+maxIndex):
		for m[1] in range(mCenter[1]-maxIndex, mCenter[1]+maxIndex):
			if (dim==2):
				pos = m[0]*latticeVectors[0]+m[1]*latticeVectors[1]
				if np.dot(pos-center, pos-center) < R*R:
					atomList.append(pos)
			if (dim==3):
				for m[2] in range(mCenter[2]-maxIndex, mCenter[2]+maxIndex):
					pos = m[0]*latticeVectors[0]+m[1]*latticeVectors[1]+m[2]*latticeVectors[2]
					if np.dot(pos-center, pos-center) < R*R:
						atomList.append(pos)
	positions = np.array(atomList)
	velocities = np.zeros(np.shape(atomList), float)
	atoms = ListOfAtoms(mass, radius, color, positions, velocities)
	ThermalizingTransformer(atoms, temperature)
	return atoms


def TriangularSphericalClusterListOfAtoms(R, latticeSpacing=None,
		 origin=None, 
		 temperature=1.0, center = None, maxIndex=None, 
		 mass=1.0, radius=0.5, color=vi.color.green):
	"""
    Fills a circular region with triangles of atoms
    """
	if latticeSpacing is None:
		latticeSpacing=2.*radius
		latticeVectors = latticeSpacing * np.array([[1.,0.] , [0.5, 0.5*np.sqrt(3.)]])
	atoms = SphericalClusterListOfAtoms(R, latticeVectors, \
	origin, temperature, center, maxIndex, mass, radius, color)
	return atoms

def FCCSphericalClusterListOfAtoms(R, latticeSpacing=None,
		 origin=None, 
		 temperature=1.0, center = None, maxIndex=None, 
		 mass=1.0, radius=0.5, color=vi.color.green):
	"""
    Fills a circular region with fcc lattice of atoms
    """
	if latticeSpacing is None:
		latticeSpacing=2.*radius
	a = latticeSpacing/np.sqrt(2.)
	latticeVectors = np.array([[a,a,0.] , [a,0.,a], [0.,a,a]])
	atoms = SphericalClusterListOfAtoms(R, latticeVectors, origin, temperature, center, maxIndex, mass, radius, color)
	return atoms
	

def TriangularSquareClusterListOfAtoms(L, latticeSpacing=None,
		 origin=None, 
		 temperature=1.0, maxIndex=None, 
		 mass=1.0, radius=0.5, color=vi.color.green):
	"""
    Fills a square [0,L]^3 with triangular lattice of atoms
    """
	if latticeSpacing is None:
		latticeSpacing=2.*radius
	a = latticeSpacing/np.sqrt(2.)
	latticeVectors = latticeSpacing * np.array([[1.,0.] , [0.5, 0.5*np.sqrt(3.)]])
	atoms = CubicalClusterListOfAtoms(L, latticeVectors, origin, temperature, maxIndex, mass, radius, color)
	return atoms


def FCCCubicalClusterListOfAtoms(L, latticeSpacing=None,
		 origin=None, 
		 temperature=1.0, maxIndex=None, 
		 mass=1.0, radius=0.5, color=vi.color.green):
	"""
    Fills a cube [0,L]^3 with fcc lattice of atoms
    """
	if latticeSpacing is None:
		latticeSpacing=2.*radius
	a = latticeSpacing/np.sqrt(2.)
	latticeVectors = np.array([[a,a,0.] , [a,0.,a], [0.,a,a]])
	atoms = CubicalClusterListOfAtoms(L, latticeVectors, origin, temperature, maxIndex, mass, radius, color)
	return atoms
	

#########################################
# Transformers
#########################################

def ThermalizingTransformer(atoms, T):
    """
    Thermalizes velocities to m v^2 / 2 = kB T/2, with kB = 1
    """
    # 
    vRMS = np.sqrt(T/atoms.mass)
    atoms.velocities = RandomArray.normal(0, vRMS, np.shape(atoms.velocities))
  

		
#########################################
# Observers
#########################################

class AtomsObserver:
	"""
    Base class for observers of the atomic simulation
    Observers should not change the state of the simulation, but only
    record, display, or analyze
    Methods which change the atoms class should call Update on a list
    of observables, which can be modified to add or subtract analysis packages
	"""
	def Update(self, atoms, time=None):
		assert False


class VelocityTrajectoryObserver:
	"""
    Keeps track of trajectories for velocities
	"""
	def __init__(self):
		self.Reset()


	def Update(self, atoms, time=None):
		self.vTrajectory.append(atoms.velocities.copy())


	def Reset(self):
		self.vTrajectory = []


	def vXvY(self, atomIndex=0):
		trajNum = np.array(self.vTrajectory)
		return trajNum[:,atomIndex,0], trajNum[:,atomIndex,1]


class UnfoldedTrajectoryObserver:
	"""
    Keeps track of atomic trajectories for periodic boundary conditions,
    unfolding them to avoid jumps at boundaries
    Assumes no atoms move more than 0.5 L
	"""
	def __init__(self, L):
		self.L = L
		self.Reset()
	
	
	def Update(self, atoms, time=None):
		if self.trajectory == []:
			self.trajectory.append(atoms.positions.copy())
			self.oldPos = atoms.positions.copy()
		else:
			dr = atoms.positions - self.oldPos
	    	# Assumes no atoms move more than 0.5 L
			dr = (((dr+1.5*self.L)%self.L)-0.5*self.L)
			self.trajectory.append(self.trajectory[-1]+ dr)
			self.oldPos = atoms.positions.copy()
    
	
	def Reset(self):
		self.trajectory = []
		self.oldPos = None


	def XY(self, atomIndex=0):
		trajNum = np.array(self.trajectory)
		return trajNum[:,atomIndex,0], trajNum[:,atomIndex,1]
	

	def r2Bar(self):
        # XXX Needs to be divided by nAtoms
	# I think nAtoms = self.trajectory.shape()[1]?
		trajDiff = np.array(self.trajectory) - self.trajectory[0]
		return np.sum(np.sum(trajDiff*trajDiff,1),1)
	

class EnergyObserver(AtomsObserver):
	def __init__(self, potential, neighborLocator, boundaryConditions):
		"""
		Stores potential, kinetic, and total energy versus time
		Usage: 	to calculate means and error bar of potential energy ignoring
			the first ten points, use 
			e.Mean(e.PEs), e.SigmaMean(e.PEs)
			to calculate temperature and specific heat from fluctuations
			in the kinetic energy, use
			e.Temperature(), e.SpecificHeat()
		Stores number of atoms and dimension of space; complains if they change;
		Uses nAtoms and dimension to compute temperature
		"""
		self.potential = potential
		self.neighborLocator = neighborLocator
		self.boundaryConditions = boundaryConditions
		self.Reset()


	def Reset(self):
		"""
		Empties stored information
		"""
		self.PEs=[]
		self.KEs=[]
		self.Es=[]
		self.nAtoms = None
		self.dimension = None


	def Update(self, atoms, time=None):
		"""
		Call each time you want to take a snapshot of the energies
		"""
		KE = atoms.KineticEnergy()
		PE = self.potential.PotentialEnergy(atoms, self.neighborLocator)
		E = KE + PE
		self.KEs.append(KE)
		self.PEs.append(PE)
		self.Es.append(E)
		nAtoms, dimension = np.shape(atoms.positions)
		if self.nAtoms is not None:
			assert(self.nAtoms == nAtoms)
			assert(self.dimension == dimension)
		self.nAtoms = nAtoms
		self.dimension = dimension


	def Mean(self, vec):
		"""
		Returns mean of Es, PEs, KEs, etc.
		"""
		assert len(vec)>0
		return sum(vec)/float(len(vec))
	

	def Sigma(self, vec):
		"""
		Returns RMS fluctuations of Es, PEs, KEs, etc.
		"""
		assert len(vec)>1
		mean = self.Mean(vec)
		diff = np.array(vec)-mean
		# Divide by N-1: corrects for imperfections in mean (N=1 has vec-mean=0)
		return np.sqrt(sum(diff*diff)/(len(vec)-1.0))
	

	def SigmaMean(self, vec):
		"""
		Returns confidence in mean for Es, PEs, KEs, etc.
		(Equals RMS fluctuations divided by sqrt(N))
		"""
		assert len(vec)>1
		return self.Sigma(vec)/np.sqrt(len(vec))
	

	def DOF(self, nUncoupledDOF=0):
		"""
		Note: if some degrees of freedom (DOF) are uncoupled 
		(e.g., center-of-mass translation modes for periodic boundary 
		conditions), their kinetic energies will not be correctly incorporated. 
		If their velocities have been set to zero, send the number of these as
		nUncoupledDOF.
		"""
		return self.nAtoms * self.dimension -nUncoupledDOF
	

	def Temperature(self, nUncoupledDOF=0):
		"""
		Calculates temperature kB T from <KE>= DOF * kB T/2
		Note: if some degrees of freedom (DOF) are uncoupled 
		(e.g., center-of-mass translation modes for periodic boundary 
		conditions), their kinetic energies will not be correctly incorporated. 
		If their velocities have been set to zero, send the number of these as
		nUncoupledDOF.
		"""
		return 2.0 * self.Mean(self.KEs) / self.DOF(nUncoupledDOF)
	

	def KineticSpecificHeat(self, nUncoupledDOF=0):
		"""
		Calculates kinetic energy contribution cKE/kB to specific heat per
		particle using
		cKE = (1/2 kB)*dimension
		Note that 
		"""
		return 0.5 * self.DOF(nUncoupledDOF)/self.nAtoms
	

	def PotentialSpecificHeat(self, nUncoupledDOF=0):
		"""
		Calculates potential energy contribution c_PE/kB = C_PE/(N kB) 
		to specific heat per particle from 
		1/C_PE + 1/C_KE = kB T^2 / sigmaKE^2
		(see exercise 3.7, Microcanonical Energy Fluctuations, in Sethna's
		"Statistical Mechanics: Entropy, Order Parameters, and Complexity")
		"""
		cKE = self.KineticSpecificHeat(nUncoupledDOF)
		T = self.Temperature(nUncoupledDOF)
		sigmaKE = self.Sigma(self.KEs)
		cPE = 1.0 / (self.nAtoms * T**2 / sigmaKE**2 - 1.0/cKE)
		return cPE
	

	def SpecificHeat(self, nUncoupledDOF=0):
		"""
		Calculates specific heat per particle from fluctuations in the
		kinetic energy
		"""
		cPE = self.PotentialSpecificHeat(nUncoupledDOF) 
		cKE = self.KineticSpecificHeat(nUncoupledDOF)
		return cPE+cKE


class VisualDisplayAtomsObserver(AtomsObserver):
	def __init__(self, atoms, L, thk=0.5, rate=200, velocityColors=False):
		"""
		Initialize vpython display, for atoms confined to (0,L)x(0,L)
		Draws box of size LxL thickness thk
		Runs at rate = 200
		"""
		self.ball_list = []
		for obj in vi.scene.objects:
			obj.visible=0
		vi.scene.title='Atoms'
		vi.scene.autoscale=True
		# if dim==3:
		# 	vi.scene.range = (L/2.+thk, L/2.+thk, L/2.+thk)
		# elif dim==2:
		# 	vi.scene.range = (L/2.+thk, L/2.+thk, L/2.+thk)
		vi.scene.center=vi.vector(L/2., L/2., L/2)
		vi.rate(rate)
		# Create Wall(s)
		sVert = L 
		sHoriz = L+2.*thk
		if dim==2:
			wallR = vi.box(pos=vi.vector(L+0.5*thk, L/2., 0), 
				  size=vi.vector(thk, sVert, thk), color=vi.color.red) 
			wallL = vi.box(pos=vi.vector(-0.5*thk, L/2., 0), 
				  size=vi.vector(thk, sVert, thk), color=vi.color.red) 
			wallB = vi.box(pos=vi.vector(L/2., -0.5*thk, 0), 
				  size=vi.vector(sHoriz, thk, thk), color=vi.color.blue) 
			wallT = vi.box(pos=vi.vector(L/2., L+0.5*thk, 0), 
				  size=vi.vector(sHoriz, thk, thk), color=vi.color.blue) 
		elif dim==3:
			wallR = vi.box(pos=vi.vector(L+0.5*thk, L/2., L/2.), 
				  size=vi.vector(thk, sVert, L), color=vi.color.red) 
			wallL = vi.box(pos=vi.vector(-0.5*thk, L/2., L/2.), 
				  size=vi.vector(thk, sVert, L), color=vi.color.red) 
			wallB = vi.box(pos=vi.vector(L/2., -0.5*thk, L/2.), 
				  size=vi.vector(sHoriz, thk, L), color=vi.color.blue) 
			wallT = vi.box(pos=vi.vector(L/2., L+0.5*thk, L/2.), 
				  size=vi.vector(sHoriz, thk, L), color=vi.color.blue) 
			wallBK= vi.box(pos=vi.vector(L/2., L/2., -0.5*thk), 
				  size=vi.vector(sHoriz, sHoriz, thk), color=vi.vector(0.7,0.7,0.7)) 
		# vi.scene.autoscale = 0
		# Create Ball(s)
		self.ball_list = []
		for n, pos in enumerate(atoms.positions):
			ball=vi.sphere(color=atoms.color,radius=atoms.radius) 
			self.ball_list.append(ball)
			self.velocityColors = velocityColors
		self.Update(atoms)
    #
	def Update(self, atoms, times=None):
		if self.velocityColors:
			vel = atoms.velocities
			vxrms = np.sqrt(np.sum(vel[:,0]*vel[:,0])/len(vel))
			svel = np.clip(np.fabs(vel)/vxrms, 0., 1.)
			for n, pos in enumerate(atoms.positions):
				if dim==2:
					self.ball_list[n].pos=vi.vector(pos[0], pos[1], 0.0)
					self.ball_list[n].color=vi.vector(svel[n][0], svel[n][1])

				if dim==3:
					self.ball_list[n].pos=vi.vector(pos[0], pos[1], pos[2])
					self.ball_list[n].color=vi.vector(svel[n][0], svel[n][1], svel[n][2])
					
		else:
			for n, pos in enumerate(atoms.positions):
				if dim==2:
					self.ball_list[n].pos=vi.vector(pos[0], pos[1], 0.0)
				if dim==3:
					self.ball_list[n].pos=vi.vector(pos[0], pos[1], pos[2])

		
#########################################
# Neighbor locators
#########################################


class NeighborLocator:
	"""
    Base class for NeighborLocator: uses boundary conditions to find
    neighbors of atoms with in distance dist
    """
	def __init__(self, dist, boundaryConditions):
		self.dist = dist
		self.boundaryConditions = boundaryConditions


	def HalfNeighbors(self, atoms, dist=None):
		assert(False)


class NoNeighborLocator(NeighborLocator):
	"""
    Useful neighbor locator which returns no neighbors
    """
	def HalfNeighbors(self, atoms, dist=None):
		"""
		Returns empty lists for n1, n2, dr and r^2
		"""
		return [], [], np.array([]), np.array([])


class SimpleNeighborLocator(NeighborLocator):
	"""
    Assembles lists by comparing all pairs of atoms. Not recommended for
    large systems: O(N^2). Cell neighbor locators and neighbor lists are
    faster for large systems, but this may compete for small ones.
    """
	
	def HalfNeighbors(self, atoms, dist=None):
		"""
		Returns n1, n2, dr = [dx, dy, dz] = x2-x1, and r^2
		for all pairs of atoms n1>n2 separated by a distance less than dist
		"""
		n1 = []
		n2 = []
		dr = []
		r2 = []
		if (dist!=None):
			distSquared=dist*dist
		else:
			distSquared = self.dist * self.dist
		for n, pos in enumerate(atoms.positions):
			dr1 = atoms.positions[0:n]-pos
			dr1 = self.boundaryConditions.DifferenceBoundaryConditions(dr1)
			dr1sq = np.sum(dr1*dr1,1)
			n1.extend([n
			for m, drsq in enumerate(dr1sq) if drsq < distSquared])
			n2.extend([m
			for m, drsq in enumerate(dr1sq) if drsq < distSquared])
			dr.extend([dr1[m]
			for m, drsq in enumerate(dr1sq) if drsq < distSquared])
			r2.extend([drsq
			for m, drsq in enumerate(dr1sq) if drsq < distSquared])
		# Conversion to np.array takes half the time (for dist->infty)!
		dr = np.array(dr)
		r2 = np.array(r2)
		return n1, n2, dr, r2


class SimpleNeighborLocator_2 (NeighborLocator):
	"""
    Assembles lists by comparing all pairs of atoms. Not recommended for
    large systems: O(N^2). Cell neighbor locators and neighbor lists are
    faster for large systems, but this may compete for small ones.
    """
	
	def HalfNeighbors(self, atoms, dist=None):
		"""
		Returns n1, n2, dr = [dx, dy, dz] = x2-x1, and r^2
		for all pairs of atoms separated by a distance less than dist
		"""
		n1 = []
		n2 = []
		dr = []
		r2 = []
		if (dist!=None):
			distSquared=dist*dist
		else:
			distSquared = self.dist * self.dist
		for n, pos in enumerate(atoms.positions):
			dr1 = atoms.positions[0:n]-pos
			dr1 = self.boundaryConditions.DifferenceBoundaryConditions(dr1)
			dr1sq = np.sum(dr1*dr1,1)
			smaller = dr1sq < distSquared
			n1.extend(n*np.ones(np.sum(smaller)))
			n2.extend(np.compress(smaller, np.arange(len(dr1sq))))
			r2.extend(np.compress(smaller, dr1sq))
			dr.extend(np.compress(smaller, dr1,0))
		# Conversion to np.array takes half the time (for dist->infty)!
		dr = np.array(dr)
		r2 = np.array(r2)
		return n1, n2, dr, r2


class SimpleNeighborLocator_3 (NeighborLocator):
	"""
    Assembles lists by comparing all pairs of atoms. Not recommended for
    large systems: O(N^2). Cell neighbor locators and neighbor lists are
    faster for large systems, but this may compete for small ones.
    """
	
	def HalfNeighbors(self, atoms, dist=None):
		"""
		Returns n1, n2, dr = [dx, dy, dz] = x2-x1, and r^2
		for all pairs of atoms separated by a distance less than dist
		"""
		n1 = []
		n2 = []
		dr = np.zeros((0,2), 'l')
		r2 = np.array([])
		if (dist!=None):
			distSquared=dist*dist
		else:
			distSquared = self.dist * self.dist
		for n, pos in enumerate(atoms.positions):
			dr1 = atoms.positions[0:n]-pos
			dr1 = self.boundaryConditions.DifferenceBoundaryConditions(dr1)
			dr1sq = np.sum(dr1*dr1,1)
			n1.extend([n
			for m, drsq in enumerate(dr1sq) if drsq < distSquared])
			n2.extend([m
			for m, drsq in enumerate(dr1sq) if drsq < distSquared])
			x = np.array([dr1[m] \
					for m, drsq in \
					enumerate(dr1sq) \
					if drsq < distSquared])
			if len(x) > 0:
				dr = np.concatenate((dr, x))
			r2 = np.concatenate((r2, \
						[drsq for m, drsq in enumerate(dr1sq) \
						if drsq < distSquared]))
		# Conversion to np.array takes half the time (for dist->infty)!
		dr = np.array(dr)
		r2 = np.array(r2)
		return n1, n2, dr, r2


#########################################
# Boundary conditions
#########################################


class BoundaryConditions:
	"""
    Boundary conditions base class
	"""
    #
	def __init__(self, L):
		self.L = L
    #
	def EnforceBoundaryConditions(self, atoms):
		assert(False)
    #
	def EnforceBoundaryConditions(self, atoms):
		assert(False)


class FreeBoundaryConditions (BoundaryConditions):
	"""
    Free boundary conditions
    """
    # Store box size for graphics, etc.
	def __init__(self, L):
		BoundaryConditions.__init__(self, L)
    #


	def EnforceBoundaryConditions(self, atoms):
		"""No change"""
		return
	

	def DifferenceBoundaryConditions(self, dr):
		"""No change"""
		return dr


class ReflectiveBoundaryConditions (BoundaryConditions):
	"""
    Reflective boundary conditions
    """
    #
	def __init__(self, L, impulseRecording = False):
		BoundaryConditions.__init__(self, L)
		self.impulseRecording = impulseRecording
		if impulseRecording:
			self.Reset()
    #


	def EnforceBoundaryConditions(self, atoms):
		"""
		If atom outside box, reflect atom inside the box, and reflect velocity
		"""
		maxpos = self.L-atoms.radius
		minpos = atoms.radius
		for n, pos in enumerate(atoms.positions):
			for d in range(dim):
				if pos[d] > maxpos:
					# save impulse if necessary
					if self.impulseRecording:
						self.impulses[d][1].append(
							2.0*atoms.mass*atoms.velocities[n][d])
					# reflect velocity
					atoms.velocities[n][d] *= -1
					# reflect positions
					atoms.positions[n][d] = maxpos-(pos[d]-maxpos)
				if pos[d] < minpos:
					if self.impulseRecording:
						self.impulses[d][0].append(
							2.0*atoms.mass*atoms.velocities[n][d])
					atoms.velocities[n][d] *= -1
					atoms.positions[n][d] = minpos+(minpos-pos[d])

		
	def DifferenceBoundaryConditions(self, dr):
		"""
		Reflective boundaries don't change distances
		"""
		return dr
	

	def Reset(self):
		"""
		Sets up / erases impulse lists
		impulses[d][0] is wall perpendicular to dimension d at zero
		impulses[d][1] is wall perpendicular to dimension d at L
		"""
		if self.impulseRecording:
			self.impulses=[]
			for d in range(dim):
				self.impulses.append([[],[]])


class PeriodicBoundaryConditions (BoundaryConditions):
	"""
    Periodic boundary conditions
    """
    #
	def __init__(self, L):
		BoundaryConditions.__init__(self, L)
    #
		

	def EnforceBoundaryConditions(self, atoms):
		"""
		Put all atoms in L^d box using periodic boundary conditions
		Assumes atoms within 100*L of box
		"""
		# Avoid negative fmod 
		if dim==2:
			shift = np.array([self.L, self.L]) 
		if dim==3:
			shift = np.array([self.L, self.L, self.L]) 
		atoms.positions = np.fmod(atoms.positions+100*shift, self.L)
    #
		
	
	def DifferenceBoundaryConditions(self, dr):
		"""
		Finds nearest copy of neighboring atom
		Shift dr by periodicity until all components of distance
		vectors have absolute value < L/2
		Assumes atoms separated by less than 100*L
		"""
		if dim==2:
			shift = 0.5*np.array([self.L, self.L]) 
		if dim==3:
			shift = 0.5*np.array([self.L, self.L, self.L]) 
		drNew = dr + 201.*shift
		drNew = np.fmod(drNew, self.L)
		drNew -= shift
		return drNew


#########################################
# Potentials: Forces, energies
#########################################

class Potential:
	"""
    No potential: ideal gas
    """
    #
	def __init__(self):
		self.cutoff = 999.e99
		self.latticeSpacing = 1.0
    #
		

	def Forces(self, atoms, neighborLocator, forcesSoFar=None):
		"""
		Base class allocates space for force arrays. If force array is 
		passed in (perhaps so that two potentials may be added), checks
		to make sure shape matches that of atoms.
		"""
		if forcesSoFar is None:
			forces = np.zeros(np.shape(atoms.positions),float)
		else:
			if np.shape(forcesSoFar)!=np.shape(atoms.positions):
				print("forcesSoFar shape is ", np.shape(forcesSoFar))
				print("atoms.positions shape is ", \
						np.shape(atoms.positions))
				print("forcesSoFar = ", forcesSoFar)
				print("atoms.positions =  ", atoms.positions)
			assert np.shape(forcesSoFar)==np.shape(atoms.positions), \
			"Wrong shape: forcesSoFar in Potential class"
			forces = forcesSoFar
		return forces
    #

	def PotentialEnergy(self, atoms, neighborLocator):
		return 0.0


class CompositePotential (Potential):
	"""
    Provides sum of a list of potentials
    """
    #
	def __init__(self, potentialList):
		"""Stores list of pointers to potentials"""
		Potential.__init__(self)
		self.potentials = potentialList[:]
    #
		

	def Forces(self, atoms, neighborLocator, forcesSoFar=None):
		"""
		Adds up forces from various potentials
		"""
		# XXX May want to associate neighborLocator with potential
		# since cutoff changes!
		neighborLocator.dist = self.potentials[0].cutoff
		forces = self.potentials[0].Forces(atoms, neighborLocator, forcesSoFar)
		for potential in self.potentials[1:]:
			neighborLocator.dist = potential.cutoff
			forces = potential.Forces(atoms, neighborLocator, forces)
		return forces
    #


	def PotentialEnergy(self, atoms, neighborLocator):
		"""
		Adds up energy from various potentials
		"""
		energy = 0.
		for potential in self.potentials:
			neighborLocator.dist = potential.cutoff
			energy += potential.PotentialEnergy(atoms, neighborLocator)
		return energy
	

class GravityPotential(Potential):
	"""
    Potential due to gravity
    """
    #
	def __init__(self, g=10.0, direction=1):
		"""Default: Acceleration in y-direction with g=10"""
		Potential.__init__(self)
		self.direction=direction
		self.g = g
		self.cutoff = 0.0
    #
		

	def Forces(self, atoms, neighborLocator, forcesSoFar=None):
		"""
		Add F = -m g to forces so far 
		(allocated if necessary by Potential base class)
		"""
		forces = Potential.Forces(self, atoms, neighborLocator, forcesSoFar)
		gravityAcceleration = np.zeros(dim)
		gravityAcceleration[self.direction]= -self.g
		forces = forces + atoms.mass * gravityAcceleration
		return forces
    #


	def PotentialEnergy(self, atoms, neighborLocator):
		"""
		Sums m g h over atoms
		"""
		energy = atoms.mass*self.g*np.sum(atoms.positions[:,self.direction])
		return energy


#class LennardJonesSimplePotential(Potential):
#    """
#    Lennard--Jones 6/12 potential
#    Slow version, not using NeighborLocator or array operations
#    """
#    #
#    def __init__(self, epsilon=1.0, sigma=1.0):
#	Potential.__init__(self)
#	self.epsilon = epsilon
#	self.sigma = sigma
#	# Lattice spacing depends slightly on second-neighbor force
#	self.latticeSpacing = 2.0**(1./6.)
#    #
#    def Forces(self, atoms, neighborLocator, boundaryConditions,
#		forcesSoFar=None):
#	"""
#	Add Lennard--Jones force to forces so far 
#	Pair force = epsilon * (48(sigma/r)^14 - 24*(sigma/r)^8)*r
#	"""
#	# Allocate forces (if necessary) in Potential base class
#	forces = Potential.Forces(self, atoms, neighborLocator,
#		boundaryConditions, forcesSoFar)
#	for n1, x1 in enumerate(atoms.positions):
#	    for n2, x2 in enumerate(atoms.positions[0:n1]): 
#		r = x2-x1
#		r = boundaryConditions.DifferenceBoundaryConditions([r])[0]
#		sOverR2 = self.sigma*self.sigma/np.dot(r,r)
#		sOverR4= sOverR2*sOverR2
#		sOverR8 = sOverR4*sOverR4
#		sOverR14 = sOverR2*sOverR4*sOverR8
#		F = self.epsilon*(48.* sOverR14 - 24.*sOverR8)*r
#		forces[n1] = forces[n1] - F
#		forces[n2] = forces[n2] + F
#	return forces
#    #
#    def PotentialEnergy(self, atoms, neighborLocator, boundaryConditions):
#	"""
#	Lennard Jones pair energy = 4*epsilon*((sigma/r)^12 - (sigma/r)^6)*r
#	"""
#	energy = 0
#	for n1, x1 in enumerate(atoms.positions):
#	    for n2, x2 in enumerate(atoms.positions[0:n1]): 
#		r = x2-x1
#		r = boundaryConditions.DifferenceBoundaryConditions([r])[0]
#		sOverR2 = self.sigma*self.sigma/np.dot(r,r)
#		sOverR4= sOverR2*sOverR2
#		sOverR6 = sOverR4*sOverR2
#		sOverR12 = sOverR6*sOverR6
#		energy += 4.0*self.epsilon*(sOverR12 - sOverR6)
#	return energy

class LennardJonesPotential(Potential):
	"""
    Lennard--Jones 6/12 potential
    Medium-fast version, no cutoff
    """
    #
	def __init__(self, epsilon=1.0, sigma=1.0):
		Potential.__init__(self)
		self.epsilon = epsilon
		self.sigma = sigma
		# Lattice spacing depends slightly on second-neighbor force
		self.latticeSpacing = 2.0**(1./6.)


    #
	def Forces(self, atoms, neighborLocator, forcesSoFar=None):
		"""
		Add Lennard--Jones force to forces so far 
		Pair force = epsilon * (48(sigma/r)^14 - 24*(sigma/r)^8)*r
		"""
		#
		# Allocate forces (if necessary) in Potential base class
		#
		forces = Potential.Forces(self, atoms, neighborLocator, forcesSoFar)
		#
		# Build lists with neighbor pairs, displacement, squared distance
		#
		n1, n2, r, r2 = neighborLocator.HalfNeighbors(atoms)
		#
		# Do things with vector operations: sOverRN = (sigma/r)^N
		#
		sOverR2 = self.sigma*self.sigma/r2
		sOverR4 = sOverR2*sOverR2
		sOverR8 = sOverR4*sOverR4
		sOverR14 = sOverR2*sOverR4*sOverR8
		fPrefactor = (48.*self.sigma)*sOverR14 - (24.*self.sigma)*sOverR8
		for i in range(len(n1)):
			f = fPrefactor[i]*r[i]
			forces[n1[i]] -= f
			forces[n2[i]] += f
		return forces
    #


	def PotentialEnergy(self, atoms=None, neighborLocator=None, r=None):
		"""
		Lennard Jones pair energy = 4*epsilon*((sigma/r)^12 - (sigma/r)^6)*r
		Can pass in distance r, vector of distances r, or atom list
		"""
		energy = 0.
		if r is not None:
			r2 = r*r
		else:
			#
			# Build lists with neighbor pairs, displacement, squared distance
			#
			n1, n2, dr, r2 = neighborLocator.HalfNeighbors(atoms)
		#
		# r2 can be a scalar or a whole vector of squared displacements
		# Do things with vector operations: sOverRN = (sigma/r)^N
		#

		sOverR2 = self.sigma*self.sigma/r2
		sOverR4 = sOverR2*sOverR2
		sOverR6 = sOverR2*sOverR4
		sOverR12 = sOverR6*sOverR6
		energy = np.sum(4.0*self.epsilon*(sOverR12-sOverR6))
		return energy


class LennardJonesCutPotential (Potential):
	"""
    Lennard--Jones 6/12 potential, cut off at r=2.7.
    Fast version. Two pieces: below cut1 energy is shifted by alpha,
	Pair Energy = 4 epsilon ( (sigma/r)^12-(sigma/r)^6 )+alpha
    between cut1 and cutoff is quadratic in r2:
	Pair Energy = beta + gamma (r/sigma)^2 + delta (r/sigma)^4
    Designed by JPS to keep two derivatives continuous and avoid square roots.
    """
    #
	def __init__(self, epsilon=1.0, sigma=1.0):
		Potential.__init__(self)
		self.epsilon = epsilon
		self.sigma = sigma
		self.cutoff = 2.7
		self.alpha = 0.01257815434637538131
		self.beta = -0.18714033296465953868
		self.gamma= 0.05134165513433732
		self.delta = -0.00352137552361713868
		self.cut1 = 2.41308778858241
		self.r2cut1 = self.cut1*self.cut1
		# Lattice spacing depends slightly on second-neighbor force
		if dim==2: self.latticeSpacing = 1.1132013
		elif dim==3: self.latticeSpacing = 1.095693678895
		else: self.latticeSpacing = 2.0**(1./6.)
    #
		

	def Forces(self, atoms, neighborLocator, forcesSoFar=None):
		"""
		Add Lennard--Jones force to forces so far 
		Inner pair force = epsilon * [48(sigma/r)^14 - 24*(sigma/r)^8] * r
		Outer pair force = epsilon * [-2 gamma - 4 delta (r/sigma)^2] * r
		"""
		#
		# Allocate forces (if necessary) in Potential base class
		#
		forces = Potential.Forces(self, atoms, neighborLocator, forcesSoFar)
		#
		# Build lists with neighbor pairs, displacement, squared distance
		#
		assert(neighborLocator.dist == self.cutoff)
		n1, n2, r, r2 = neighborLocator.HalfNeighbors(atoms)
		#
		# Do things with vector operations: sOverRN = (sigma/r)^N
		#
		rOverSigma2 = r2/(self.sigma*self.sigma)
		sOverR2 = self.sigma*self.sigma/r2
		sOverR4 = sOverR2*sOverR2
		sOverR8 = sOverR4*sOverR4
		sOverR14 = sOverR2*sOverR4*sOverR8
		fInnerPrefactor =  \
			(48.*self.epsilon)*sOverR14 - (24.*self.epsilon)*sOverR8
		fOuterPrefactor = -2.0*self.epsilon*self.gamma \
				-(4.0*self.epsilon*self.delta)*rOverSigma2
		for i in range(len(n1)):
			if r2[i] < self.r2cut1:
				f = fInnerPrefactor[i]*r[i]
			else:
				f = fOuterPrefactor[i]*r[i]
			forces[n1[i]] -= f
			forces[n2[i]] += f
		return forces
    #


	def PotentialEnergy(self, atoms=None, neighborLocator=None, r=None):
		"""
		Lennard Jones pair energy = 4*epsilon*((sigma/r)^12 - (sigma/r)^6)*r
		Can pass in distance r, vector of distances r, or atom list
		"""
		energy = 0.
		if r is not None:
			r2 = r*r
		else:
			assert(neighborLocator.dist == self.cutoff)
			#
			# Build lists with neighbor pairs, displacement, squared distance
			#
			n1, n2, dr, r2 = neighborLocator.HalfNeighbors(atoms)
		#
		# r2 can be a scalar or a whole vector of squared displacements
		# Do things with vector operations: sOverRN = (sigma/r)^N
		#
		rOverSigma2 = r2/(self.sigma*self.sigma)
		sOverR2 = self.sigma*self.sigma/r2
		sOverR6 = sOverR2*sOverR2*sOverR2
		sOverR12 = sOverR6*sOverR6
		energyInner = 4.0*self.epsilon*(sOverR12-sOverR6) + self.alpha
		energyOuter = self.epsilon * self.beta \
				+ self.epsilon * self.gamma * rOverSigma2 \
				+ self.epsilon * self.delta * rOverSigma2 * rOverSigma2
		inside = (r2<self.r2cut1)  # Array of ones and zeros based on truth
		energy = inside * energyInner + (1-inside)*energyOuter
		return energy


#########################################
# Movers
#########################################

def RunVelocityVerlet(atoms, observers, potential, neighborLocator, 
			boundaryConditions, nSteps=100, timeStep=0.01, 
			lastForces = None, initialTime = 0.0):
	"""
    Runs nSteps of velocity Verlet.
       velocities += 0.5 * force * dt
       positions += velocities * dt
       enforce boundary conditions
       re-calculate forces
       velocities += 0.5 * force * dt
       observe, display, diagnose
    Re-uses last force calculation if available. (Don't pass it in
    if anything relevant has been changed!)
    """
    #
    # If nothing has changed since last force calculation, re-use it.
    #
	if lastForces == None:
		forces = potential.Forces(atoms, neighborLocator)
	else:
		assert np.shape(lastForces) == np.shape(atoms.positions)
		forces = lastForces
    #
    # Verlet
    #
	time = initialTime
	for step in range(nSteps): 
		time += timeStep
		#
		# Accellerate ball half a timeStep
		#
		atoms.velocities += 0.5*forces*timeStep
		atoms.positions += atoms.velocities*timeStep
		#
		# Enforce boundary condition
		#
		boundaryConditions.EnforceBoundaryConditions(atoms) 
		#
		# Calculate new forces
		#
		forces = potential.Forces(atoms, neighborLocator)
		#
		# Accellerate ball half a timeStep
		#
		atoms.velocities += 0.5*forces*timeStep
		#
		# Observer
		#
		for observer in observers:
			observer.Update(atoms, time)
		#
		# Return forces for re-use
		#
	return forces, time

#########################################
# MDSystem: set up and manage runs for
#
# ListOfAtoms
# Initializers
# Transformers
# Observers
# NeighborLocators
# BoundaryConditions
# Potentials
# Movers
#
################################

class MDSystem:
	"""
    Base class for DigitalMaterial molecular dynamics systems
    """
	def __init__(self, L = 10., atoms = None, 
		observers = [], neighborLocator = None, boundaryConditions = None,
		potential = None, mover = None):
		#
		self.L = L
		#
		if atoms is not None:
			self.atoms = atoms
		else:
			self.atoms = ListOfAtoms()
		#
		if potential is not None:
			self.potential = potential
		else:
			self.potential = Potential()	# No potential default
		#
		if boundaryConditions is not None:
			self.boundaryConditions = boundaryConditions
		else:
			self.boundaryConditions = FreeBoundaryConditions(L)
		#
		if neighborLocator is not None:
			self.neighborLocator = neighborLocator
		else:
			self.neighborLocator = NoNeighborLocator(self.L,
						self.boundaryConditions)
		#
		self.observers = observers
		#
		if mover is not None:
			self.mover = mover
		else:
			self.mover = RunVelocityVerlet
    #
			
	
	def Run(self, nSteps=500, timeStep=0.01):
		"""
		Run an equilibrium distribution of a Lennard Jones gas in gravity;
		equilbrate to temprature T.
		Return final state of atoms.
		"""
		self.mover(self.atoms, self.observers, self.potential, 
				self.neighborLocator, self.boundaryConditions,
				nSteps, timeStep)
		
	
    #
	def Equilibrate(self, nCoolSteps=10, temperature=0.3,
		nStepsPerCool=50, timeStep=0.01):
		"""
		Tries to equilibrate system to given temperature 
		"""
		for step in range(nCoolSteps):
			self.Run(nSteps = nStepsPerCool)
			ThermalizingTransformer(self.atoms, temperature) 



# Copyright (C) Cornell University
# All rights reserved.
# Apache License, Version 2.0


