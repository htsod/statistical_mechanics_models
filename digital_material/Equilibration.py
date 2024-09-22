################################
# Equilibration exercise
################################

################################
# Import Digital Material 
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

import DigitalMaterial as DM
import importlib as imp
imp.reload(DM)

################################
# Dimension of space
################################

DM.dim = 3

################################
# Transfer other libraries from DM
################################
np = DM.np
vi = DM.vi 			# Visual Python
pylab = DM.pylab

################################
# Set up pieces of simulation
################################


class EquilibrationTest(DM.MDSystem):
	"""
    Appropriate potentials and stuff for watching random walks of gas molecules
    """
    #
	def __init__(self, atoms=None, nAtoms=10, T=1.0, L=4., minDist=1.0):
		potential = DM.LennardJonesCutPotential()
		boundaryConditions = DM.PeriodicBoundaryConditions(L)
		neighborLocator = DM.SimpleNeighborLocator(potential.cutoff, boundaryConditions)
		if atoms is None:
			atoms = DM.RandomNonoverlappingListOfAtoms(L, neighborLocator, minDist=1.0, nAtoms=nAtoms, temperature=T)
		# Remove center-of-mass momentum from gas
		atoms.velocities -= DM.scipy.sum(atoms.velocities)/len(atoms.velocities)
		self.displayObserver = DM.VisualDisplayAtomsObserver(atoms,L)
		self.velocityTrajectoryObserver = DM.VelocityTrajectoryObserver()
		observers = [self.displayObserver, self.velocityTrajectoryObserver]
		mover = DM.RunVelocityVerlet
		DM.MDSystem.__init__(self, L, atoms, observers, neighborLocator,
					boundaryConditions, potential, mover)
    #
	def PlotHistograms(self, nHists=6, factorBetweenHists=4):
		"""
		Plots velocity histograms for nHists time intervals growing by factor
		(Equals momentum histogram if m=1)
		"""
		for h in range(nHists):
			nSteps = factorBetweenHists**h
			print(nSteps)
			self.Run(nSteps=nSteps)
			pylab.figure(h+1)
			velocity_component = [x
						 for xxs in self.velocityTrajectoryObserver.vTrajectory
						 for xs in xxs
						 for x in xs]
			pylab.hist(velocity_component, density=True, bins=50)
			pylab.show()
			self.velocityTrajectoryObserver.Reset()


def demo():
    """Demonstrates solution for exercise: example of usage"""
    print("Equilibration demo")
    sys = EquilibrationTest()
    sys.PlotHistograms()
    

if __name__=="__main__":
    demo()


# Copyright (C) Cornell University
# All rights reserved.
# Apache License, Version 2.0


