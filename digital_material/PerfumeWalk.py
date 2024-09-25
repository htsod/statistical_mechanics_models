################################
# Perfume Walk exercise
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
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("qtAgg")
imp.reload(DM)

################################
# Dimension of space
################################

DM.dim = 3

################################
# Transfer other libraries from DM
################################

vi = DM.vi 			# Visual Python
np = DM.np

################################
# Set up pieces of simulation
################################


class PerfumeWalk(DM.MDSystem):
    """
    Appropriate potentials and stuff for watching random walks of gas molecules
    """
    #
    def __init__(self, atoms=None, nAtoms=10, T=10.0, L=4., minDist=1.0):
        potential = DM.LennardJonesCutPotential()
        boundaryConditions = DM.PeriodicBoundaryConditions(L)
        neighborLocator = DM.SimpleNeighborLocator(potential.cutoff, boundaryConditions)
        if atoms is None:
            atoms = DM.RandomNonoverlappingListOfAtoms(L, neighborLocator, 
                                                       minDist=1.0, nAtoms=nAtoms,
                                                       temperature = T)
        # Remove center-of-mass momentum from gas
        atoms.velocities -= DM.np.sum(atoms.velocities)/len(atoms.velocities)
        self.displayObserver = DM.VisualDisplayAtomsObserver(atoms, L)
        self.displayObserver.ball_list[0].color=DM.vi.color.yellow
        self.unfoldedTrajectoryObserver = DM.UnfoldedTrajectoryObserver(L)
        observers = [self.displayObserver, self.unfoldedTrajectoryObserver]
        mover = DM.RunVelocityVerlet
        DM.MDSystem.__init__(self, L, atoms, observers, neighborLocator,
                    boundaryConditions, potential, mover)
    #


    def PlotRandomWalks(self):
        for at in range(len(self.atoms.positions)):
            x, y = self.unfoldedTrajectoryObserver.XY(at)
            plt.plot(x-x[0],y-y[0])
        plt.show()

        for at in range(len(self.atoms.positions)):
            x, y = self.unfoldedTrajectoryObserver.XY(at)
            plt.plot(np.arange(len(x)), x-x[0])
        plt.show()
            

    def PlotR2Bar(self):
        r2 = self.unfoldedTrajectoryObserver.r2Bar()
        plt.plot(r2)
        plt.show()


def RunCrystal(T=0.5):
    pot = DM.LennardJonesCutPotential()
    latticeSpacing = pot.latticeSpacing
    L=4.
    atoms = DM.FCCSphericalClusterListOfAtoms(R=L/2., 
    			latticeSpacing = latticeSpacing,
			center = vi.vector(L/2., L/2., L/2.), temperature=T)
    sys = PerfumeWalk(atoms, T=T, L=L)
    sys.Run()
    return sys

def yesno():
    response = input('    Continue? (y/n) ')
    if len(response)==0:        # [CR] returns true
        return True
    elif response[0] == 'n' or response[0] == 'N':
        return False
    else:                       # Default
        return True
    

def demo():
    """Demonstrates solution for exercise: example of usage"""
    print("Perfume Walk Demo: Gas State")
    sys = PerfumeWalk()
    sys.Run(2000)
    if not yesno(): return
    print("  Random Gas Atom Walks")
    sys.PlotRandomWalks()
    if not yesno(): return
    print("R^2 in gas versus number of timesteps (Change to time?)")
    sys.PlotR2Bar()
    if not yesno(): return
    print("Crystal of atoms, slowly heated")
    L = 4.
    T = 0.5
    atoms = DM.FCCSphericalClusterListOfAtoms(R=L/2., 
    			latticeSpacing = sys.potential.latticeSpacing,
			center = [L/2., L/2., L/2.], temperature = T)
    sys = PerfumeWalk(atoms, T=T, L=4)
    sys.Run()
    if not yesno(): return
    print("  Random Gas Atom Walks")
    sys.PlotRandomWalks()
    

if __name__=="__main__":
    demo()


# Copyright (C) Cornell University
# All rights reserved.
# Apache License, Version 2.0


