################################
# Pair distribution exercise
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
from importlib import reload
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("qtAgg")

reload(DM)

################################
# Dimension of space
################################

DM.dim = 2

################################
# Transfer other libraries from DM
################################


np = DM.np
vi = DM.vi 			# Visual Python
	# histograms without graphics



################################
# Set up pieces of simulation
################################


class PairDistribution (DM.MDSystem):
    """
    Appropriate potentials and stuff for measuring pair distribution function
    """
       
    def __init__(self, atoms=None, nAtoms=20, T=5.0, L=10.):
        potential = DM.LennardJonesCutPotential()
        boundaryConditions = DM.PeriodicBoundaryConditions(L)
        neighborLocator = DM.SimpleNeighborLocator(potential.cutoff,
							boundaryConditions)
        if atoms is None:
            atoms = DM.RandomNonoverlappingListOfAtoms(L, neighborLocator,
				 nAtoms=nAtoms, temperature=T)
        self.energyObserver = DM.EnergyObserver(potential, neighborLocator, 
                                                boundaryConditions)
        self.displayObserver = DM.VisualDisplayAtomsObserver(atoms,L)
        observers = [self.displayObserver, self.energyObserver]
        mover = DM.RunVelocityVerlet
        DM.MDSystem.__init__(self, L, atoms, observers, neighborLocator,
			     boundaryConditions, potential, mover)
        self.Reset()


    def PlotPairDistribution(self, nAverages=10, nStepsBetweenAverages=200,
				nBins=30, rMax=3., r=None, showPlot=True, savePlot=False):
        """
	    Plots histogram of pair distribution function
	    g(r) = (V/N) * (# atoms in [r, r+Delta r]) / (4 pi r^2 Delta r)
	    """
        if r is None:
            self.Reset()
            for av in range(nAverages):
                n1, n2, dr, r2 = self.neighborLocator.HalfNeighbors( \
							self.atoms,rMax)
                self.r.extend(np.sqrt(r2))
                self.Run()
            r = self.r
        n, bins,patchlist = plt.hist(r, nBins)
        rBar = 0.5*(bins[1:]+bins[:-1])
        deltaR = bins[1]-bins[0]
        nAtoms = len(self.atoms.positions)
        V = (self.boundaryConditions.L**DM.dim)
        if DM.dim==3:
            nOver4PiR2 = n[:] / (4.*np.pi*rBar*rBar)
            patches = plt.bar(bins[:-1], (2.*V/(nAtoms*nAtoms))* \
                    nOver4PiR2/(deltaR*nAverages), width=deltaR)
        elif DM.dim==2:
            nOver2PiR = n[:] / (2.*np.pi*rBar)
            patches = plt.bar(bins[:-1], (2.*V/(nAtoms*nAtoms))* \
                    nOver2PiR/(deltaR*nAverages), width=deltaR)
        plt.axis(xmin=0.0)
        if showPlot:
            plt.show()



    def Reset(self):
        self.r = []


    def Equilibrate(self, nCoolSteps=10, temperature=0.0, nStepsPerCool=100,
                timeStep=0.01):
        self.energyObserver.Reset()
        for step in range(nCoolSteps):
            DM.MDSystem.Equilibrate(self, nCoolSteps=1,
                temperature=temperature, nStepsPerCool=nStepsPerCool, 
                timeStep=timeStep)
            print("Thermalizing: kinetic temperature = ", \
                    self.energyObserver.Temperature())
            self.energyObserver.Reset()

def yesno():
    response = input('    Continue? (y/n) ')
    if len(response)==0:        # [CR] returns true
        return True
    elif response[0] == 'n' or response[0] == 'N':
        return False
    else:                       # Default
        return True
    
# def demo(nAverages=4, nCoolSteps=10):
def demo(nAverages=1, nCoolSteps=1):
    """Demonstrates solution for exercise: example of usage"""
    print("Pair distribution function Demo")
    #
    # Gas
    #
    T = 0.5
    L_Gas = 20.
    nGasAtoms=20
    sysGas = PairDistribution(L=L_Gas, nAtoms=nGasAtoms, T=T)
    print("Gas")
    rho=len(sysGas.atoms.positions)/(L_Gas*L_Gas)
    print("Gas Target temperature, actual density rho = ", T, rho)
    sysGas.Equilibrate(temperature=T, nStepsPerCool=1000, nCoolSteps=nCoolSteps)
    sysGas.Run()
    print("Gas final kinetic temperature=", \
		sysGas.energyObserver.Temperature())
    sysGas.Run()
    if not yesno(): return
    print("  Measuring pair distribution function: gas")
    plt.figure()
    sysGas.PlotPairDistribution(nAverages=nAverages, showPlot=False)
    r = np.arange(0.9,3.0,0.01)
    gTheory = np.exp(-sysGas.potential.PotentialEnergy(r=r)/T)
    plt.plot(r, gTheory, "r-", linewidth=2)
    plt.show()
 
    if not yesno(): return
    #
    # Liquid
    #
    # Find nice square for triangular lattice for crystal 
    # Continued fraction expansion for double layer
    # sqrt[3] = (1,1,2,1,2,1,2...) = 1/(1+2/(1+1/(...)))
    #
    crystalLatticeSpacing = DM.LennardJonesCutPotential().latticeSpacing
    #
    # Want liquid initial condition to have density 0.75:
    # crystal density as given is 0.96874
    #
    liquidLatticeSpacing = crystalLatticeSpacing * np.sqrt(0.96874/0.75)
    #
    # 6 layers upward = 3 Sqrt(3) = L = 5.19615 layers sideways 
    # (Alternative: 8 layers upward = 4 Sqrt(3) = L = 6.9282 layers sideways)
    #
    L_Crystal = 4.999*crystalLatticeSpacing
    L_Liquid = 4.999*liquidLatticeSpacing
    T = 0.5
    # XXX Atoms have to be set up at final temperature
    # XXX PlotPairDistribution does not thermalize them when passed in
    # XXX Need to pack atoms more tightly: fix box size?
    atoms = DM.TriangularSquareClusterListOfAtoms(L_Liquid, 
    		latticeSpacing=liquidLatticeSpacing, 
		temperature=T, origin = np.random.random(DM.dim))
    sysLiq = PairDistribution(atoms=atoms, L=L_Liquid)
    print("Liquid")
    rho=len(sysLiq.atoms.positions)/(L_Liquid*L_Liquid)
    print("Liquid Target temperature, actual density rho = ", T, rho)
    sysLiq.Equilibrate(temperature=T, nCoolSteps=nCoolSteps)
    sysLiq.Run()
    print("Liquid final kinetic temperature=", \
		sysLiq.energyObserver.Temperature())
    if not yesno(): return
    print("  Measuring pair distribution function: liquid")
    sysLiq.PlotPairDistribution(nAverages=nAverages, nBins=50)
    if not yesno(): return
    print("Crystal")
    T = 0.1
    atoms = DM.TriangularSquareClusterListOfAtoms(L_Crystal, 
    		latticeSpacing=crystalLatticeSpacing, 
		temperature=T, origin = np.random.random(DM.dim))
    sysSol = PairDistribution(atoms=atoms, L=L_Crystal)
    sysSol.Equilibrate(temperature=T, nCoolSteps=nCoolSteps)
    sysSol.Run()
    print("Solid final kinetic temperature=", \
		sysSol.energyObserver.Temperature())
    if not yesno(): return
    print("  Measuring pair distribution function: solid")
    sysSol.PlotPairDistribution(nAverages=nAverages, nBins=50)
    return sysGas, sysLiq, sysSol
   
if __name__=="__main__":
    demo()


# Copyright (C) Cornell University
# All rights reserved.
# Apache License, Version 2.0

