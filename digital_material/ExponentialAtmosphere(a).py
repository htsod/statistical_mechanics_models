################################
# Exponential Atmosphere exercise
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
scipy = DM.scipy
RandomArray = DM.RandomArray
vi = DM.vi 			# Visual Python
pylab = DM.pylab

################################
# Set up pieces of simulation
################################

def ExponentialAtmosphereListOfAtoms(g=10.0, L=20, nAtoms=1000, 
			temperature=30.0,
                        mass=1.0, radius=0.5, color=vi.color.green, 
			heightIndex=1, dimension = DM.dim):
    """
    Put atoms of radius r in box (0,L)^dim, with probability
    decaying like Exp(-m g z / kB T)
    """
    nAtomsToTry = int(2.0*nAtoms/(1.0-np.exp(-mass*g*L/temperature)) + 20)
    #
    # Build lots of extra atoms, prune those with z>L
    #
    positionsT = []
    min = radius		# Don't overlap with walls
    max = L-radius
    for i in range(dimension):
        if (i != heightIndex):
            positionsT.append(RandomArray.uniform(min, max, nAtoms))
        else:
            atomsToTry=int(2*nAtoms/(1-np.exp(-mass*g*L/temperature)))+20
        
            if atomsToTry < 1000000:	# Truncate exponential distribution
    		#
		# Exponential distribution, with minimum
		# height equal to atoms.radius (atoms repelled from floor)
    		#
                posToTry = RandomArray.exponential(
                    temperature/(mass*g),nAtomsToTry)+radius
                posOK = [pos for i, pos in enumerate(posToTry) if pos<max]
                #
                # If not enough atoms, repeat
                #
                while len(posOK)<nAtoms:
                    posToTry = RandomArray.exponential(
                        temperature/(mass*g),nAtomsToTry)+radius
                    posOK = [pos for i, pos in enumerate(posToTry) if pos<max]
            else:	# Prune uniform distribution
                atomsToTry = int(2*nAtoms/np.exp(-mass*g*L/temperature))+20
                posToTry = RandomArray.uniform(min, max, nAtomsToTry) 
                posOK = [pos for i, pos in enumerate(posToTry) 
                        if RandomArray.random()<np.exp(mass*g*pos/temperature)]
                while len(posOK)<nAtoms:
                    posToTry = RandomArray.uniform(min, max,nAtomsToTry)
                    posOK = [pos for i, pos in enumerate(posToTry) 
                            if RandomArray.random() < 
                            np.exp(mass*g*pos/temperature)]
            positionsT.append(np.array(posOK[0:nAtoms]))
    positions = np.transpose(np.array(positionsT))
    velocities = np.zeros(np.shape(positions))
    atoms = DM.ListOfAtoms(mass, radius, color, positions, velocities)
    DM.ThermalizingTransformer(atoms, temperature)
    return atoms


def IdealGasInGravity(atoms=None, nAtoms=1000, T=30., L=20., g=10.0,
		      nSteps=500, timeStep=0.01):
    """
    Run an equilibrium distribution of an ideal gas in gravity. 
    Return final configuration of atoms.
    Rerun with IdealGasInGravity(atoms)
    """
    potential = DM.GravityPotential(g=g)
    boundaryConditions = DM.ReflectiveBoundaryConditions(L)
    neighborLocator = DM.NoNeighborLocator(potential.cutoff, boundaryConditions)
    if atoms is None:
        atoms = ExponentialAtmosphereListOfAtoms(
			L=L,g=g,nAtoms=nAtoms,temperature=T)

    displayObserver = DM.VisualDisplayAtomsObserver(atoms,
    				boundaryConditions.L, velocityColors=True)
    observers = [displayObserver]

    mover = DM.RunVelocityVerlet

    mover(atoms, observers, potential, neighborLocator, boundaryConditions,
    			nSteps, timeStep)
    return atoms
 

def IdealGasDistributions(atoms, g, T, L):
    import pylab
    pylab.figure(1)
    # Atoms repelled from floor by distance equal to radius
    pylab.hist([pos[1]-atoms.radius for pos in atoms.positions], 
			bins=50, density=True)
    zs = np.arange(0,L-2.*atoms.radius)	
    rho =(atoms.mass*g/T) * np.exp(-atoms.mass*g*(zs)/T)
    pylab.plot(zs, rho, "k-", linewidth=2)
    pylab.title('Height distribution')
    pylab.xlabel('Height')
    pylab.ylabel('Probability density rho(v)')
    pylab.figure(2)
    pylab.hist(atoms.velocities, bins=50, density=True)
    sigmaV=np.sqrt(T/atoms.mass)
    vs = np.arange(-3.0*sigmaV, 3.0*sigmaV)
    rho = (1.0/(np.sqrt(2.0*np.pi)*sigmaV)) \
		* np.exp(-0.5*atoms.mass*vs*vs/T)
    pylab.plot(vs, rho, "k-", linewidth=2)
    pylab.title('Velocity distribution: all velocities')
    pylab.xlabel('Velocity v')
    pylab.ylabel('Probability density rho(v)')
    pylab.figure(3)
    upperHalfVelocities = [vel for i, vel in enumerate(atoms.velocities)
				if atoms.positions[i][1]>L/2.]
    pylab.hist(upperHalfVelocities, bins=10, density=True)
    pylab.plot(vs, rho, "k-", linewidth=2)
    pylab.title('Top half of box: Velocity distribution')
    pylab.xlabel('Velocity v')
    pylab.ylabel('Probability density rho(v)')
    pylab.show()


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
    print("Exponential Atmosphere Demo: Ideal gas")
    print("  Ideal gas under gravity")
    atoms = IdealGasInGravity(nAtoms=500, T=30., L=20., g=10.0, atoms=None,
		      nSteps=500, timeStep=0.01)
    if not yesno(): return
    print("  Position and velocity distributions")
    IdealGasDistributions(atoms, g=10, T=30., L=20.)


if __name__=="__main__":
    demo()


# Copyright (C) Cornell University
# All rights reserved.
# Apache License, Version 2.0

