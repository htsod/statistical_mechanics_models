# Sometimes gives interactive new windows
# Must show() after plot, figure() before new plot
# %matplotlib

# Adds static figures to notebook: good for printing
# %matplotlib inline 

# Interactive windows inside notebook! Must include plt.figure() between plots


# Better than from numpy import *, but need np.sin(), np.array(), plt.plot(), etc.
import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import odeint



EarthMass = 1.;
TrueJupiterMass = 317.83;
SunMass = 332830.
RJ = 5.20316
G = 4*np.pi**2 / SunMass


def SetOmega(JupiterMass):
    """
    As Jupiter's mass is varied, its orbital angular frequency changes
    """
    return np.sqrt(G*(JupiterMass+SunMass) / RJ**3)

def SetJupiterMass(JupiterMass):
    """
    Usage: 
    JupiterMass, omega, RJcm, RScm, GMeMj, GMeMs = SetJupiterMass(newJupiterMass)
    
    omega is Jupiter's orbital frequency (which in the original app wasn't changed with the mass)
    RJcm is the circular orbit radius Jupiter makes about the center of mass
    RScm is the same for the Sun
    GMeMj and GMeMs are just convenient for Newton's laws
    """
    omega = SetOmega(JupiterMass)
    RJcm = RJ * (SunMass/(JupiterMass+SunMass))
    RScm = RJ * (JupiterMass/(JupiterMass+SunMass))
    GMeMj = G*EarthMass*JupiterMass
    GMeMs = G*EarthMass*SunMass
    return JupiterMass, omega, RJcm, RScm, GMeMj, GMeMs


RxEarth0 = 1.
RyEarth0 = 0.
VxEarth0 = 0.
VyEarth0 = 2*np.pi
t0 = 0.



def Newton(RV,t,Mj):
    """
    Given time and RV = (Rx, Ry, Vx, Vy) for Earth, returns derivatives of RV wrt time
    """
    JupiterMass, omega, RJcm, RScm, GMeMj, GMeMs = SetJupiterMass(Mj)
    Rx, Ry, Vx, Vy = RV
    # Solar and Jovian positions at time t
    RxJupiter = RJcm * np.cos(omega * t);
    RyJupiter = RJcm * np.sin(omega * t);
    RxSun = - RScm * np.cos(omega * t);
    RySun = - RScm * np.sin(omega * t);
    # vx, vy derivatives are F/M
    REJ = np.sqrt( (RxJupiter-Rx)**2 + (RyJupiter-Ry)**2 );
    RSJ = np.sqrt( (RxSun-Rx)**2 + (RySun-Ry)**2 );
    Fx = GMeMj*(RxJupiter-Rx) / REJ**3 + GMeMs*(RxSun-Rx)/RSJ**3;
    Fy = GMeMj*(RyJupiter-Ry)/REJ**3 + GMeMs*(RySun-Ry)/RSJ**3;
    return np.array([Vx, Vy, Fx/EarthMass, Fy/EarthMass])

def PlotTrajectory(Mj,tMax,dt=0.01):
    plt.figure()
    times = np.arange(0.,tMax,dt)
    Rx, Ry, Vx, Vy = odeint(Newton, [RxEarth0, RyEarth0, VxEarth0, VyEarth0],\
                            times, args=(Mj,)).transpose()
    plt.axes().set_aspect('equal')
    plt.plot(Rx,Ry);
    plt.show()
    




def EnergySurfaceInitial(r_radial, v_radial, Current_Energy, Mj):
    """
        Given Jupiter, Earth, and Sun are in a line with the Earth between Jupiter and the Sun, find initial conditions
        that will give a trajectory at a given rotating reference frame energy.
        
        Assumes new trajectory is launched with the three bodies along the x-axis, so solves for Vy given Rx, Vx.
        
        Doesn't do anything if v_radial too large for solution (too much radial kinetic and potential energy already)
        prints 'verboten' if so.
    """
    JupiterMass, omega, RJcm, RScm, GMeMj, GMeMs = SetJupiterMass(Mj)
    Rx = r_radial
    Vx = v_radial
    Ry = 0.
    RxJupiter = RJcm
    RxSun = -RScm
    Potential_Energy = -GMeMj/abs(RxJupiter-Rx) -GMeMs/abs(RxSun-Rx)
    C = -Current_Energy + Potential_Energy + 0.5 * EarthMass*Vx**2
    discriminant = omega**2 * Rx**2 - 2*C/EarthMass
    if discriminant > 0:
        Vy = omega * Rx + np.sqrt(discriminant)
        return np.array([Rx, Ry, Vx, Vy])
    else:
        print('verboten')
        return None


from scipy.interpolate import InterpolatedUnivariateSpline

def PoincareIntersections(Rx,Ry,Vx,Vy,times,Mj):
    """
    Returns coordinates where Earth passes through Jupiter's line through the origin, 
    with positive y velocity in rotating frame
    """
    omega = SetOmega(Mj)
    SJ = np.sin(omega*times) # Note these operations work on the whole array of times at once
    CJ = np.cos(omega*times)
    check = Rx*SJ - Ry*CJ   # Zero if Earth and Jupiter are in same line with the origin (and hence with the Sun)
    checkI = InterpolatedUnivariateSpline(times, check)
    tZeroBothDirections = checkI.roots() # Finds interpolated times where they cross
    # direction = True for those where Earth passes with positive velocity
    direction = (checkI.derivative()(tZeroBothDirections) < 0.)     
    tZero = np.array([tZ for tZ, d in zip(tZeroBothDirections, direction) if d])
    # Interpolate four coordinates. (Could interpolate r_radial and v_radial instead.)
    RxI = InterpolatedUnivariateSpline(times, Rx)
    RyI = InterpolatedUnivariateSpline(times, Ry)
    VxI = InterpolatedUnivariateSpline(times, Vx)
    VyI = InterpolatedUnivariateSpline(times, Vy)
    r_radial = RxI(tZero)*np.cos(omega*tZero) + RyI(tZero)*np.sin(omega*tZero)
    v_radial = VxI(tZero)*np.cos(omega*tZero) + VyI(tZero)*np.sin(omega*tZero)
    
    return r_radial, v_radial


def RotatingFrameEnergy(RV, t, Mj):
    """
    Given time and RV = (Rx, Ry, Vx, Vy) for Earth, returns the conserved 'rotating reference frame' energy
    """
    JupiterMass, omega, RJcm, RScm, GMeMj, GMeMs = SetJupiterMass(Mj)
    Rx, Ry, Vx, Vy = RV
    # Kinetic Energy in rotating frame
    KineticEnergy = 0.5 * EarthMass * (Vx**2+Vy**2) + EarthMass*omega*(Vx*Ry-Vy*Rx)
    # Solar and Jovian positions at time t
    RxJupiter = RJcm * np.cos(omega * t);
    RyJupiter = RJcm * np.sin(omega * t);
    RxSun = - RScm * np.cos(omega * t);
    RySun = - RScm * np.sin(omega * t);
    REJ = np.sqrt( (RxJupiter-Rx)**2  + (RyJupiter-Ry)**2 );
    RES = np.sqrt( (RxSun-Rx)**2 + (RySun-Ry)**2 );
    PotentialEnergy = -GMeMj/REJ - GMeMs/RES
    return PotentialEnergy+KineticEnergy
    


class PoincareLauncher:
    def __init__(self, Mj, tMax=50.):
        self.Mj = Mj
        self.tMax = tMax
        omega = SetOmega(Mj)
        self.omega = omega
        self.Current_Energy = RotatingFrameEnergy([RxEarth0,RyEarth0,VxEarth0,VyEarth0],\
                                                  t0,Mj)
        self.times = np.arange(t0,tMax,0.02)
        Rx, Ry, Vx, Vy = odeint(Newton,[RxEarth0, RyEarth0, VxEarth0, VyEarth0],\
                                self.times,args=(self.Mj,)).transpose()
        r, v = PoincareIntersections(Rx,Ry,Vx,Vy,self.times,self.Mj)
        fig = plt.figure("Poincaré Section, Mj=%s"%self.Mj)
        plt.plot(r,v,'.')
        plt.show()
        self.cid = fig.canvas.mpl_connect('button_press_event', self);

    def __call__(self, event):
        Rx = event.xdata
        Vx = event.ydata
        initial = EnergySurfaceInitial(Rx, Vx, self.Current_Energy, self.Mj)
        if initial is not None:
            Rx,Ry,Vx,Vy = odeint(Newton,initial,self.times,args=(self.Mj,)).transpose()
            r, v = PoincareIntersections(Rx,Ry,Vx,Vy, self.times, self.Mj)
            plt.figure("Poincaré Section, Mj=%s"%self.Mj)
            plt.plot(r,v,'.')
            plt.draw();



launcher = PoincareLauncher(Mj=22000.,tMax=1000.)