import numpy as np
import matplotlib.pyplot as plt


class RandomWalker:
    def __init__(self, dim, step=1/2):
        """
        Define a walker object
        """
        self.dim = dim
        self.step = step
        self.pos = np.array([0 for i in range(self.dim)])
        
    
    def UpdatePosition(self):
        """
        Update the position of the walker object
        """
        self.pos = self.pos + np.array([np.random.uniform(low=-self.step, high=self.step) for i in range(self.dim)])
    

def CreateWalkers(N, dim, step=1/2):
    WalksList = []
    for n in range(N):
        WalksList.append(RandomWalker(dim, step))
    return WalksList


def EvolveWalkers(WalksList, TotalSteps):
    trajectory = [[] for i in range(len(WalksList))]
    for s in range(TotalSteps):
        for n in range(len(WalksList)):
            WalksList[n].UpdatePosition()
            trajectory[n].append(WalksList[n].pos)
    return trajectory


def ListOfEndPoints(WalksList, TotalSteps):
    end_points = []
    for s in range(TotalSteps):
        for n in range(len(WalksList)):
            WalksList[n].UpdatePosition()
    for w in WalksList:
        end_points.append(w.pos)
    return np.array(end_points)


def gaussian(x, sigma):
    denominator = np.sqrt(2*np.pi)*sigma
    factor = np.exp(-x**(2)/(2*sigma**(2)))
    return factor / denominator


def plot1D(N, steps):
    WalksList = CreateWalkers(N, 1)
    trajectory = EvolveWalkers(WalksList, steps)

    fig, ax = plt.subplots()
    for p in range(len(trajectory)):
        ax.plot(trajectory[p])
    ax.set_title(f"{N} one-dimensional Random Walks with {steps} steps")
    ax.set_ylabel("Distance")
    ax.set_xlabel("Steps")
    plt.show()


def plot2D(N, steps):
    WalksList = CreateWalkers(N, 2)
    trajectory = EvolveWalkers(WalksList, steps)

    fig, ax = plt.subplots()
    for p in range(len(trajectory)):
        ax.plot((np.array(trajectory[p]).T)[0], (np.array(trajectory[p]).T)[1])
    ax.set_title(f"{N} two-dimensional Random Walks with {steps} steps")
    ax.set_ylabel("y")
    ax.set_xlabel("x")
    plt.show()


def plot2DScatter(N1, N2, s1, s2):
    WalksList = CreateWalkers(N2, 2)
    trajectory = EvolveWalkers(WalksList, s2)
    fig, ax = plt.subplots()
    for w in range(len(trajectory)):
        ax.scatter(trajectory[w][-1][0], trajectory[w][-1][1], color="b")
    
    WalksList = CreateWalkers(N1, 2)
    trajectory = EvolveWalkers(WalksList, s1)
    for w in range(len(trajectory)):
        ax.scatter(trajectory[w][-1][0], trajectory[w][-1][1], color="orange")

    ax.set_title(f"Scatter plot for {s1} steps and {s2} steps for {N1} and {N2} random walks")
    ax.set_ylabel("y")
    ax.set_xlabel("x")
    plt.show()


def plot2DScatter(N1, N2, s1, s2):
    WalksList = CreateWalkers(N2, 2)
    end_points = ListOfEndPoints(WalksList, TotalSteps=s2)
    fig, ax = plt.subplots()
    ax.scatter(end_points.T[0], end_points.T[1], color="b")
    
    WalksList = CreateWalkers(N1, 2)
    end_points = ListOfEndPoints(WalksList, TotalSteps=s1)
    ax.scatter(end_points.T[0], end_points.T[1], color="orange")

    ax.set_title(f"Scatter plot for {s1} steps and {s2} steps for {N1} and {N2} random walks")
    ax.set_ylabel("y")
    ax.set_xlabel("x")
    plt.show()


def plotHist(N, s):
    WalksList = CreateWalkers(N, 1)
    end_points = ListOfEndPoints(WalksList, s)
    fig, ax = plt.subplots()
    ax.hist(end_points[:], bins=50, density=True, label="RandomWalks")

    a = 0.289
    sigma = np.sqrt(s) * a
    x = np.linspace(-3*sigma, 3*sigma, N)
    y = gaussian(x, sigma)
    ax.plot(x, y, label="Gaussian")

    ax.set_title(f"one-dimensional randomwalk for {s} steps")
    plt.show()




