import numpy as np
import matplotlib.pyplot as plt


class grid:
    # Initialized two dimensional grid
    def __init__(self, size):
        self.size = size
        self.grid = np.zeros([self.size, self.size])



    def UpdateSteps(self, steps):
        """
            Update the steps not allowing walking back
            return the trajectory
        """

        # start the walks at the center of the grid
        position_list = [[] for i in range(steps)]
        x, y = int(self.size/2), int(self.size/2)
        position_list[0].append(x)
        position_list[0].append(y)


        self.grid[x][y] = 0
        past_step = [0, 0]
        for i in range(1, steps):
            # optinal steps to take in a 2d array
            options = [[0, 1], [1, 0], [-1, 0], [0, -1]]
            # self-avoiding random walks, no walking back
            if past_step in options:
                options.remove(past_step)
            step_index = np.random.choice(len(options))
            random_step = options[step_index]
            x += random_step[0]
            y += random_step[1]
            if [x, y] in position_list:
                return None
            self.grid[x][y]
            position_list[i].append(x)
            position_list[i].append(y)
            past_step = random_step
        self.pos = position_list[-1]
        return position_list[-1]


    def MeasureDistance(self):
        """
            Measure the distance away from the initial position
        """
        return np.sqrt((int(self.size/2) - self.pos[0])**2 + (int(self.size/2) - self.pos[1])**2)

    
