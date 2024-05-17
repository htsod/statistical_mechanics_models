import numpy as np
import matplotlib.pyplot as plt


class grid:
    def __init__(self, size):
        self.size = size
        self.grid = np.zeros([self.size, self.size])



    def UpdateSteps(self, steps):
        position_list = [[] for i in range(steps)]
        x, y = int(self.size/2), int(self.size/2)
        position_list[0].append(x)
        position_list[0].append(y)

        self.grid[x][y] = 0
        past_step = [0, 0]
        for i in range(1, steps):
            options = [[0, 1], [1, 0], [-1, 0], [0, -1]]
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
        return np.sqrt((int(self.size/2) - self.pos[0])**2 + (int(self.size/2) - self.pos[1])**2)

        


    
grid_points = grid(1000)

steps, repeatition = 10, 1000
rms_distance = []
for s in range(1, steps):
    r_square = []
    for rep in range(repeatition):
        last_position = grid_points.UpdateSteps(s)
        if last_position is not None:
            d = grid_points.MeasureDistance()
            r_square.append(d**2)
    rms_distance.append(np.sqrt(np.array(r_square).mean()))


fig, ax = plt.subplots()

ax.scatter(np.arange(1, steps), rms_distance, label="random walks")
ax.plot(np.arange(1, steps), np.arange(1, steps) - 1, label="R = L")
ax.plot(np.arange(1, steps), np.arange(1, steps)**(3/4) - 1, label=r"R = $L^{3/4}$")

ax.legend()
ax.set_title("RMS distance versus Number of steps")
ax.set_xlabel("Number of steps")
ax.set_ylabel("RMS distance")
plt.show()