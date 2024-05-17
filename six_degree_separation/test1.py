import network_model
import matplotlib.pyplot as plt

L = 1000
Z = 2
p = 0.1
net = network_model.network(L, Z, p)
net.RandomAdd()



net.Display()

dist_list = net.FindAllPathLengths()

fig, ax = plt.subplots()
ax.hist(dist_list, bins=20, density=True)
ax.set_title(f"Separation Distribution for L = {L}, Z = {Z}, p = {p}")
ax.set_xlabel("Distance of Separation")
ax.set_ylabel("Fractional Counts")
plt.show()