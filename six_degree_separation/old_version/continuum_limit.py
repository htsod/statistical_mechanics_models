import network_model
import numpy as np
import matplotlib.pyplot as plt
import csv


L = 40
Z = 2

N = int(input("Please enter number of iteration "))
p = np.linspace(0.0001, 1, N)

mean_path_list = []

net = network_model.network(L, Z, p[0])
net.RandomAdd()
_, mean_path = net.FindAveragePathLength()
norm_mean_path = mean_path
mean_path_list.append(mean_path / norm_mean_path)

for i in range(1, len(p)):
    print(f"p is equal to {p[i]}")
    net = network_model.network(L, Z, p[i])
    net.RandomAdd()
    _, mean_path = net.FindAveragePathLength()
    mean_path_list.append(mean_path / norm_mean_path)



# with open('variation_p.csv', 'w', newline='') as csvfile:
#     spamwriter = csv.writer(csvfile, delimiter=' ',
#                             quotechar='|', quoting=csv.QUOTE_MINIMAL)
#     spamwriter.writerow(mean_path_list)


fig, ax = plt.subplots()

ax.scatter(p, mean_path_list)
ax.set_xscale("log")
ax.set_title("Mean Path Variation")
ax.set_ylabel("Mean Path / Steps")
ax.set_xlabel("p")

plt.show()


