import network_model
import matplotlib.pyplot as plt
import csv


L = 100
Z = 2
p = 0.2

N = int(input("Please enter number of iteration "))

mean_path_list = []
dist_list = []

for i in range(N):
    net = network_model.network(L, Z, p)
    net.RandomAdd()
    dist, mean_path = net.FindAveragePathLength()
    mean_path_list.append(mean_path)
    dist_list.append(dist)


with open('mean_path_list.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    spamwriter.writerow(mean_path_list)


fig, ax = plt.subplots(N)

for i in range(N):
    ax[i].hist(dist_list[i])
    ax[i].set_title("Separation Distribution")
    ax[i].set_xlabel("Distance of Separation")
    ax[i].set_ylabel("Counts")

plt.tight_layout()
plt.show()


