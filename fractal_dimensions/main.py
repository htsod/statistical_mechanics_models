import numpy as np
import matplotlib.pyplot as plt


def LogisticMapping(x, mu):
    return 4*mu * x * (1-x)


def main(mu, epsilon_list):
    """
    Collect points about the attractor at mu
    Return the set of vectors P_{j} of bin size epsilon
    """
    x = np.random.random()
    while abs(x-0.5) <= 0.01:
        x = np.random.random()
    print(x)
        
    
    # repeat mapping until getting to the attractor
    rep = 1000
    x_list = []
    for r in range(rep):
        x = LogisticMapping(x, mu)
    
    
    # Now collect the data points near the attractor
    N_tot = 2**15
    x_list = [x]
    for n in range(N_tot):
        x = LogisticMapping(x, mu)
        x_list.append(LogisticMapping(x, mu))

    

    count_list = []
    bins_list = []
    for e in epsilon_list:
        # counts, edges, bars = plt.hist(x_list, bins=int(1/e), density=True)
        # p_list.append(counts)
        counts, bins = np.histogram(x_list, bins=int(1/e), density=True)
        count_list.append(counts)
        bins_list.append(bins)
    return count_list, bins_list


def D_capacity(N_e, e_list):
    d_cap = []
    x_axis = []
    for i in range(len(N_e)-1):
        d_cap.append((N_e[i+1] - N_e[i]) / (e_list[i+1] - e_list[i]))
        x_axis.append((np.log(e_list[i+1]) + np.log(e_list[i])/2))
    return d_cap, x_axis


mu = 0.8
e_list = np.linspace(0.001, 0.01, 10)

count_list, bins_list = main(mu, e_list)

N_e = []
for l in count_list:
    n = np.count_nonzero(np.array(l))
    N_e.append(n)

# print(N_e)
# print(e_list)

d_cap, x_axis = D_capacity(N_e, e_list)

# print(d_cap, x_axis)

plt.plot(x_axis, d_cap)
plt.show()



# for i in range(len(count_list)):
#     plt.stairs(count_list[i], bins_list[i])
#     plt.show()

