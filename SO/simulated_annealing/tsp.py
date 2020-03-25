import numpy as np
import matplotlib.pyplot as plt


def calc_energy(state, cities):
    '''Calculate energy as sum of distances between cities'''
    e = 0
    for i in range(len(state)):
        e += np.linalg.norm(cities[state[i]] - cities[state[i-1]], ord=2)
    return e


def is_trans(p):
    return np.random.rand() <= p


def trans_prob(dE, t):
    return np.exp(-dE/t)


def make_trans(path):
    path = path.copy()
    i, j = np.random.randint(len(path), size=2)
    if i < j:
        path[i:j] = np.flipud(path[i:j])
    else:
        path[j:i] = np.flipud(path[j:i])
    #path[i], path[j] = path[j], path[i]
    return path


def annealing(cities, init_temp, min_temp, decay=1.2, max_iter=1e5):
    n = len(cities)
    # initilization
    state = np.random.randint(n, size=n)
    T = init_temp
    E = calc_energy(state, cities)
    i = 1
    while T > min_temp:
        new_state = make_trans(state)
        new_energy = calc_energy(new_state, cities)
        dE = new_energy - E
        if dE < 0:
            E = new_energy
            state = new_state
        else:
            if is_trans(trans_prob(dE, T)):
                state = new_state
                E = new_energy
        if i == max_iter:
            break
        T -= T/(1+i*decay)
        i += 1
    return state


# 63.7 (n_cities=100, n_iter=1e5)
def main():
    t_min = 1e-3
    init_t = 10
    max_axis = 10
    n_cities = int(input())
    np.random.seed(42)
    cities = np.random.rand(n_cities, 2)*max_axis
    #plt.scatter(cities[:, 0], cities[:, 1])
    # plt.show()
    optimal_path = annealing(cities, init_t, t_min)
    print(calc_energy(optimal_path, cities))
    plt.plot(cities[optimal_path][:, 0], cities[optimal_path][:, 1])
    plt.show()


if __name__ == "__main__":
    main()
