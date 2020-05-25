import numpy as np
import matplotlib.pyplot as plt


def calc_energy(state, cities):
    """Calculate energy as sum of distances between cities."""
    return np.linalg.norm(cities[state[1:]] - cities[state[:-1]], ord=2) \
            + np.linalg.norm(cities[state[0]] - cities[state[-1]], ord=2) 


def is_trans(p):
    """Determine whether to make a transition or not."""
    return np.random.rand() <= p


def trans_prob(dE, t):
    """Probability of switching to a state with more energy."""
    return np.exp(-dE/t)


def make_trans(path):
    """Reverse a random part of the path."""
    path = path.copy()
    i = np.random.randint(len(path))
    j = np.random.randint(2, max(len(path) - i, 3))
    path[i : i + j] = np.flipud(path[i : i + j])
    return path


def annealing(cities, T, min_temp, decay=1, max_iter=1e5):
    """Simulate annealing"""
    # Random initilization
    rng = np.random.default_rng()
    state = rng.choice(len(cities), len(cities), replace=False)
    energies = []
    E = calc_energy(state, cities)
    i = 1
    while T > min_temp and i <= max_iter:
        energies.append(E)
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
        T *= i / (i + decay)
        i += 1
    return state, energies


def main():
    t_min = 1e-4
    init_t = 5
    axis_scale = 100
    n_cities = int(input('number of cities: '))
    np.random.seed(42)
    cities = np.random.rand(n_cities, 2)*axis_scale
    optimal_path, energies = annealing(cities, init_t, t_min, decay=0.4)

    plt.plot(np.arange(0, min(10000, len(energies))), energies[:min(10000, len(energies))])
    plt.show()
    print(calc_energy(optimal_path, cities))
    plt.scatter(cities[:, 0], cities[:, 1], c='r')
    plt.plot(cities[optimal_path][:, 0], cities[optimal_path][:, 1])
    plt.plot([cities[optimal_path[0]][0], cities[optimal_path[-1]][0]], 
        [cities[optimal_path[0]][1], cities[optimal_path[-1]][1]], c='b')
    plt.show()


if __name__ == "__main__":
    main()
