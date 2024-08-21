import numpy as np
from scipy import optimize


def tanh_function(x, A, x0, d):
    return A * (1 - np.tanh(2 * (x - x0) / d))


def fit_tanh(r, density):
    """Fit a tanh function to the density profile"""
    A_guess = (np.max(density)) / 2
    x0_guess = np.mean(r)
    d_guess = (np.max(r) - np.min(r)) / 10

    p0 = [A_guess, x0_guess, d_guess]

    try:
        popt, _ = optimize.curve_fit(tanh_function, r, density, p0, method="lm")
        return popt
    except RuntimeError:
        print("Error: Tanh fitting failed.")
        return None


def find_tanh_midpoint(r_centers, density_slice):
    """Find the midpoint of the tanh function fit"""
    popt = fit_tanh(r_centers, density_slice)
    if popt is not None:
        _, x0, _ = popt
        return x0
    return None


def find_gibbs_dividing_plane(density, r_centers, z_centers):
    """Find the Gibbs equimolar dividing plane"""
    tanh_midpoints = [
        find_tanh_midpoint(r_centers, density[:, i]) for i in range(density.shape[1])
    ]

    # Filter out None values
    valid_data = [(r, z) for r, z in zip(tanh_midpoints, z_centers) if r is not None]

    if valid_data:
        valid_r, valid_z = zip(*valid_data)
        return np.array(valid_r), np.array(valid_z)
    else:
        return None, None
