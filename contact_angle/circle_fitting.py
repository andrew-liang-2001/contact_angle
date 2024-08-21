import numpy as np
from scipy import optimize


def calc_R(yc, x, y):
    """Calculate the distance of each 2D point from the center (0, yc)"""
    return np.sqrt(x**2 + (y - yc) ** 2)


def f_2(c, x, y):
    """Calculate the algebraic distance between the data points and the mean circle centered at (0, yc)"""
    yc, R = c
    Ri = calc_R(yc, x, y)
    return Ri - R


def fit_circle(x, y):
    """Fit a circle to the given points with x_c constrained to 0"""
    y_m = np.mean(y)
    R_guess = np.mean(np.sqrt(x**2 + y**2))
    center_estimate = y_m, R_guess
    center, _ = optimize.leastsq(f_2, center_estimate, args=(x, y))
    yc, R = center
    return 0, yc, R  # Return x_c = 0, y_c, and R


def calculate_contact_angle(R, h):
    """Calculate the contact angle given the radius and height of the droplet"""
    return np.degrees(np.arccos((R - h) / R))
