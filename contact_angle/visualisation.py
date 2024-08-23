import matplotlib.pyplot as plt
import numpy as np
from .gibbs_dividing_plane import tanh_function, fit_tanh, find_tanh_midpoint


def plot_droplet_analysis(
    all_r, all_z, valid_r, valid_z, xc, yc, R, baseline, max_radius, plot_every=1000
):
    """Plot the droplet analysis results"""
    fig, ax = plt.subplots(figsize=(8, 8))

    x_min, x_max = -max_radius, max_radius
    ax.set_xlim(x_min, x_max)

    y_range = x_max - x_min
    y_mid = (np.min(all_z) + np.max(all_z)) / 2
    y_min = y_mid - y_range / 2
    y_max = y_mid + y_range / 2
    ax.set_ylim(y_min, y_max)

    ax.set_aspect("equal", adjustable="box")
    ax.scatter(
        all_r[::plot_every],
        all_z[::plot_every],
        alpha=0.1,
        s=0.04,
        color="blue",
        rasterized=True,
    )
    ax.scatter(valid_r, valid_z, color="red", s=10, label="Gibbs dividing plane")
    ax.scatter([-r for r in valid_r], valid_z, color="red", s=10)

    theta = np.linspace(0, np.pi, 100)
    x_circle = xc + R * np.cos(theta)
    y_circle = yc + R * np.sin(theta)
    ax.plot(x_circle, y_circle, "g-", label="Fitted circle")

    ax.set_xlabel("r (Å)")
    ax.set_ylabel("z (Å)")
    ax.legend()
    ax.grid(True)
    ax.axhline(y=baseline, color="k", linestyle="--", linewidth=0.5)
    ax.axvline(x=0, color="k", linestyle="--", linewidth=0.5)
    plt.tight_layout()

    return fig, ax


def plot_radial_density_profiles(
    density, r_centers, z_centers, num_profile_plots=1, custom_slices=None
):
    """
    Plot radial density profiles with tanh fits.

    Parameters:
    - density: 2D array of density values
    - r_centers, z_centers: Coordinates for density plot
    - num_profile_plots: Number of radial density profile plots to display
    - custom_slices: List of specific slice indices to plot
    """
    if custom_slices is None:
        plot_indices = np.linspace(
            0, density.shape[1] - 1, min(num_profile_plots, density.shape[1]), dtype=int
        )
    else:
        if any(idx >= density.shape[1] for idx in custom_slices):
            raise ValueError(
                f"One or more custom slice indices exceed the number of available slices ({density.shape[1]})."
            )
        plot_indices = custom_slices
        num_profile_plots = len(plot_indices)

    rows = int(np.ceil(np.sqrt(num_profile_plots)))
    cols = int(np.ceil(num_profile_plots / rows))
    fig, axs = plt.subplots(rows, cols, figsize=(cols * 5, rows * 5))
    axs = axs.ravel() if num_profile_plots > 1 else [axs]

    for i, idx in enumerate(plot_indices):
        midpoint = find_tanh_midpoint(r_centers, density[:, idx])

        ax = axs[i]
        ax.plot(r_centers, density[:, idx], "b.")

        if midpoint is not None:
            r_dense = np.linspace(r_centers[0], r_centers[-1], 1000)
            popt = fit_tanh(r_centers, density[:, idx])
            ax.plot(r_dense, tanh_function(r_dense, *popt), "r-", label="Tanh fit")
            ax.axvline(x=midpoint, color="g", linestyle="--", label="Midpoint")

        ax.set_xlabel("Radial distance (Å)")
        ax.set_ylabel("Number Density (a.u.)")
        ax.set_title(f"Slice {idx+1}, z = {z_centers[idx]:.2f} Å")
        ax.legend()

    # Remove any unused subplots
    for j in range(num_profile_plots, len(axs)):
        fig.delaxes(axs[j])

    fig.tight_layout()
    return fig, axs
