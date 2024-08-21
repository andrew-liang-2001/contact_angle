import matplotlib.pyplot as plt
import numpy as np


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
