import numpy as np
import MDAnalysis as mda


def cartesian_to_cylindrical(positions):
    """Convert cartesian to cylindrical coordinates"""
    xy = positions[:, :2]
    r = np.linalg.norm(xy, axis=1)
    z = positions[:, 2]
    return r, z


def compute_normalized_density(r, z, cell_size, max_radius, min_molecules):
    """Compute normalized density from cylindrical coordinates"""
    z_min, z_max = np.min(z), np.max(z)
    r_bins = np.arange(0, max_radius + cell_size, cell_size)
    z_bins = np.arange(z_min, z_max + cell_size, cell_size)

    hist, _, _ = np.histogram2d(r, z, bins=(r_bins, z_bins))

    r_centers = (r_bins[1:] + r_bins[:-1]) / 2
    z_centers = (z_bins[1:] + z_bins[:-1]) / 2
    shell_volumes = np.pi * (r_bins[1:] ** 2 - r_bins[:-1] ** 2) * cell_size

    mask = np.sum(hist, axis=0) >= min_molecules
    hist = hist[:, mask]
    z_centers = z_centers[mask]

    density = hist / shell_volumes[:, np.newaxis]

    return density, r_centers, z_centers


def aggregate_data(
    universe, selection, start_frame, num_frames, cell_size, max_radius, min_molecules
):
    """Aggregate data from multiple frames and compute density"""
    all_r, all_z_cyl = [], []

    for frame in range(start_frame, start_frame + num_frames):
        atoms = universe.select_atoms(selection)
        universe.trajectory[frame]
        positions = atoms.positions
        com = atoms.center_of_mass()
        centered_positions = positions - com
        r, z_cyl = cartesian_to_cylindrical(centered_positions)
        all_r.extend(r)
        all_z_cyl.extend(z_cyl)

    all_r, all_z_cyl = map(np.array, (all_r, all_z_cyl))

    density, r_centers, z_centers = compute_normalized_density(
        all_r, all_z_cyl, cell_size, max_radius, min_molecules
    )

    return density, r_centers, z_centers, all_r, all_z_cyl
