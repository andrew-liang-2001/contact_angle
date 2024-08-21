import MDAnalysis as mda
import numpy as np
import yaml
import matplotlib.pyplot as plt
from contact_angle import (
    aggregate_data,
    find_gibbs_dividing_plane,
    fit_circle,
    calculate_contact_angle,
    plot_droplet_analysis,
)


def load_input_params(file_path):
    """Load input parameters from a YAML file."""
    with open(file_path, "r") as file:
        return yaml.safe_load(file)


def main():
    # Load input parameters
    params = load_input_params("input_params.yaml")

    # Initialize the MDAnalysis universe
    universe = mda.Universe(params["trajectory_file"])

    print("Starting water droplet contact angle analysis...")
    print(f"Analyzing trajectory: {params['trajectory_file']}")
    print(
        f"Frames: {params['start_frame']} to {params['start_frame'] + params['num_frames'] - 1}"
    )

    # Step 1: Aggregate data
    print("\nStep 1: Aggregating data...")
    density, r_centers, z_centers, all_r, all_z = aggregate_data(
        universe,
        params["water_selection"],
        params["start_frame"],
        params["num_frames"],
        params["cell_size"],
        params["max_radius"],
        params["min_molecules"],
    )
    print("Data aggregation complete.")

    # Step 2: Find Gibbs dividing plane
    print("\nStep 2: Finding Gibbs dividing plane...")
    valid_r, valid_z = find_gibbs_dividing_plane(density, r_centers, z_centers)

    if valid_r is not None and valid_z is not None:
        print("Gibbs dividing plane found successfully.")

        # Step 3: Fit circle and calculate contact angle
        print("\nStep 3: Fitting circle and calculating contact angle...")
        xc, yc, R = fit_circle(valid_r, valid_z)
        baseline = np.min(all_z)
        h = R + yc - baseline
        contact_angle = calculate_contact_angle(R, h)

        # Step 4: Visualize results
        print("\nStep 4: Generating visualization...")
        fig, ax = plot_droplet_analysis(
            all_r, all_z, valid_r, valid_z, xc, yc, R, baseline, params["max_radius"]
        )
        plt.savefig(params["output_plot"])
        print(f"Plot saved as: {params['output_plot']}")

        # Step 5: Save results
        print("\nStep 5: Saving results...")
        with open(params["output_data"], "w") as f:
            f.write(f"Fitted circle: Center=(0, {yc:.2f}), Radius={R:.2f}\n")
            f.write(f"Height above baseline: {h:.2f}\n")
            f.write(f"Estimated contact angle: {contact_angle:.2f}°\n")
            f.write(f"\nAnalysis parameters:\n")
            f.write(f"Cell size: {params['cell_size']} Å\n")
            f.write(f"Radial range: 0 to {params['max_radius']} Å\n")
            f.write(f"Height range: {z_centers[0]:.2f} to {z_centers[-1]:.2f} Å\n")
            f.write(f"Minimum molecules per height slice: {params['min_molecules']}\n")
            f.write(f"Number of valid height slices: {len(z_centers)}\n")
        print(f"Results saved as: {params['output_data']}")

        print("\nAnalysis complete. Check the output files for results.")
    else:
        print(
            "Error: Not enough valid data points to fit a circle. Try adjusting the analysis parameters."
        )


if __name__ == "__main__":
    main()
