import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Directory for storing results
results_dir = "results"

# Create main directory if it doesn't exist
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

def load_data(test_file):
    """Load data from the specified file."""
    try:
        data = np.loadtxt(test_file)
        return data
    except Exception as e:
        print(f"Error loading file {test_file}: {e}")
        return None

def save_plot(fig, file_path):
    """Save plot to the specified file path."""
    try:
        plt.savefig(file_path)
    except Exception as e:
        print(f"Error saving plot to {file_path}: {e}")
    finally:
        plt.close(fig)

def plot_3d_surface(X, Y, Z, title, xlabel, ylabel, zlabel):
    """Create a 3D surface plot."""
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')
    try:
        ax.plot_surface(X, Y, Z, cmap='viridis')
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_zlabel(zlabel)
        return fig
    except Exception as e:
        print(f"Error plotting surface: {e}")
        plt.close(fig)
        return None

def solve_poisson_case_1():
    """
    Function to process test associated with the analytical solution u(x, y) = sin^2(x + y)
    and right-hand side f(x, y) = 8*cos^2(x + y) - 4.
    """
    data = load_data("poisson_solution.txt")
    if data is None or data.shape[1] < 3:
        print("Insufficient data in file. Expected columns: x y z.")
        return

    x, y, z_numerical = data[:, 0], data[:, 1], data[:, 2]

    # Filter data to include ranges 0 <= x <= 2 and 0 <= y <= 3
    valid_indices = (x >= 0) & (x <= 2) & (y >= 0) & (y <= 3)
    df = pd.DataFrame({'x': x[valid_indices], 'y': y[valid_indices], 'z': z_numerical[valid_indices]})
    df = df.drop_duplicates(subset=['x', 'y'])

    # Create pivot table
    try:
        pivot_table = df.pivot(index='y', columns='x', values='z')
    except ValueError as e:
        print(f"Error creating pivot table: {e}")
        return

    # Get unique x and y values from the pivot table
    x_unique = pivot_table.columns.values
    y_unique = pivot_table.index.values
    X, Y = np.meshgrid(x_unique, y_unique)
    Z_numerical = pivot_table.values

    # Handle NaN values in Z_numerical
    if np.isnan(Z_numerical).any():
        print("Warning: Missing values (NaN) in numerical solution data.")
        Z_numerical = np.nan_to_num(Z_numerical)

    # Compute exact values using the analytical solution u(x, y) = sin^2(x + y)
    Z_exact = np.sin(X + Y)**2

    # Create subdirectory for this test
    test_name = "poisson_case_1"
    test_dir = os.path.join(results_dir, test_name)
    if not os.path.exists(test_dir):
        os.makedirs(test_dir)

    # Save exact values to file
    exact_data = np.column_stack((X.flatten(), Y.flatten(), Z_exact.flatten()))
    exact_values_file_path = os.path.join(test_dir, "exact_values.txt")
    np.savetxt(exact_values_file_path, exact_data, fmt="%.6f", header="x y z_exact", comments='')

    # Plot numerical solution and save the result
    fig = plot_3d_surface(X, Y, Z_numerical, 'Numerical Poisson Solution for Case 1', 'X Axis', 'Y Axis', 'Z Axis')
    if fig:
        save_plot(fig, os.path.join(test_dir, "numerical_solution_plot.png"))

    # Plot exact solution and save the result
    fig = plot_3d_surface(X, Y, Z_exact, 'Exact Solution for Case 1', 'X Axis', 'Y Axis', 'Z Axis')
    if fig:
        save_plot(fig, os.path.join(test_dir, "exact_solution_plot.png"))

    print(f"Results for test '{test_name}' have been saved in directory '{test_dir}'.")

def solve_poisson_case_2():
    """
    Function to process test associated with solving the Poisson equation
    in a circular domain with boundary condition u(r=R) = 0 and right-hand side f(x, y) = -x * y.
    """
    data = load_data("poisson_solution_circular.txt")
    if data is None or data.shape[1] < 3:
        print("Insufficient data in file. Expected columns: x y z.")
        return

    x, y, z_numerical = data[:, 0], data[:, 1], data[:, 2]

    # Define grid dimensions
    num_radial, num_angular = 100, 360
    expected_size = num_radial * num_angular

    # Verify data size matches expected grid size
    if x.size != expected_size:
        print(f"Data size mismatch. Expected {expected_size} points, got {x.size}.")
        return

    # Reshape data
    try:
        Z_numerical = z_numerical.reshape((num_radial, num_angular))
        X = x.reshape((num_radial, num_angular))
        Y = y.reshape((num_radial, num_angular))
    except ValueError as e:
        print(f"Error reshaping data: {e}")
        return

    # Compute r and theta from x and y
    R = np.sqrt(X**2 + Y**2)
    Theta = np.mod(np.arctan2(Y, X), 2 * np.pi)

    # Compute the exact solution
    Z_exact = (R**2 * (1 - R**2) * np.sin(2 * Theta)) / 24

    # Create subdirectory for this test
    test_name = "poisson_case_2"
    test_dir = os.path.join(results_dir, test_name)
    if not os.path.exists(test_dir):
        os.makedirs(test_dir)

    # Save exact values to file
    exact_values = np.column_stack((X.flatten(), Y.flatten(), Z_exact.flatten()))
    exact_values_file_path = os.path.join(test_dir, "exact_values.txt")
    np.savetxt(exact_values_file_path, exact_values, fmt="%.6f", header="x y z_exact", comments='')

    # Handle NaN values in Z_numerical
    if np.isnan(Z_numerical).any():
        print("Warning: Missing values (NaN) in numerical solution data.")
        Z_numerical = np.nan_to_num(Z_numerical, nan=0.0)

    # Plot numerical solution and save the result
    fig = plot_3d_surface(X, Y, Z_numerical, 'Numerical Poisson Solution for Case 2', 'X Axis', 'Y Axis', 'Z Axis')
    if fig:
        save_plot(fig, os.path.join(test_dir, "numerical_solution_plot.png"))

    # Plot exact solution and save the result
    fig = plot_3d_surface(X, Y, Z_exact, 'Exact Solution for Case 2', 'X Axis', 'Y Axis', 'Z Axis')
    if fig:
        save_plot(fig, os.path.join(test_dir, "exact_solution_plot.png"))

    # Compute and visualize error
    error = np.abs(Z_numerical - Z_exact)
    fig, ax = plt.subplots(figsize=(8, 6))
    try:
        contour = ax.contourf(X, Y, error, levels=50, cmap='hot')
        fig.colorbar(contour)
        ax.set_xlabel('X Axis')
        ax.set_ylabel('Y Axis')
        ax.set_title('Error between Numerical and Exact Solutions')
        save_plot(fig, os.path.join(test_dir, "error_contour_plot.png"))
    except Exception as e:
        print(f"Error creating error contour plot: {e}")
        plt.close(fig)

    print(f"Results for test '{test_name}' have been saved in directory '{test_dir}'.")

if __name__ == "__main__":
    solve_poisson_case_2()
