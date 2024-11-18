# Poisson2D-SOR: A Numerical Solver for the 2D Poisson Equation

## Overview

This project implements a numerical solver for the 2D Poisson equation using the Successive Over-Relaxation (SOR) method. The project supports solving the Poisson equation on both rectangular and circular grids and includes a comparison between numerical and exact solutions.

### Features:
- Supports solving the Poisson equation on both **rectangular** and **circular** domains.
- Implements **multithreading** to accelerate computation.
- Provides visualization of **numerical vs. exact solutions** for verification.
- Generates performance metrics including **convergence times**.

## Project Structure

The project consists of three main components:
- **Solver Implementation** (C++): Contains the solver logic for rectangular and circular domains.
- **Visualization Script** (Python): Processes the output data to compare numerical and exact solutions.
- **Data and Plots**: Saves performance logs, numerical solution plots, and comparison plots.

## Dependencies

- **C++17** or later for the solver.
- **Python 3** for visualization, with required libraries:
  - `numpy`
  - `pandas`
  - `matplotlib`
