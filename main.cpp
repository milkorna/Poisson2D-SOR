#include "Functions.h"
#include "PoissonSolver.h"
#include "RectangularGrid.h"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>

namespace fs = std::filesystem;

// Function to save the computed grid data to a file
// This function saves the solution for a rectangular grid to a text file
void saveData(const GridBase& grid, size_t numRows, size_t numCols, double xStart, double xEnd, double yStart,
              double yEnd)
{
    std::ofstream file("poisson_solution.txt");

    // Calculate the step sizes in x and y directions
    double dx = (xEnd - xStart) / (numCols - 1);
    double dy = (yEnd - yStart) / (numRows - 1);

    // Iterate through each grid point and save the x, y, and corresponding value
    for (size_t i = 0; i < numRows; ++i) {
        double y = yStart + i * dy;
        for (size_t j = 0; j < numCols; ++j) {
            double x = xStart + j * dx;
            file << x << " " << y << " " << grid(i, j) << "\n";
        }
    }

    file.close();
}

// Function to run the rectangular grid Poisson solver test
void runTest()
{
    // Create directory for storing results, if it does not already exist
    fs::create_directories("results/poisson_case_1");

    // Open log file for recording performance metrics
    std::ofstream logFile("results/poisson_case_1/performance_log.txt");
    if (!logFile.is_open()) {
        std::cerr << "Failed to open log file for writing" << std::endl;
        return;
    }

    // Define the grid parameters
    size_t numRows = 1000;           // number of rows
    size_t numCols = 1000;           // number of columns
    double xStart = 0.0, xEnd = 2.0; // x-axis range
    double yStart = 0.0, yEnd = 3.0; // y-axis range

    // Create a rectangular grid with the specified parameters
    auto grid = std::make_shared<RectangularGrid>(numRows, numCols, xStart, xEnd, yStart, yEnd);

    // Define boundary conditions for the grid
    Test1::BoundaryConditionLeft boundaryLeft;
    Test1::BoundaryConditionRight boundaryRight;
    Test1::BoundaryConditionBottom boundaryBottom;
    Test1::BoundaryConditionTop boundaryTop;

    // Set boundary conditions for each side of the grid
    grid->setBoundaryConditions(boundaryLeft, boundaryRight, boundaryBottom, boundaryTop);

    // Define the right-hand side function for the Poisson equation
    std::shared_ptr<TestFunction> testFunction = std::make_shared<Test1::RightHandSideFunction>();

    double omega = 1.5; // relaxation parameter for iterative solver

    // Solve the problem without multithreading
    {
        bool useMultithreading = false;
        PoissonSolver solver(omega, grid, testFunction, useMultithreading);

        // Record start time
        auto start = std::chrono::high_resolution_clock::now();
        solver.initialize();
        int iterations = solver.solve();
        auto end = std::chrono::high_resolution_clock::now(); // Record end time

        // Calculate and log the duration of the solver
        std::chrono::duration<double> duration = end - start;
        logFile << "Without multithreading: " << duration.count() << " seconds" << std::endl;
        if (iterations >= solver.maxIterations()) {
            logFile << "Reason for stopping: Reached maximum iterations (" << solver.maxIterations() << ")"
                    << std::endl;
        } else {
            logFile << "Reason for stopping: Convergence reached with tolerance " << solver.tolerance() << std::endl;
        }

        // Save solution to file (without multithreading)
        saveData(*solver.grid_, numRows, numCols, xStart, xEnd, yStart, yEnd);
    }

    // Solve the problem with multithreading
    {
        bool useMultithreading = true;
        PoissonSolver solver(omega, grid, testFunction, useMultithreading);

        // Record start time
        auto start = std::chrono::high_resolution_clock::now();
        solver.initialize();
        int iterations = solver.solve();
        auto end = std::chrono::high_resolution_clock::now(); // Record end time

        // Calculate and log the duration of the solver
        std::chrono::duration<double> duration = end - start;
        logFile << "With multithreading: " << duration.count() << " seconds" << std::endl;
        if (iterations >= solver.maxIterations()) {
            logFile << "Reason for stopping: Reached maximum iterations (" << solver.maxIterations() << ")"
                    << std::endl;
        } else {
            logFile << "Reason for stopping: Convergence reached with tolerance " << solver.tolerance() << std::endl;
        }

        // Save solution to file (with multithreading)
        saveData(*solver.grid_, numRows, numCols, xStart, xEnd, yStart, yEnd);
    }

    logFile.close();
}

// Function to save the computed data for a circular grid to a file
void saveDataCircular(const CircularGrid& grid, size_t numRadial, size_t numAngular)
{
    std::ofstream file("poisson_solution_circular.txt");

    // Iterate through each grid point in polar coordinates and save the Cartesian (x, y) and corresponding value
    for (size_t r = 0; r < numRadial; ++r) {
        double radius = grid.rCoord(r);
        for (size_t theta = 0; theta < numAngular; ++theta) {
            double angle = grid.thetaCoord(theta);
            file << radius * cos(angle) << " " << radius * sin(angle) << " " << grid(r, theta) << "\n";
        }
    }

    file.close();
}

// Function to run the circular grid Poisson solver test
void runCircularTest()
{
    // Create directory for storing results, if it does not already exist
    fs::create_directories("results/poisson_case_2");

    // Open log file for recording performance metrics
    std::ofstream logFile("results/poisson_case_2/performance_log.txt");
    if (!logFile.is_open()) {
        std::cerr << "Failed to open log file for writing" << std::endl;
        return;
    }

    // Define parameters for circular grid
    size_t numRadial = 100;  // number of radial divisions
    size_t numAngular = 360; // number of angular divisions
    double radius = 1.0;     // radius of the circle

    // Create the circular grid with the specified parameters
    auto circularGrid = std::make_shared<CircularGrid>(numRadial, numAngular, radius);

    // Set boundary conditions for the circular grid (e.g., u(r=R) = 0)
    Test2::BoundaryConditionCircular boundaryFunc;
    circularGrid->setBoundaryConditions(boundaryFunc);

    // Define the right-hand side function for the Poisson equation
    std::shared_ptr<TestFunction> circularTestFunction = std::make_shared<Test2::RightHandSideFunction>();
    double omega = 1.5; // relaxation parameter for iterative solver

    // Solve the problem without multithreading (Circular Grid)
    {
        bool useMultithreading = false;
        PoissonSolver solver(omega, circularGrid, circularTestFunction, useMultithreading);

        // Record start time
        auto start = std::chrono::high_resolution_clock::now();
        solver.initialize();
        int iterations = solver.solveCircular();
        auto end = std::chrono::high_resolution_clock::now(); // Record end time

        // Calculate and log the duration of the solver
        std::chrono::duration<double> duration = end - start;
        logFile << "Without multithreading (Circular Grid): " << duration.count() << " seconds" << std::endl;
        if (iterations >= solver.maxIterations()) {
            logFile << "Reason for stopping: Reached maximum iterations (" << solver.maxIterations() << ")"
                    << std::endl;
        } else {
            logFile << "Reason for stopping: Convergence reached with tolerance " << solver.tolerance() << std::endl;
        }

        // Save solution to file (without multithreading)
        saveDataCircular(*circularGrid, numRadial, numAngular);
    }

    // Solve the problem with multithreading (Circular Grid)
    {
        bool useMultithreading = true;
        PoissonSolver solver(omega, circularGrid, circularTestFunction, useMultithreading);

        // Record start time
        auto start = std::chrono::high_resolution_clock::now();
        solver.initialize();
        int iterations = solver.solveCircular();
        auto end = std::chrono::high_resolution_clock::now(); // Record end time

        // Calculate and log the duration of the solver
        std::chrono::duration<double> duration = end - start;
        logFile << "With multithreading (Circular Grid): " << duration.count() << " seconds" << std::endl;
        if (iterations >= solver.maxIterations()) {
            logFile << "Reason for stopping: Reached maximum iterations (" << solver.maxIterations() << ")"
                    << std::endl;
        } else {
            logFile << "Reason for stopping: Convergence reached with tolerance " << solver.tolerance() << std::endl;
        }

        // Save solution to file (with multithreading)
        saveDataCircular(*circularGrid, numRadial, numAngular);
    }

    logFile.close();
}

// Main function
int main()
{
    // Run the test for rectangular grid
    // runTest();

    // Run the test for circular grid
    runCircularTest();

    return 0;
}
