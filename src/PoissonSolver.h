#ifndef POISSON_SOLVER_H
#define POISSON_SOLVER_H

#include "CircularGrid.h"
#include "Functions.h"
#include "RectangularGrid.h"
#include <memory>
#include <omp.h>

/**
 * @brief Class that solves the Poisson equation using the Successive Over-Relaxation (SOR) method.
 */
class PoissonSolver {
public:
    /**
     * @brief Constructor for the PoissonSolver class.
     *
     * @param omega The relaxation factor used in the SOR algorithm (should be between 1 and 2).
     * @param grid The computational grid on which the equation is solved.
     * @param testFunction A pointer to the function used to initialize the right-hand side of the equation.
     * @param useMultithreading Flag to enable or disable multithreading.
     */
    PoissonSolver(double omega, const std::shared_ptr<GridBase>& grid, std::shared_ptr<TestFunction> testFunction,
                  bool useMultithreading)
        : omega_(omega), grid_(grid), testFunction_(testFunction), useMultithreading_(useMultithreading),
          maxIterations_(10000), tolerance_(1e-6)
    {
    }

    /**
     * @brief Initializes the grid by setting initial values and applying boundary conditions.
     *
     * Initializes the internal grid values to zero at all interior points.
     * Boundary conditions are applied separately through the grid's boundary condition functions.
     */
    void initialize();

    /**
     * @brief Solves the Poisson equation on a rectangular grid using the SOR method until convergence.
     *
     * This function implements the Successive Over-Relaxation (SOR) method to solve the Poisson equation
     * on a rectangular grid. It supports both sequential and parallel execution based on the useMultithreading flag.
     *
     * @return The number of iterations performed.
     */
    int solve();

    /**
     * @brief Solves the Poisson equation on a circular grid using the SOR method until convergence.
     *
     * This function is specialized to handle the SOR solution of the Poisson equation on a circular grid.
     * The circular grid uses polar coordinates, and periodic boundary conditions are applied in the angular direction.
     *
     * @return The number of iterations performed.
     */
    int solveCircular();

    /**
     * @brief Calculates the next value in the SOR update sequence for a rectangular grid.
     *
     * This function calculates the updated value for a grid point at position (i, j) using neighboring points.
     * It takes into account the relaxation factor omega_ and computes the value according to the SOR algorithm.
     *
     * @param i The row index of the grid element to update.
     * @param j The column index of the grid element to update.
     * @return The updated value for the grid element at position (i, j).
     */
    double calculateSORUpdate(const RectangularGrid& grid_values, size_t i, size_t j);

    /**
     * @brief Getter for the maximum number of iterations.
     *
     * @return The maximum number of iterations allowed for convergence.
     */
    int maxIterations() const
    {
        return maxIterations_;
    }

    /**
     * @brief Getter for the convergence tolerance.
     *
     * @return The convergence tolerance.
     */
    double tolerance() const
    {
        return tolerance_;
    }

    std::shared_ptr<GridBase> grid_; ///< The computational grid used for solving the equation.

private:
    double omega_;                               ///< The relaxation factor used in the SOR method.
    std::shared_ptr<TestFunction> testFunction_; ///< Function representing the right-hand side of the equation.
    bool useMultithreading_;                     ///< Flag to control the use of multithreading.
    int maxIterations_;                          ///< The maximum number of iterations allowed for convergence.
    double tolerance_; ///< The convergence tolerance to determine when the solution has converged.
};

#endif // POISSON_SOLVER_H