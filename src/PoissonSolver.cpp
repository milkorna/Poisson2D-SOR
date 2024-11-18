#include "PoissonSolver.h"
#include <cmath>
#include <stdexcept>

void PoissonSolver::initialize()
{
    size_t numRows = grid_->numRows();
    size_t numCols = grid_->numCols();

    // Initialize grid values to zero for interior points
    if (useMultithreading_) {
#pragma omp parallel for collapse(2)
        for (size_t i = 1; i < numRows - 1; ++i) {
            for (size_t j = 1; j < numCols - 1; ++j) {
                (*grid_)(i, j) = 0.0;
            }
        }
    } else {
        for (size_t i = 1; i < numRows - 1; ++i) {
            for (size_t j = 1; j < numCols - 1; ++j) {
                (*grid_)(i, j) = 0.0;
            }
        }
    }
}

int PoissonSolver::solve()
{
    const double tolerance = 1e-6;
    bool converged = false;
    int maxIterations = 10000;
    int iteration = 0;

    auto gridRect = std::dynamic_pointer_cast<RectangularGrid>(grid_);
    if (!gridRect) {
        throw std::runtime_error("Grid is not of type Grid");
    }

    auto solution_old = std::make_shared<RectangularGrid>(*gridRect);
    auto solution_new = std::make_shared<RectangularGrid>(*gridRect);

    while (!converged && iteration < maxIterations) {
        double err_sum = 0.0;

        if (useMultithreading_) {
            // Update red nodes in parallel
#pragma omp parallel for schedule(dynamic) reduction(+ : err_sum)
            for (size_t i = 1; i < solution_old->numRows() - 1; ++i) {
                for (size_t j = 1; j < solution_old->numCols() - 1; ++j) {
                    if ((i + j) % 2 == 0) { // Red nodes
                        double old_value = (*solution_old)(i, j);
                        double new_value = calculateSORUpdate(*solution_old, i, j);
                        (*solution_new)(i, j) = new_value;
                        err_sum += pow(new_value - old_value, 2);
                    }
                }
            }

            // Update black nodes in parallel
#pragma omp parallel for schedule(dynamic) reduction(+ : err_sum)
            for (size_t i = 1; i < solution_old->numRows() - 1; ++i) {
                for (size_t j = 1; j < solution_old->numCols() - 1; ++j) {
                    if ((i + j) % 2 != 0) { // Black nodes
                        double old_value = (*solution_old)(i, j);
                        double new_value = calculateSORUpdate(*solution_new, i, j);
                        (*solution_new)(i, j) = new_value;
                        err_sum += pow(new_value - old_value, 2);
                    }
                }
            }
        } else {
            // Sequential update of red nodes
            for (size_t i = 1; i < solution_old->numRows() - 1; ++i) {
                for (size_t j = 1; j < solution_old->numCols() - 1; ++j) {
                    if ((i + j) % 2 == 0) { // Red nodes
                        double old_value = (*solution_old)(i, j);
                        double new_value = calculateSORUpdate(*solution_old, i, j);
                        (*solution_new)(i, j) = new_value;
                        err_sum += pow(new_value - old_value, 2);
                    }
                }
            }

            // Sequential update of black nodes
            for (size_t i = 1; i < solution_old->numRows() - 1; ++i) {
                for (size_t j = 1; j < solution_old->numCols() - 1; ++j) {
                    if ((i + j) % 2 != 0) { // Black nodes
                        double old_value = (*solution_old)(i, j);
                        double new_value = calculateSORUpdate(*solution_new, i, j);
                        (*solution_new)(i, j) = new_value;
                        err_sum += pow(new_value - old_value, 2);
                    }
                }
            }
        }

        // Update boundary conditions in solution_new
        if (useMultithreading_) {
#pragma omp parallel for
            for (size_t i = 0; i < solution_new->numRows(); ++i) {
                (*solution_new)(i, 0) = (*gridRect)(i, 0);
                (*solution_new)(i, gridRect->numCols() - 1) = (*gridRect)(i, gridRect->numCols() - 1);
            }

#pragma omp parallel for
            for (size_t j = 0; j < solution_new->numCols(); ++j) {
                (*solution_new)(0, j) = (*gridRect)(0, j);
                (*solution_new)(gridRect->numRows() - 1, j) = (*gridRect)(gridRect->numRows() - 1, j);
            }
        } else {
            for (size_t i = 0; i < solution_new->numRows(); ++i) {
                (*solution_new)(i, 0) = (*gridRect)(i, 0);
                (*solution_new)(i, gridRect->numCols() - 1) = (*gridRect)(i, gridRect->numCols() - 1); // Ïðàâàÿ ãðàíèöà
            }

            for (size_t j = 0; j < solution_new->numCols(); ++j) {
                (*solution_new)(0, j) = (*gridRect)(0, j);
                (*solution_new)(gridRect->numRows() - 1, j) = (*gridRect)(gridRect->numRows() - 1, j);
            }
        }

        // Update the previous solution
        *solution_old = *solution_new;

        iteration++;
        double norm = sqrt(err_sum);
        converged = norm < tolerance;
    }

    // Copy the final solution to the main grid
    *gridRect = *solution_new;

    return iteration;
}

double PoissonSolver::calculateSORUpdate(const RectangularGrid& grid_values, size_t i, size_t j)
{
    double dx = grid_->xCoord(1) - grid_->xCoord(0);
    double dy = grid_->yCoord(1) - grid_->yCoord(0);

    double right_hand_side = testFunction_->valueAt(grid_->xCoord(j), grid_->yCoord(i));

    double coeff = omega_ / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy)));

    // Values from neighboring nodes
    double u_ip1j = grid_values(i + 1, j);
    double u_im1j = grid_values(i - 1, j);
    double u_ijp1 = grid_values(i, j + 1);
    double u_ijm1 = grid_values(i, j - 1);

    double new_value = (1.0 - omega_) * grid_values(i, j) +
                       coeff * ((u_ip1j + u_im1j) / (dx * dx) + (u_ijp1 + u_ijm1) / (dy * dy) - right_hand_side);

    return new_value;
}

int PoissonSolver::solveCircular()
{
    const double tolerance = 1e-6;
    bool converged = false;
    int maxIterations = 10000;
    int iteration = 0;

    auto circularGrid = std::dynamic_pointer_cast<CircularGrid>(grid_);
    if (!circularGrid) {
        throw std::runtime_error("Grid is not of type CircularGrid");
    }

    double dr = circularGrid->radialStepSize();
    double dtheta = circularGrid->angularStepSize();

    size_t numRadial = circularGrid->numRows();
    size_t numAngular = circularGrid->numCols();

    double dr2 = dr * dr;
    double dtheta2 = dtheta * dtheta;

    while (!converged && iteration < maxIterations) {
        double err_sum = 0.0;

        if (useMultithreading_) {
// Parallel version
#pragma omp parallel for reduction(+ : err_sum)
            for (size_t r = 1; r < numRadial - 1; ++r) {
                double currentR = circularGrid->rCoord(r);
                double r_inv = 1.0 / currentR;

                for (size_t theta = 0; theta < numAngular; ++theta) {
                    // Indices for neighboring points with periodic boundary conditions in theta
                    size_t rp1 = r + 1;
                    size_t rm1 = r - 1;
                    size_t tp1 = (theta + 1) % numAngular;
                    size_t tm1 = (theta + numAngular - 1) % numAngular;

                    // Values at neighboring points
                    double u_rp1 = (*circularGrid)(rp1, theta);
                    double u_rm1 = (*circularGrid)(rm1, theta);
                    double u_tp1 = (*circularGrid)(r, tp1);
                    double u_tm1 = (*circularGrid)(r, tm1);
                    double u_ij = (*circularGrid)(r, theta);

                    // Right-hand side function value at current point
                    double x = currentR * cos(circularGrid->thetaCoord(theta));
                    double y = currentR * sin(circularGrid->thetaCoord(theta));
                    double rhs = testFunction_->valueAt(x, y);

                    // Discretized radial derivatives
                    double u_rr = (u_rp1 - 2.0 * u_ij + u_rm1) / dr2;
                    double u_r = (u_rp1 - u_rm1) / (2.0 * dr);

                    // Discretized angular second derivative
                    double u_tt = (u_tp1 - 2.0 * u_ij + u_tm1) / dtheta2;

                    // Discretized Poisson equation
                    double denom = 2.0 / dr2 + 2.0 / (currentR * currentR * dtheta2);
                    double numer = (u_rp1 + u_rm1) / dr2 + (1.0 / currentR) * u_r +
                                   (u_tp1 + u_tm1) / (currentR * currentR * dtheta2) - rhs;
                    double u_new = numer / denom;

                    // Apply relaxation
                    double updated_value = (1.0 - omega_) * u_ij + omega_ * u_new;

                    // Accumulate error
                    double diff = updated_value - u_ij;
                    err_sum += diff * diff;

                    // Update grid value
                    (*circularGrid)(r, theta) = updated_value;
                }
            }
        } else {
            // Sequential version
            for (size_t r = 1; r < numRadial - 1; ++r) {
                double currentR = circularGrid->rCoord(r);
                double r_inv = 1.0 / currentR;

                for (size_t theta = 0; theta < numAngular; ++theta) {
                    // Indices for neighboring points with periodic boundary conditions in theta
                    size_t rp1 = r + 1;
                    size_t rm1 = r - 1;
                    size_t tp1 = (theta + 1) % numAngular;
                    size_t tm1 = (theta + numAngular - 1) % numAngular;

                    // Values at neighboring points
                    double u_rp1 = (*circularGrid)(rp1, theta);
                    double u_rm1 = (*circularGrid)(rm1, theta);
                    double u_tp1 = (*circularGrid)(r, tp1);
                    double u_tm1 = (*circularGrid)(r, tm1);
                    double u_ij = (*circularGrid)(r, theta);

                    // Right-hand side function value at current point
                    double x = currentR * cos(circularGrid->thetaCoord(theta));
                    double y = currentR * sin(circularGrid->thetaCoord(theta));
                    double rhs = testFunction_->valueAt(x, y);

                    // Discretized radial derivatives
                    double u_rr = (u_rp1 - 2.0 * u_ij + u_rm1) / dr2;
                    double u_r = (u_rp1 - u_rm1) / (2.0 * dr);

                    // Discretized angular second derivative
                    double u_tt = (u_tp1 - 2.0 * u_ij + u_tm1) / dtheta2;

                    // Discretized Poisson equation
                    double denom = 2.0 / dr2 + 2.0 / (currentR * currentR * dtheta2);
                    double numer = (u_rp1 + u_rm1) / dr2 + (1.0 / currentR) * u_r +
                                   (u_tp1 + u_tm1) / (currentR * currentR * dtheta2) - rhs;
                    double u_new = numer / denom;

                    // Apply relaxation
                    double updated_value = (1.0 - omega_) * u_ij + omega_ * u_new;

                    // Accumulate error
                    double diff = updated_value - u_ij;
                    err_sum += diff * diff;

                    // Update grid value
                    (*circularGrid)(r, theta) = updated_value;
                }
            }
        }

        iteration++;
        double norm = sqrt(err_sum);
        converged = norm < tolerance;
    }

    return iteration;
}
