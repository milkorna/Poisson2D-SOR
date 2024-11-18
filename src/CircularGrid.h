#ifndef CIRCULAR_GRID_H
#define CIRCULAR_GRID_H

#include "Functions.h"
#include "GridBase.h"
#include <vector>

/**
 * @brief Represents a circular numerical grid.
 *
 * This class provides a grid representation for problems defined in a circular domain.
 * It includes operations for accessing grid values, setting boundary conditions, and
 * querying the radial and angular coordinates.
 */
class CircularGrid : public GridBase {
public:
    /**
     * @brief Construct a CircularGrid object.
     *
     * @param numRadial Number of divisions in the radial direction.
     * @param numAngular Number of divisions in the angular direction.
     * @param radius Radius of the circular domain.
     */
    CircularGrid(size_t numRadial, size_t numAngular, double radius);

    /**
     * @brief Access grid value at position (r, theta).
     *
     * @param r Radial index.
     * @param theta Angular index.
     * @return Reference to the grid value at (r, theta).
     */
    double& operator()(size_t r, size_t theta) override;

    /**
     * @brief Access grid value at position (r, theta) (const version).
     *
     * @param r Radial index.
     * @param theta Angular index.
     * @return Const reference to the grid value at (r, theta).
     */
    const double& operator()(size_t r, size_t theta) const override;

    /**
     * @brief Get the radial coordinate for the specified radial index.
     *
     * @param r Radial index.
     * @return Radial coordinate value.
     */
    double rCoord(size_t r) const;

    /**
     * @brief Get the angular coordinate for the specified angular index.
     *
     * @param theta Angular index.
     * @return Angular coordinate value in radians.
     */
    double thetaCoord(size_t theta) const;

    /**
     * @brief Get the radius of the circular domain.
     *
     * @return Radius of the circular domain.
     */
    double radius() const;

    /**
     * @brief Set the boundary conditions for the grid.
     *
     * @param boundaryFunc Function representing the boundary conditions.
     */
    void setBoundaryConditions(const TestFunction& boundaryFunc) override;

    /**
     * @brief Get the step size in the radial direction.
     *
     * @return Radial step size.
     */
    double radialStepSize() const;

    /**
     * @brief Get the step size in the angular direction.
     *
     * @return Angular step size.
     */
    double angularStepSize() const;

    /**
     * @brief Get the x-coordinate of the specified angular index (not applicable for CircularGrid).
     *
     * @param j Angular index.
     * @return x-coordinate value.
     */
    double xCoord(size_t j) const override;

    /**
     * @brief Get the y-coordinate of the specified radial index (not applicable for CircularGrid).
     *
     * @param i Radial index.
     * @return y-coordinate value.
     */
    double yCoord(size_t i) const override;

    /**
     * @brief Get the number of radial divisions.
     *
     * @return Number of radial divisions.
     */
    size_t numRows() const override;

    /**
     * @brief Get the number of angular divisions.
     *
     * @return Number of angular divisions.
     */
    size_t numCols() const override;

private:
    size_t numRadial_;                      ///< Number of radial divisions.
    size_t numAngular_;                     ///< Number of angular divisions.
    double radius_;                         ///< Radius of the circular domain.
    double dr_;                             ///< Radial step size.
    double dtheta_;                         ///< Angular step size.
    std::vector<std::vector<double>> data_; ///< 2D data storage for grid values.
};

#endif // CIRCULAR_GRID_H