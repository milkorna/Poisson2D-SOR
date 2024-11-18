#ifndef GRID_BASE_H
#define GRID_BASE_H

#include <Functions.h>
#include <cstddef> // For size_t definition

/**
 * @brief Abstract base class for a numerical grid.
 *
 * This class provides a common interface for different types of numerical grids,
 * allowing operations like accessing grid values, querying dimensions, and setting boundary conditions.
 */
class GridBase {
public:
    virtual ~GridBase() = default;

    /**
     * @brief Access grid value at position (i, j).
     *
     * @param i Row index.
     * @param j Column index.
     * @return Reference to the grid value at (i, j).
     */
    virtual double& operator()(size_t i, size_t j) = 0;

    /**
     * @brief Access grid value at position (i, j) (const version).
     *
     * @param i Row index.
     * @param j Column index.
     * @return Const reference to the grid value at (i, j).
     */
    virtual const double& operator()(size_t i, size_t j) const = 0;

    /**
     * @brief Get the x-coordinate of the specified column index.
     *
     * @param j Column index.
     * @return x-coordinate value.
     */
    virtual double xCoord(size_t j) const = 0;

    /**
     * @brief Get the y-coordinate of the specified row index.
     *
     * @param i Row index.
     * @return y-coordinate value.
     */
    virtual double yCoord(size_t i) const = 0;

    /**
     * @brief Get the number of rows in the grid.
     *
     * @return Number of rows.
     */
    virtual size_t numRows() const = 0;

    /**
     * @brief Get the number of columns in the grid.
     *
     * @return Number of columns.
     */
    virtual size_t numCols() const = 0;

    /**
     * @brief Set the boundary conditions for the grid.
     *
     * @param boundaryFunc Function representing the boundary conditions.
     */
    virtual void setBoundaryConditions(const TestFunction& boundaryFunc) = 0;
};

#endif // GRID_BASE_H