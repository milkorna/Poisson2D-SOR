#ifndef RECTANGULAR_GRID_H
#define RECTANGULAR_GRID_H

#include "GridBase.h"
#include <cstddef>
#include <vector>

/**
 * @brief Represents a two-dimensional numerical grid.
 *
 * This class provides a grid abstraction, which is commonly used for numerical
 * simulations in two dimensions, including operations for accessing and setting
 * grid values, as well as querying its dimensions and coordinates.
 */
class RectangularGrid : public GridBase {
public:
    /**
     * @brief Construct a RectangularGrid object.
     *
     * @param rows Number of rows in the grid.
     * @param cols Number of columns in the grid.
     * @param xStart Starting value of the x-axis.
     * @param xEnd Ending value of the x-axis.
     * @param yStart Starting value of the y-axis.
     * @param yEnd Ending value of the y-axis.
     */
    RectangularGrid(size_t rows, size_t cols, double xStart, double xEnd, double yStart, double yEnd);

    /**
     * @brief Set the boundary conditions for the grid using a single boundary function.
     *
     * @param boundaryFunc Function representing the boundary conditions.
     */
    void setBoundaryConditions(const TestFunction& boundaryFunc) override;

    /**
     * @brief Set the boundary conditions for all four sides of the grid.
     *
     * @param leftFunc Function for the left boundary.
     * @param rightFunc Function for the right boundary.
     * @param bottomFunc Function for the bottom boundary.
     * @param topFunc Function for the top boundary.
     */
    void setBoundaryConditions(const TestFunction& leftFunc, const TestFunction& rightFunc,
                               const TestFunction& bottomFunc, const TestFunction& topFunc);

    /**
     * @brief Access grid value at position (i, j).
     *
     * @param i Row index.
     * @param j Column index.
     * @return Reference to the grid value at (i, j).
     */
    double& operator()(size_t i, size_t j) override;

    /**
     * @brief Access grid value at position (i, j) (const version).
     *
     * @param i Row index.
     * @param j Column index.
     * @return Const reference to the grid value at (i, j).
     */
    const double& operator()(size_t i, size_t j) const override;

    /**
     * @brief Get the number of rows in the grid.
     *
     * @return Number of rows.
     */
    size_t numRows() const override;

    /**
     * @brief Get the number of columns in the grid.
     *
     * @return Number of columns.
     */
    size_t numCols() const override;

    /**
     * @brief Get the x-coordinate for the specified column index.
     *
     * @param j Column index.
     * @return x-coordinate value.
     */
    double xCoord(size_t j) const override;

    /**
     * @brief Get the y-coordinate for the specified row index.
     *
     * @param i Row index.
     * @return y-coordinate value.
     */
    double yCoord(size_t i) const override;

private:
    size_t rows_;                           ///< Number of rows in the grid.
    size_t cols_;                           ///< Number of columns in the grid.
    std::vector<std::vector<double>> data_; ///< 2D data storage for grid values.
    double xStart_;                         ///< Starting value of the x-axis.
    double xEnd_;                           ///< Ending value of the x-axis.
    double yStart_;                         ///< Starting value of the y-axis.
    double yEnd_;                           ///< Ending value of the y-axis.
    double dx_;                             ///< Step size in the x direction.
    double dy_;                             ///< Step size in the y direction.
};

#endif // RECTANGULAR_GRID_H