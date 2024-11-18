#include "RectangularGrid.h"

RectangularGrid::RectangularGrid(size_t rows, size_t cols, double xStart, double xEnd, double yStart, double yEnd)
    : rows_(rows), cols_(cols), xStart_(xStart), xEnd_(xEnd), yStart_(yStart), yEnd_(yEnd)
{
    data_.resize(rows_, std::vector<double>(cols_, 0.0));
    dx_ = (xEnd_ - xStart_) / (cols_ - 1);
    dy_ = (yEnd_ - yStart_) / (rows_ - 1);
}

void RectangularGrid::setBoundaryConditions(const TestFunction& boundaryFunc)
{
    for (size_t i = 0; i < rows_; ++i) {
        data_[i][0] = boundaryFunc.valueAt(xStart_, yCoord(i));
        data_[i][cols_ - 1] = boundaryFunc.valueAt(xEnd_, yCoord(i));
    }

    for (size_t j = 0; j < cols_; ++j) {
        data_[0][j] = boundaryFunc.valueAt(xCoord(j), yStart_);
        data_[rows_ - 1][j] = boundaryFunc.valueAt(xCoord(j), yEnd_);
    }
}

void RectangularGrid::setBoundaryConditions(const TestFunction& leftFunc, const TestFunction& rightFunc,
                                            const TestFunction& bottomFunc, const TestFunction& topFunc)
{
    for (size_t i = 0; i < rows_; ++i) {
        double y = yStart_ + i * dy_;
        (*this)(i, 0) = leftFunc.valueAt(xStart_, y);
        (*this)(i, cols_ - 1) = rightFunc.valueAt(xEnd_, y);
    }

    for (size_t j = 0; j < cols_; ++j) {
        double x = xStart_ + j * dx_;
        (*this)(0, j) = bottomFunc.valueAt(x, yStart_);
        (*this)(rows_ - 1, j) = topFunc.valueAt(x, yEnd_);
    }
}

double& RectangularGrid::operator()(size_t i, size_t j)
{
    return data_[i][j];
}

const double& RectangularGrid::operator()(size_t i, size_t j) const
{
    return data_[i][j];
}

size_t RectangularGrid::numRows() const
{
    return rows_;
}

size_t RectangularGrid::numCols() const
{
    return cols_;
}

double RectangularGrid::xCoord(size_t j) const
{
    return xStart_ + j * dx_;
}

double RectangularGrid::yCoord(size_t i) const
{
    return yStart_ + i * dy_;
}
