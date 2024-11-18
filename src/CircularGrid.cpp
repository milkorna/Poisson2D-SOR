#include "CircularGrid.h"
#include <cmath>

CircularGrid::CircularGrid(size_t numRadial, size_t numAngular, double radius)
    : numRadial_(numRadial), numAngular_(numAngular), radius_(radius)
{
    dr_ = radius_ / (numRadial_ - 1);
    dtheta_ = 2.0 * M_PI / numAngular_;
    data_.resize(numRadial_, std::vector<double>(numAngular_, 0.0));
}

double& CircularGrid::operator()(size_t r, size_t theta)
{
    return data_[r][theta];
}

const double& CircularGrid::operator()(size_t r, size_t theta) const
{
    return data_[r][theta];
}

double CircularGrid::rCoord(size_t r) const
{
    return r * dr_;
}

double CircularGrid::thetaCoord(size_t theta) const
{
    return theta * dtheta_;
}

double CircularGrid::radius() const
{
    return radius_;
}

void CircularGrid::setBoundaryConditions(const TestFunction& boundaryFunc)
{
    for (size_t theta = 0; theta < numAngular_; ++theta) {
        double x = radius_ * cos(thetaCoord(theta));
        double y = radius_ * sin(thetaCoord(theta));
        (*this)(numRadial_ - 1, theta) = boundaryFunc.valueAt(x, y);
    }
}

double CircularGrid::radialStepSize() const
{
    return dr_;
}

double CircularGrid::angularStepSize() const
{
    return dtheta_;
}

double CircularGrid::xCoord(size_t j) const
{
    return radius_ * std::cos(static_cast<double>(j) * angularStepSize());
}

double CircularGrid::yCoord(size_t i) const
{
    return radius_ * std::sin(static_cast<double>(i) * angularStepSize());
}

size_t CircularGrid::numRows() const
{
    return numRadial_;
}

size_t CircularGrid::numCols() const
{
    return numAngular_;
}