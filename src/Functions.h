#ifndef FUNCTIONS_H // Ensures that the header file is included only once in a single compilation.
#define FUNCTIONS_H

#include <cmath> // Include the cmath library for access to mathematical functions such as sin and cos.

/**
 * @brief Abstract base class for test functions.
 *
 * This class serves as an interface for defining functions that can return a value based
 * on x and y coordinates. It is intended for use in numerical methods that require
 * evaluation of functions at discrete grid points.
 */
class TestFunction {
public:
    /**
     * @brief Compute the function value at a given point (x, y).
     *
     * This is a pure virtual function that must be implemented by derived classes to
     * define specific function behaviors.
     *
     * @param x The x-coordinate of the point.
     * @param y The y-coordinate of the point.
     * @return double The function value at the given point.
     */
    virtual double valueAt(double x, double y) const = 0;

    /**
     * @brief Virtual destructor.
     *
     * Ensures derived classes can be cleaned up correctly through base class pointers.
     */
    virtual ~TestFunction()
    {
    }
};

/**
 * @namespace Test1
 *
 * @brief Contains test functions for use in rectangular grid simulations.
 *
 * This namespace includes classes that define the right-hand side of the Poisson equation
 * as well as boundary conditions for the four sides of a rectangular domain.
 */
namespace Test1 {

    /**
     * @brief Defines the right-hand side function of the Poisson equation.
     *
     * This class provides an implementation of a specific right-hand side function, which
     * is used in the Poisson solver to represent the source term.
     */
    class RightHandSideFunction : public TestFunction {
    public:
        /**
         * @brief Compute the value of the right-hand side function at a given point (x, y).
         *
         * @param x The x-coordinate of the point.
         * @param y The y-coordinate of the point.
         * @return double The function value at the given point.
         */
        double valueAt(double x, double y) const override
        {
            return 8.0 * std::cos(x + y) * std::cos(x + y) - 4.0;
        }
    };

    /**
     * @brief Defines the boundary condition on the left side of the rectangular domain.
     *
     * This class provides an implementation for the left boundary condition of the domain,
     * where the value of the function is defined as u(0, y) = sin^2(y).
     */
    class BoundaryConditionLeft : public TestFunction {
    public:
        /**
         * @brief Compute the boundary value at a given point (x, y).
         *
         * @param x The x-coordinate of the point.
         * @param y The y-coordinate of the point.
         * @return double The function value at the given point.
         */
        double valueAt(double x, double y) const override
        {
            return std::sin(y) * std::sin(y); // u(0, y) = sin^2(y)
        }
    };

    /**
     * @brief Defines the boundary condition on the right side of the rectangular domain.
     *
     * This class provides an implementation for the right boundary condition of the domain,
     * where the value of the function is defined as u(2, y) = sin^2(y + 2).
     */
    class BoundaryConditionRight : public TestFunction {
    public:
        /**
         * @brief Compute the boundary value at a given point (x, y).
         *
         * @param x The x-coordinate of the point.
         * @param y The y-coordinate of the point.
         * @return double The function value at the given point.
         */
        double valueAt(double x, double y) const override
        {
            return std::sin(y + 2) * std::sin(y + 2); // u(2, y) = sin^2(y + 2)
        }
    };

    /**
     * @brief Defines the boundary condition on the bottom side of the rectangular domain.
     *
     * This class provides an implementation for the bottom boundary condition of the domain,
     * where the value of the function is defined as u(x, 0) = sin^2(x).
     */
    class BoundaryConditionBottom : public TestFunction {
    public:
        /**
         * @brief Compute the boundary value at a given point (x, y).
         *
         * @param x The x-coordinate of the point.
         * @param y The y-coordinate of the point.
         * @return double The function value at the given point.
         */
        double valueAt(double x, double y) const override
        {
            return std::sin(x) * std::sin(x); // u(x, 0) = sin^2(x)
        }
    };

    /**
     * @brief Defines the boundary condition on the top side of the rectangular domain.
     *
     * This class provides an implementation for the top boundary condition of the domain,
     * where the value of the function is defined as u(x, 3) = sin^2(x + 3).
     */
    class BoundaryConditionTop : public TestFunction {
    public:
        /**
         * @brief Compute the boundary value at a given point (x, y).
         *
         * @param x The x-coordinate of the point.
         * @param y The y-coordinate of the point.
         * @return double The function value at the given point.
         */
        double valueAt(double x, double y) const override
        {
            return std::sin(x + 3) * std::sin(x + 3); // u(x, 3) = sin^2(x + 3)
        }
    };
}

/**
 * @namespace Test2
 *
 * @brief Contains test functions for use in circular grid simulations.
 *
 * This namespace includes classes that define the right-hand side of the Poisson equation
 * as well as the boundary condition for a circular domain.
 */
namespace Test2 {

    /**
     * @brief Defines the right-hand side function of the Poisson equation for a circular domain.
     *
     * This class provides an implementation of a specific right-hand side function, which
     * is used in the Poisson solver to represent the source term for a circular grid.
     */
    class RightHandSideFunction : public TestFunction {
    public:
        /**
         * @brief Compute the value of the right-hand side function at a given point (x, y).
         *
         * @param x The x-coordinate of the point.
         * @param y The y-coordinate of the point.
         * @return double The function value at the given point.
         */
        double valueAt(double x, double y) const override
        {
            return -x * y;
        }
    };

    /**
     * @brief Defines the boundary condition for the circular domain.
     *
     * This class provides an implementation for the boundary condition at the outer edge
     * of the circular domain, where the value of the function is set to zero.
     */
    class BoundaryConditionCircular : public TestFunction {
    public:
        /**
         * @brief Compute the boundary value at a given point (x, y).
         *
         * @param x The x-coordinate of the point.
         * @param y The y-coordinate of the point.
         * @return double The function value at the given point.
         */
        double valueAt(double x, double y) const override
        {
            return 0.0; // u(r=R) = 0
        }
    };
}

#endif // FUNCTIONS_H