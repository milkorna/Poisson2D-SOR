# Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.10)

# Project's name
project(PoissonSolverProject)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

# Set the installation directory relative to the build directory
set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install")

# Add the option to enable OpenMP if the code uses it
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Include the directory with header files
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# Compile the executable from the source files
add_executable(PoissonSolver
    main.cpp
    src/RectangularGrid.cpp
    src/CircularGrid.cpp
    src/PoissonSolver.cpp
)

# Link the executable with OpenMP if available
if(OpenMP_CXX_FOUND)
    target_link_libraries(PoissonSolver PRIVATE OpenMP::OpenMP_CXX)
endif()

# Install the executable in the 'bin' directory inside 'install'
install(TARGETS PoissonSolver DESTINATION bin)
