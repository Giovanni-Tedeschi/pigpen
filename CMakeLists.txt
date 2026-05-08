cmake_minimum_required(VERSION 3.15)
project(pigpen VERSION 1.0 LANGUAGES CXX)

# ---------------------------------------------------------------
# C++ standard
# ---------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------
find_package(HDF5 REQUIRED COMPONENTS CXX)
find_package(OpenMP REQUIRED)

# ---------------------------------------------------------------
# Executable
# ---------------------------------------------------------------
add_executable(pigpen
    run.cpp
    RiemannSolvers.cpp
    IO.cpp
    DragIntegrators.cpp
    BoundaryConditions.cpp
    Reconstruction.cpp
)

target_include_directories(pigpen PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(pigpen PRIVATE
    ${HDF5_CXX_LIBRARIES}
    OpenMP::OpenMP_CXX
)

# ---------------------------------------------------------------
# Compiler flags
# ---------------------------------------------------------------
target_compile_options(pigpen PRIVATE
    -O3                  # -O7 is silently capped at -O3 by gcc anyway
    -Wall                # useful warnings
)

# ---------------------------------------------------------------
# Optional: build type (Release/Debug)
# ---------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()