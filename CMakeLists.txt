cmake_minimum_required(VERSION 3.12.0...4.1.2)

project(optimization LANGUAGES CXX)

set(${PROJECT_NAME}_VERSION 1.0)

# Set your C++ standard
# =====================
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Default build type if not specified
# ===================================
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the build type." FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")


# Set build flags
# ===============
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

option(OPTIMIZATION_USE_MPI "Build OPTIMIZATION with MPI support" ON)
# option(USE_MPI "Build DISCHARGES with MPI support" ON)

# === FIND MPI ===
if(OPTIMIZATION_USE_MPI)
    find_package(MPI REQUIRED)
    add_definitions(-DOPTIMIZATION_USE_MPI)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    message(STATUS "MPI found: ${MPI_CXX_COMPILER}")
endif()


# MFEM, HYPRE and Metis install directories
# =========================================
if(NOT OPTIMIZATION_USE_MPI)
    set(MFEM_DIR "$ENV{MFEM_DIR_SERIAL}")
else()
    set(MFEM_DIR "$ENV{MFEM_DIR_PARALLEL}")
    set(HYPRE_DIR "$ENV{HYPRE_DIR}")
    set(METIS_DIR "$ENV{METIS_DIR}")
endif()
if(NOT MFEM_DIR)
    message(FATAL_ERROR "Please set MFEM_DIR to the MFEM install directory")
endif()
# if(USE_MPI AND (NOT HYPRE_DIR))
if(OPTIMIZATION_USE_MPI AND (NOT HYPRE_DIR))
    message(FATAL_ERROR "Please set HYPRE_DIR to the HYPRE install directory")
endif()
# if(USE_MPI AND (NOT METIS_DIR))
if(OPTIMIZATION_USE_MPI AND (NOT METIS_DIR))
    message(FATAL_ERROR "Please set METIS_DIR to the METIS install directory")
endif()

# Include MFEM headers
# ====================
include_directories(${MFEM_DIR}/mfem/include)
include_directories(${MFEM_DIR}/mfem/include/mfem)
# if(USE_MPI)
if(OPTIMIZATION_USE_MPI)
    include_directories(${HYPRE_DIR}/include)
    include_directories(${METIS_DIR})
    include_directories(${METIS_DIR}/Lib)
endif()

# message(STATUS "THE FOLDER MFEM found : ${MFEM_DIR}")
# message(FATAL_ERROR "MY_REQUIRED_VAR is not set. Exiting.")

# Link MFEM, HYPRE and Metis static libraries
# ===========================================
add_library(mfem STATIC IMPORTED)
set_target_properties(mfem PROPERTIES 
    IMPORTED_LOCATION ${MFEM_DIR}/libmfem.a 
    INTERFACE_INCLUDE_DIRECTORIES ${MFEM_DIR}/mfem/include)

# if(USE_MPI)
if(OPTIMIZATION_USE_MPI)
    add_library(HYPRE STATIC IMPORTED)
    set_target_properties(HYPRE PROPERTIES 
        IMPORTED_LOCATION ${HYPRE_DIR}/lib/libHYPRE.a 
        INTERFACE_INCLUDE_DIRECTORIES ${HYPRE_DIR}/include)

    add_library(metis STATIC IMPORTED)
    set_target_properties(metis PROPERTIES 
        IMPORTED_LOCATION ${METIS_DIR}/libmetis.a 
        INTERFACE_INCLUDE_DIRECTORIES ${METIS_DIR}/Lib)
endif()

# Set the directory for executables 
# =================================
message(STATUS "Binary dir is: ${CMAKE_BINARY_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Add subdirectories
# ==================
add_subdirectory(src)
# <<<<<<< HEAD
add_subdirectory(examples)
# =======
# <<<<<<< HEAD
# add_subdirectory(sandbox/discharges)

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# add executables
# ===============
#add_executable( main.cpp)
#add_executable(ppm main.cpp)

# Link against MFEM
# =================
#target_link_libraries(MyApp PRIVATE mfem)
#target_link_libraries(ppm PRIVATE mfem)

# Optional: If MFEM was compiled with MPI, LAPACK, METIS, etc., link those too
# (You might need to manually link these depending on how MFEM was built)
# Example:
# target_link_libraries(MyApp PRIVATE mfem -lmetis -llapack -lblas -lHYPRE -lmpi)

# Optional: Add additional include or lib directories
#set(PPM_DIR "/home/innocent/softwares/src/ppm") 
#include_directories(/home/innocent/softwares/src/ppm/build/include)
#link_directories(/home/innocent/softwares/src/ppm/build/lib)
# =======
# add_subdirectory(examples)
# >>>>>>> 78a09d1 (The optimization using SIMP and a constant step is working!)
# >>>>>>> f85cf13 (One more commit after rebase)
