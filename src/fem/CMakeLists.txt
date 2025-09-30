# Create the static library for the fem directory
# ===============================================
add_library(optimization_fem
    datacollection.cpp
)

# Make headers accessible to other targets
# ========================================
target_include_directories(optimization_fem PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
