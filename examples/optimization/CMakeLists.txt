# Create examples 
# ===============
if(NOT OPTMIZATION_USE_MPI)
	add_executable(density_optimization density_optimization.cpp)
	target_link_libraries(density_optimization PUBLIC optimization_fem mfem)
endif()