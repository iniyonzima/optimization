# Create examples 
# ===============
if(NOT OPTIMIZATION_USE_MPI)
	add_executable(density_optimization density_optimization.cpp)
	target_link_libraries(density_optimization PUBLIC optimization_fem mfem)
else()
	add_executable(density_optimization_p density_optimization_p.cpp)
	target_link_libraries(density_optimization_p PUBLIC optimization_fem mfem HYPRE metis MPI::MPI_CXX)
endif()
