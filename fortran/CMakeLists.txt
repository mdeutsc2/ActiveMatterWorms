cmake_minimum_required (VERSION 3.22.1)
#set(CMAKE_Fortran_COMPILER "/usr/bin/gfortran")
#set(CMAKE_C_COMPILER "/usr/bin/gcc")
#set(CMAKE_CXX_COMPILER "/usr/bin/g++")
project(
	amatter
	VERSION 1.0
	LANGUAGES Fortran)

 # find packages go here
 find_package(OpenMP)
if(OPENMP_FOUND AND NOT(CMAKE_BUILD_TYPE MATCHES Debug))
	#setting flags for fast run with openmp
	message("OpenMP Build")
	set(HAVE_OPENMP 1)
    # set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -Wall -cpp -optimize")# -fopenmp")
	set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -Wall -cpp -O3 -fopenmp")
endif()
if(CMAKE_BUILD_TYPE MATCHES Debug)
	#setting flags for debug run
	message("Debug Build")
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -cpp -Og")
endif()
message("Fortran Flags: ${CMAKE_Fortran_FLAGS}")
#add_library(MyLibExample simple_lib.cpp simplelib.hpp)

# link each target with other targets or add options

# adding something we can run
add_executable(amatter src/amatter.f90
					   src/active_module.f90
					   src/fluid_module.f90)

#link targets here
# if(OPENMP_FOUND)
# 	target_link_libraries(amatter PUBLIC OpenMP::OpenMP)
# endif()
#target_link_libraries(amatter PRIVATE MyLibExample)
