cmake_minimum_required(VERSION 2.8)

if(INTEL)
	message("Intel!")
endif()

PROJECT(qgibbs Fortran)
if (CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
	set(LAPACK_LINKER_FLAGS -mkl=sequential)
	set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xHost -ipo")
	set(CMAKE_C_FLAGS_RELEASE "-O3 -xHost -ipo")

	find_library(CHOLMOD_LIBS NAMES cholmod PATHS /opt/hpc/intel-2011/lib NO_DEFAULT_PATH)

elseif(MAKE_Fortran_COMPILER_ID STREQUAL GNU)
	find_package(LAPACK REQUIRED)
endif()

include(FindCHOLMOD.cmake)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

if (APPLE)
	set(LINK_FLAGS -Wl,-no_pie)
endif()

#Targets!
ADD_EXECUTABLE(qgibbs utils.f90 dlsode.f vgw.f90 qgibbs.f90)

include_directories(qgibbs ${CHOLMOD_INC})
target_link_libraries(qgibbs ${CHOLMOD_LIBS} ${LAPACK_LIBRARIES})

set_target_properties(qgibbs PROPERTIES LINK_FLAGS ${LAPACK_LINKER_FLAGS})
