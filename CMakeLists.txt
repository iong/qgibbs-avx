cmake_minimum_required(VERSION 2.8)

if (NOT ENV{FC})
	set(ENV{CC} icc)
	set(ENV{FC} ifort)
endif()

PROJECT(qgibbs Fortran C)
if (CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
	set(MKL -mkl=sequential)
	set(LAPACK_LIBRARIES ${MKL})
	set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xHost -ipo")
	set(CMAKE_Fortran_FLAGS_DEBUG "-fpe0 -traceback -check all -ftrapuv -warn unused -O0 -g")
	set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O1 -g")
	set(CMAKE_C_FLAGS_RELEASE "-O3 -xHost -ipo")

	find_library(CHOLMOD_LIBS NAMES cholmod PATHS /opt/hpc/intel-2011/lib NO_DEFAULT_PATH)

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
	find_package(LAPACK REQUIRED)
	if (APPLE)
		message(WARNING "\nMKL for MacOSX has no gfortran interface.\n"
			"qgibbs_fastvgw will not be build!")
	else()
		find_library(MKL_GF NAMES mkl_gf_lp64 PATH_SUFFIXES intel64 ia32 HINTS $ENV{MKLROOT}/lib)
		get_filename_component(MKL_LIBDIR ${MKL_GF} PATH)
		set(MKL -Bstatic -Wl,--start-group -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -Bdynamic)
	endif()
endif()

include(FindCHOLMOD.cmake)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

if (APPLE)
	set(LINK_FLAGS -Wl,-no_pie)
endif()

#Targets!

# Quantum Gibbs with VGW-SP
ADD_EXECUTABLE(qgibbs utils.f90 dlsode.f vgw.f90 qgibbs.f90)
set_target_properties(qgibbs PROPERTIES COMPILE_FLAGS -fpp)

target_link_libraries(qgibbs ${LAPACK_LIBRARIES})
if (LAPACK_LINKER_FLAGS)
	set_target_properties(qgibbs PROPERTIES LINK_FLAGS ${LAPACK_LINKER_FLAGS})
endif()

if(NOT (CMAKE_Fortran_COMPILER_ID STREQUAL GNU AND APPLE))
# Quantum Gibbs with Fast-VGW
	ADD_EXECUTABLE(qgibbs_fastvgw cholmod_logdet.c utils.f90 dlsode.f vgwspfm.f90 qgibbs.f90)
	set_target_properties(qgibbs_fastvgw PROPERTIES COMPILE_FLAGS "-fpp -DVGWSPFM")
	include_directories(qgibbs_fastvgw ${CHOLMOD_INC})
	if (MKL_LIBDIR)
		link_directories(qgibbs_fastvgw ${MKL_LIBDIR})
	endif()
	target_link_libraries(qgibbs_fastvgw ${CHOLMOD_LIBS} ${MKL})
endif()
