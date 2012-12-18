cmake_minimum_required(VERSION 2.8)

# Change the default behavior by setting the environment variables CC, FC and
# BLA_VENDOR
# CC: icc, gcc, pgcc, etc.
# FC: ifort, gfortran, pgf90 etc.
# BLA_VENDOR: Intel10_lp64, ACML_MP, ACML_FM4_MP, ACML, etc.
set(CMAKE_C_COMPILER_INIT "icc")
set(CMAKE_Fortran_COMPILER_INIT "ifort")
set(OPENMP OFF CACHE BOOL "Enables OpenMP.")
set(ARCH "" CACHE STRING "GCC argument for -march=")
set(BLA_VENDOR "mkl_sequential")

if (APPLE)
    string(REPLACE ":" ";" CMAKE_LIBRARY_PATH $ENV{DYLD_LIBRARY_PATH})
else()
    string(REPLACE ":" ";" CMAKE_LIBRARY_PATH $ENV{LD_LIBRARY_PATH})
endif()

PROJECT(qgibbs C CXX Fortran)
if(OPENMP)
	find_package(OpenMP)
endif()

if (CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
    set(_arch "-xHost")
    if (ARCH)
        set(_arch "-march=${ARCH}")
    endif()
    set(CMAKE_C_FLAGS ${_arch})
    set(CMAKE_Fortran_FLAGS "${_arch} -heap-arrays")

    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ipo -no-prec-div")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -ipo -no-prec-div")

    set(CMAKE_Fortran_FLAGS_DEBUG "-fpe0 -traceback -ftrapuv -check all  -warn unused -O0 -g")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -g")
    set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -g")

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
    set(_arch "-march=native")
    if (ARCH)
        set(_arch "-march=${ARCH}")
    endif()

    set(CMAKE_C_FLAGS ${_arch})
    set(CMAKE_Fortran_FLAGS "${_arch} -fimplicit-none")

    set(CMAKE_Fortran_FLAGS_DEBUG "-fbacktrace -fbounds-check -ffpe-trap=invalid,zero,overflow,denormal -fimplicit-none -O0 -g")
    set(CMAKE_Fortran_FLAGS_DEBUG "-fbacktrace -fbounds-check -fimplicit-none -O0 -g")

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL PGI OR CMAKE_C_COMPILER_ID STREQUAL PGI)
    set(_arch "")
    if (ARCH)
        set(_arch "-tp=${ARCH}")
    endif()
    set(CMAKE_Fortran_FLAGS_RELEASE "-fast -O3 -Mvect=sse,simd:256 -Mipa=fast -Minline -Mfprelaxed ${_arch}")
    set(CMAKE_C_FLAGS_RELEASE "-fast -O3 -Mvect=sse,simd:256 -Mipa=fast -Minline -Mfprelaxed ${_arch}")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -Ktrap=fp -traceback")
endif()


if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


set(BLA_STATIC TRUE)
include(cmake/FindBLAS2.cmake)

#set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

if (APPLE)
    set(LINK_FLAGS -Wl,-no_pie)
endif()

#Targets!

# Quantum Gibbs with VGW-SP
set(EXE_EXT "")
if (ARCH)
    set(EXE_EXT ".${ARCH}")
endif()
ADD_EXECUTABLE(qgibbs${EXE_EXT} utils.f90 dlsode.f vgw.f90 gaussian_average_avx.c index_sort.c qgibbs.f90)
set_target_properties(qgibbs${EXE_EXT} PROPERTIES LINK_FLAGS "-static")
#ADD_EXECUTABLE(cluster utils.f90 dlsode.f vgw.f90 gaussian_average_avx.c index_sort.c vgwref.f90 xyz.f90 cluster.f90)

#add_executable(qsort qsort_test.f90 utils.f90 index_sort.c)

target_link_libraries(qgibbs${EXE_EXT} ${BLAS_LIBRARIES})
#target_link_libraries(cluster ${BLAS_LIBRARIES})
if (BLAS_LINKER_FLAGS)
    set_target_properties(qgibbs${EXE_EXT} PROPERTIES LINK_FLAGS ${BLAS_LINKER_FLAGS})
endif()
