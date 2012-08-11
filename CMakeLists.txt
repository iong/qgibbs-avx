cmake_minimum_required(VERSION 2.8)

# Change the default behavior by setting the environment variables CC, FC and
# BLA_VENDOR
# CC: icc, gcc, pgcc, etc.
# FC: ifort, gfortran, pgf90 etc.
# BLA_VENDOR: Intel10_lp64, ACML_MP, ACML_FM4_MP, ACML, etc.
set(CMAKE_C_COMPILER_INIT "icc")
set(CMAKE_Fortran_COMPILER_INIT "ifort")
set(BLA_VENDOR "mkl_sequential")

if (APPLE)
    set(CMAKE_LIBRARY_PATH string(REPLACE ":" " " $ENV{DYLD_LIBRARY_PATH}))
else()
    set(CMAKE_LIBRARY_PATH string(REPLACE ":" " " $ENV{LD_LIBRARY_PATH}))
endif()

PROJECT(qgibbs C Fortran)

if (CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
    set(CMAKE_C_FLAGS "-mavx")
    set(CMAKE_Fortran_FLAGS "-heap-arrays")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ipo -no-prec-div")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -ipo -no-prec-div")

    set(CMAKE_Fortran_FLAGS_DEBUG "-fpe0 -traceback -ftrapuv -check all  -warn unused -O0 -g")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -g")
    set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -g")

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
    set(CMAKE_C_FLAGS "-march=native")
    set(CMAKE_Fortran_FLAGS "-march=native -fimplicit-none")

    set(CMAKE_Fortran_FLAGS_DEBUG "-fbacktrace -fbounds-check -ffpe-trap=invalid,zero,overflow,denormal -fimplicit-none -O0 -g")
    set(CMAKE_Fortran_FLAGS_DEBUG "-fbacktrace -fbounds-check -fimplicit-none -O0 -g")

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL PGI OR CMAKE_C_COMPILER_ID STREQUAL PGI)
    set(CMAKE_Fortran_FLAGS_RELEASE "-fast -O3 -Mvect=sse,simd:256 -Mipa=fast -Minline -tp=bulldozer-64 -Mfprelaxed")
    set(CMAKE_C_FLAGS_RELEASE "-fast -O3 -Mvect=sse,simd:256 -Mipa=fast -Minline -tp=bulldozer-64 -Mfprelaxed")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -Ktrap=fp -traceback")
endif()

set(BLA_STATIC TRUE)
include(cmake/FindBLAS2.cmake)

#set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

if (APPLE)
    set(LINK_FLAGS -Wl,-no_pie)
endif()

#Targets!

# Quantum Gibbs with VGW-SP
ADD_EXECUTABLE(qgibbs utils.f90 dlsode.f vgw.f90 pinvdetinv.c gaussian_average_avx.c index_sort.c total_energy_vgw.f90 qgibbs.f90)
ADD_EXECUTABLE(cluster utils.f90 dlsode.f vgw.f90 pinvdetinv.c gaussian_average_avx.c index_sort.c vgwref.f90 xyz.f90 cluster.f90)

add_executable(qsort qsort_test.f90 utils.f90 index_sort.c)

target_link_libraries(qgibbs ${BLAS_LIBRARIES})
target_link_libraries(cluster ${BLAS_LIBRARIES})
if (BLAS_LINKER_FLAGS)
    set_target_properties(qgibbs PROPERTIES LINK_FLAGS ${BLAS_LINKER_FLAGS})
endif()
