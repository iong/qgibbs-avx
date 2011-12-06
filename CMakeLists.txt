cmake_minimum_required(VERSION 2.8)

if (NOT DEFINED ENV{FC})
    message(STATUS "NOT DEFINED ENV{FC}")
    set(ENV{CC} icc)
    set(ENV{FC} ifort)
endif()

PROJECT(qgibbs C Fortran)

if (CMAKE_Fortran_COMPILER_ID STREQUAL Intel)

    set(MKL_LIBRARIES -mkl=parallel)
    set(LAPACK_LIBRARIES -mkl=sequential)
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xHost -ipo -fpp")
    set(CMAKE_Fortran_FLAGS_DEBUG "-fpe0 -traceback -check all -ftrapuv -warn unused -O0 -g -fpp")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O1 -g -fpp")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -xHost -ipo")

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)

        set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fimplicit-none")
        set(CMAKE_Fortran_FLAGS_DEBUG "-fbounds-check -ffpe-trap=invalid,zero,overflow,denormal -fimplicit-none -O0 -g")
        set(CMAKE_C_FLAGS_RELEASE "-O3")

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL PGI OR CMAKE_C_COMPILER_ID STREQUAL PGI)

        set(CMAKE_Fortran_FLAGS_RELEASE "-fast")
        set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -Ktrap=fp -Mbounds")
        set(CMAKE_C_FLAGS_RELEASE "-fast")

endif()

include(FindMKL.cmake)

if (NOT DEFINED LAPACK_LIBRARIES)
	find_package(LAPACK REQUIRED)
endif()

include(FindCHOLMOD.cmake)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

if (APPLE)
    set(LINK_FLAGS -Wl,-no_pie)
endif()

#Targets!

# Quantum Gibbs with VGW-SP
ADD_EXECUTABLE(qgibbs utils.f90 dlsode.f vgw.f90 total_energy_vgw.f90 qgibbs.f90)

target_link_libraries(qgibbs ${LAPACK_LIBRARIES})
if (LAPACK_LINKER_FLAGS)
    set_target_properties(qgibbs PROPERTIES LINK_FLAGS ${LAPACK_LINKER_FLAGS})
endif()

if(CHOLMOD_LIBRARIES AND MKL_LIBRARIES)
# Quantum Gibbs with Fast-VGW
    ADD_EXECUTABLE(qgibbs_fastvgw cholmod_logdet.c utils.f90
            dlsode.f vgwspfm.f90 total_energy_vgwspfm.f90 qgibbs.f90)
    set_target_properties(qgibbs_fastvgw PROPERTIES COMPILE_FLAGS "-DVGWSPFM")
    include_directories(qgibbs_fastvgw ${CHOLMOD_INC})
    target_link_libraries(qgibbs_fastvgw ${CHOLMOD_LIBRARIES} ${MKL_LIBRARIES})
endif()
