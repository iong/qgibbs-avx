if (APPLE)
	set(_libdir ENV DYLD_LIBRARY_PATH)
else()
	set(_libdir ENV LD_LIBRARY_PATH)
endif()

if(NOT CHOLMOD_LIBS)
	find_library(CHOLMOD_LIBS NAMES cholmod PATHS ${_libdir})
	if (NOT CHOLMOD_LIBS)
		message(FATAL_ERROR "CHOLMOD is needed!")
	endif()
endif()

get_filename_component(CHOLMOD_LIB ${CHOLMOD_LIBS} PATH)
get_filename_component(CHOLMOD_HOME ${CHOLMOD_LIB} PATH)

foreach(lib amd camd colamd ccolamd)
	find_library(${lib}_lib NAMES "${lib}" PATHS ${CHOLMOD_HOME}/lib NO_DEFAULT_PATH)
	set(CHOLMOD_LIBS ${CHOLMOD_LIBS} ${${lib}_lib})
endforeach()

find_library(metis_lib NAMES metis PATHS ${CHOLMOD_HOME}/lib NO_DEFAULT_PATH)
if (NOT metis_lib)
	find_library(metis_lib NAMES metis PATHS ${_libdir})
endif()
if (metis_lib)
	set(CHOLMOD_LIBS ${CHOLMOD_LIBS} ${metis_lib})
endif()

find_path(CHOLMOD_INC cholmod.h PATH_SUFFIXES ufsparse suitesparse HINTS ${CHOLMOD_HOME}/include)

message(STATUS "Found CHOLMOD in ${CHOLMOD_HOME}")
message(STATUS "\theaders in ${CHOLMOD_INC}")


