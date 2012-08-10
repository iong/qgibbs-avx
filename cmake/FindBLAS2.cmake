include (FindLibraryList.cmake)

if($ENV{BLA_VENDOR} MATCHES ".+")
  set(BLA_VENDOR $ENV{BLA_VENDOR})
endif()

if(NOT BLA_VENDOR)
  set(BLA_VENDOR "All")
endif()

set(_blas_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})

if (BLA_STATIC)
	set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

if (BLA_VENDOR MATCHES "ACML.*" OR BLA_VENDOR STREQUAL "All")
  if(NOT BLAS_LIBRARIES)
      set(_lib "acml")
      if( BLA_VENDOR MATCHES "ACML.*MP" )
	set(_lib "acml_mp")
      endif()

      find_library(BLAS_LIBRARIES "${_lib}")
  endif()
endif()

# Apple BLAS library?
if (BLA_VENDOR STREQUAL "Apple" OR BLA_VENDOR STREQUAL "All")
  if(NOT BLAS_LIBRARIES)
    find_library(BLAS_LIBRARIES "Accelerate")
  endif()
endif()

if (BLA_VENDOR MATCHES "mkl_.*" OR BLA_VENDOR STREQUAL "All")
  if( NOT BLAS_LIBRARIES)

    if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
      set(_MKL_LP64 "mkl_intel_lp64")
      set(_MKL_THREAD "mkl_intel_thread;mkl_core;iomp5")
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      set(_MKL_LP64 "mkl_gf_lp64")
      set(_MKL_THREAD "mkl_gnu_thread;mkl_core;gomp")
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL PGI OR CMAKE_C_COMPILER_ID STREQUAL PGI)
      set(_MKL_LP64 "mkl_intel_lp64")
      set(_MKL_THREAD "mkl_pgi_thread;mkl_core")
    endif()

    if (BLA_VENDOR STREQUAL "mkl_sequential")
      set(_MKL_THREAD "mkl_sequential;mkl_core")
    endif()

    find_library_list(BLAS_LIBRARIES "${_MKL_LP64};${_MKL_THREAD}")
  endif()
endif()

if(BLAS_LIBRARIES)
	message("Found BLAS " ${BLAS_LIBRARIES})
endif()

set(CMAKE_FIND_LIBRARY_SUFFIXES ${_blas_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
