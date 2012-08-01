ifndef BLAS
	ifeq $(COMPILER) "intel"
		BLAS:=mkl_sequential
	else ifeq $(COMPILER) "GNU"
		BLAS:=blas
		ifeq $(OS) "darwin"
			BLAS:="Apple"
		endif
	else ifeq $(COMPILER) "PGI"
		BLAS:=acml
	endif
endif

ifeq $(BLAS) "mkl_sequential"
	ifeq $(COMPILER) "intel"
		BLAS_LIBRARIES:=-mkl=sequential
	else ifeq $(COMPILER) "GNU"
		BLAS_LIBRARIES=-lmkl_gf_lp64 -lmkl_sequential -lmkl_core
	else ifeq $(COMPILER) "PGI"
		BLAS_LIBRARIES=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core
	endif
else ifeq $(BLAS) "mkl_parallel"
	ifeq $(COMPILER) "intel"
		BLAS_LIBRARIES:=-mkl=parallel
	else ifeq $(COMPILER) "GNU"
		BLAS_LIBRARIES=-lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp
	else ifeq $(COMPILER) "PGI"
		BLAS_LIBRARIES=-lmkl_intel_lp64 -lmkl_pgi_thread -lmkl_core
	endif
else ifeq $(BLAS) "Apple"
	BLAS_LIBRARIES:=-framework Accelerate
else
	BLAS_LIBRARIES:=-l$(BLAS)
endif
