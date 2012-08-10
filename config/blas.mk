blas_default_intel:=mkl_sequential
blas_default_GNU:=$(if $(findstring "darwin",$(OS)),"Apple","blas")
blas_default_PGI:=acml

blas_mkl_sequential_intel:=-mkl=sequential
blas_mkl_sequential_GNU:=-lmkl_gf_lp64 -lmkl_sequential -lmkl_core
blas_mkl_sequential_PGI:=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core

blas_mkl_parallel_intel:=-mkl=parallel
blas_mkl_parallel_GNU:=-lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp
blas_mkl_parallel_PGI:=-lmkl_intel_lp64 -lmkl_pgi_thread -lmkl_core

ifndef BLAS
	BLAS:=default
endif

ifdef blas_$(BLAS)_$(COMPILER)
	BLAS_LIBRARIES:=$(blas_$(BLAS)_$(COMPILER))
else ifeq $(BLAS) "Apple"
	BLAS_LIBRARIES:=-framework Accelerate
else
	BLAS_LIBRARIES:=-l$(BLAS)
endif


ifeq (mkl_intel,$(findstring mkl,$(BLAS))_$(COMPILER))
	BLAS_LIBRARIES:=-Wl,-Bstatic -Wl,--start-group $(BLAS_LIBRARIES) -Wl,--end-group -Wl,-Bdynamic
	ifeq ($(OS),Darwin)
		BLAS_LDFLAGS:= -L$(MKLROOT)/lib
	else
		BLAS_LDFLAGS:= -L$(MKLROOT)/lib/intel64
	endif
endif
