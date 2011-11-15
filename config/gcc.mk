ifeq ($(OS),Darwin)
	CC:=gcc-mp-4.6
	FC:=gfortran-mp-4.6
	CPPFLAGS:=-I/opt/local/include/ufsparse
	LDFLAGS += -L$(MKLROOT)/lib -L/opt/local/lib
	LAPACK:=-lmkl_gf_lp64 -lmkl_sequential -lmkl_core
	LIBS += -lmetis
else
	CC:=gcc-4.6
	FC:=gfortran-4.6
	CPPFLAGS:=-I/usr/include/suitesparse
	LDFLAGS += -L$(MKLROOT)/lib/intel64
	LAPACK:=-Bstatic -Wl,--start-group -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -Bdynamic
endif


OPTFLAGS:=-O3
DBGFLAGS:=-O0 -ggdb -Wall

CFLAGS:=
FFLAGS:=-ffree-line-length-0 -ffixed-line-length-0 -fimplicit-none
FDBG:=-fbounds-check -ffpe-trap=invalid,zero,overflow,denormal \
	-finit-real=SNAN -finit-integer=-1


