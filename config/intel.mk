CC:=icc
FC:=ifort
#Generic AMD
ifeq ($(CPU),amd)
    OPTFLAGS:=-ipo -O3 -msse3 -no-prec-div
#Generic Intel
else ifeq ($(CPU),intel)
    OPTFLAGS:=-ipo -O3 -no-prec-div -xSSE4.2
# Host
else
    OPTFLAGS:=-ipo -O3 -xHost
endif
OPTFLAGS += -heap-arrays
DBGFLAGS:=-O0 -g 
FDBG:=-fpe0 -traceback -check all -ftrapuv -warn unused
FFLAGS:=
LDFLAGS += -L/opt/hpc//intel-2011/lib
CPPFLAGS += -I/opt/hpc//intel-2011/include

LAPACK := -mkl=sequential
LIBS += -lcholmod -lamd -lcamd -lcolamd -lccolamd
