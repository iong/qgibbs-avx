CC:=icc
FC:=ifort
#Generic AMD
ifeq ($(CPU),amd)
    OPTFLAGS=-ipo -O3 -msse3 -no-prec-div
#Generic Intel
else ifeq ($(CPU),intel)
    OPTFLAGS=-ipo -O3 -no-prec-div -xSSE4.2
# Host
else
    OPTFLAGS=-fast
endif
DBGFLAGS=-O0 -g -warn unused
FDBG:=-fpe0 -traceback -check all -ftrapuv
FFLAGS:=-openmp
LDFLAGS:=$(LDFLAGS) -openmp

