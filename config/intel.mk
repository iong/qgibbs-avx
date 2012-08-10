CC:=icc
FC:=ifort

#Generic AMD
ifeq ($(TARGET),amd)
    COPT:=-ipo -O3 -msse3 -no-prec-div
else ifeq ($(CPU),intel)
    COPT:=-ipo -O3 -no-prec-div -xSSE4.2
else
    COPT:=-ipo -O3 -xHost
endif
FOPT:=$(COPT)

FFLAGS += -heap-ararys

CDBG:=-O0 -g 
FDBG:=-fpe0 -traceback -check all -ftrapuv -warn unused
