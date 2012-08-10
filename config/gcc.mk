CC:=gcc
FC:=gfortran

FFLAGS+=-ffree-line-length-0 -ffixed-line-length-0 -fimplicit-none

#Generic AMD
ifeq ($(TARGET),host)
    CFLAGS += -march=native
    FFLAGS += -march=native
else ifdef TARGET
    CFLAGS += -march=$(TARGET)
    FFLAGS += -march=$(TARGET)
endif

COPT:=-O3 -ffast-math
FOPT:=$(COPT)

CDBG:=-O0 -g -Wall
FDBG:=-fbounds-check -ffpe-trap=invalid,zero,overflow,denormal \
	-finit-real=SNAN -finit-integer=-1
