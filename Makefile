THIS_MAKEFILE:=$(realpath $(firstword $(MAKEFILE_LIST)))
VPATH:=$(dir $(THIS_MAKEFILE))
ifeq "$(VPATH)" ""
    VPATH:=.
endif


#DBG=1
COMPILER:=intel
BLAS:=mkl_sequential
BT:=opt
TARGET:=host

OS=$(shell uname -s)
include $(VPATH)/config/$(COMPILER).mk
include $(VPATH)/config/blas.mk
#include $(VPATH)/config/cholmod.mk

ifeq $(BT), "debug"
CFLAGS   += $(CDBG) 
FFLAGS   += $(FDBG) 
LDFLAGS  += $(FDBG)
else
FFLAGS += $(OPTFLAGS)
CFLAGS += $(OPTFLAGS)
LDFLAGS += $(OPTFLAGS)
endif


LIBS += $(BLAS)
LDFLAGS += $(BLAS_LDFLAGS)

clustergs:=cholmod_logdet.c  dlsode.f  utils.f90 vgw.f90 vgwspfm.f90  vgwfm.f90 clustergs.f90
qgibbs:=cholmod_logdet.c utils.f90 dlsode.f vgw.f90 qgibbs.f90
objects=$(addsuffix .o,$(basename $(1)))


all: pimcsc

ifneq ($(wildcard $(VPATH)/deps.mk),)
include $(VPATH)/deps.mk
endif

deps:
	$(RM) deps.mk
	gfortran -MM -cpp $(FFLAGS) $(gibbs4) $(gibbs3h) > deps.mk
	#$(VPATH)/f90deps $(ALL_SRC) > deps.mk

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

%.o : %.f
	$(FC) $(FFLAGS) -c $<


qgibbs: $(call objects,$(qgibbs))
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

clustergs: $(call objects,$(clustergs))
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS) 

%.ps: %.f90
	a2ps -1 --borders=no -f 10 -o $@ $^

vgw.o: vgw0.f90 rhss0.f90 interaction_lists.f90 potential_energy.f90 species.f90

vgwfm.o: vgw0fm.f90 rhssfm.f90 species.f90
vgwspfm.o: vgw0spfm.f90 rhssspfm.f90 species.f90

clean:
	$(RM) *.o *.mod

.PHONY: deps debug clean
