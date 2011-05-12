THIS_MAKEFILE:=$(realpath $(firstword $(MAKEFILE_LIST)))
VPATH:=$(dir $(THIS_MAKEFILE))
ifeq "$(VPATH)" ""
    VPATH:=.
endif

#DBG=1
COMPILER:=gcc

OS=$(shell uname -s)
include $(VPATH)/config/$(COMPILER).mk
EXT:=
ifdef CPU
    EXT:=.$(CPU)
endif

ifdef DBG
	OPTFLAGS := $(DBGFLAGS)
	FFLAGS   += $(FDBG)
	LDFLAGS  += $(FDBG)
endif

LIBS:= $(LIBS) $(LAPACK)


gibbs3:=utils.f90 pairint.f90 gibbs3.f90
gibbs4:=gibbs4.f90
gibbs3h:=gibbs3h.f90
qgibbs:=utils.f90 vgw.f90 qgibbs.f90
objects=$(addsuffix .o,$(basename $(1)))


#all: $(addsuffix $(EXT),gibbs4 gibbs3h gibbs3)
all: qgibbs

ifneq ($(wildcard $(VPATH)/deps.mk),)
include $(VPATH)/deps.mk
endif

deps:
	$(RM) deps.mk
	gfortran -MM -cpp $(FFLAGS) $(gibbs4) $(gibbs3h) > deps.mk
	#$(VPATH)/f90deps $(ALL_SRC) > deps.mk

%.o : %.f90
	$(FC) $(FFLAGS) $(OPTFLAGS) -c $<

%.o : %.f
	$(FC) $(FFLAGS) $(OPTFLAGS) -c $<

debug:
	@echo $(PWD)
	$(MAKE) -f $(THIS_MAKEFILE) DBG=1

gibbs4$(EXT): $(call objects,$(gibbs4))
	$(FC) $(LDFLAGS) $(OPTFLAGS) -o $@ $^ $(LIBS)

gibbs3h$(EXT): $(call objects,$(gibbs3h))
	$(FC) $(LDFLAGS) $(OPTFLAGS) -o $@ $^ $(LIBS)

gibbs3$(EXT): $(call objects,$(gibbs3))
	$(FC) $(LDFLAGS) $(OPTFLAGS) -o $@ $^ $(LIBS)

qgibbs: $(call objects,$(qgibbs))
	$(FC) $(LDFLAGS) $(OPTFLAGS) -o $@ $^ $(LIBS)

gibbs4.ps: $(gibbs4)
	a2ps -1 --borders=no -f 10 -o $@ $^

vgw.f90: vgw0.f90 rhss0.f90 interaction_lists.f90 potential_energy.f90

clean:
	$(RM) *.o *.mod *.mod.F90 *.opari.inc opari.rc opari.tab.c *.mod.F

.PHONY: deps debug clean
