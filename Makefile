THIS_MAKEFILE:=$(realpath $(firstword $(MAKEFILE_LIST)))
VPATH:=$(dir $(THIS_MAKEFILE))
ifeq "$(VPATH)" ""
    VPATH:=.
endif


#DBG=1
COMPILER:=intel

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

FFLAGS += $(OPTFLAGS)
CFLAGS += $(OPTFLAGS)
LDFLAGS += $(OPTFLAGS)

LIBS += -lcholmod -lamd -lcamd -lcolamd -lccolamd $(LAPACK)

pimc:=utils.f90 pairint.f90 pimc.f90
dimer:=utils.f90 mvgwmodeffpot.f90 vgw.f90 vgwfm.f90 dlsode.f vgwspb_H2_4G_Rc_Q_tau_SqrtMeff_Mar03.f dimer.f90
cluster:=utils.f90 mvgwmodeffpot.f90 vgwspfm.f90  vgwfm.f90 dlsode.f vgwspb_H2_4G_Rc_Q_tau_SqrtMeff_Mar03.f cluster.f90
mccluster:=det_sparse_g.c utils.f90 mvgwmodeffpot.f90 vgw.f90  vgwspfm.f90 dlsode.f vgwspb_H2_4G_Rc_Q_tau_SqrtMeff_Mar03.f mccluster.f90
gibbs3:=utils.f90 pairint.f90 gibbs3.f90
gibbs4:=gibbs4.f90
gibbs3h:=gibbs3h.f90
qgibbs:=cholmod_logdet.c utils.f90 dlsode.f vgwspfm.f90 qgibbs.f90
objects=$(addsuffix .o,$(basename $(1)))


#all: $(addsuffix $(EXT),gibbs4 gibbs3h gibbs3)
all: cluster

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

debug:
	@echo $(PWD)
	$(MAKE) -f $(THIS_MAKEFILE) DBG=1

gibbs4$(EXT): $(call objects,$(gibbs4))
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

gibbs3h$(EXT): $(call objects,$(gibbs3h))
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

gibbs3$(EXT): $(call objects,$(gibbs3))
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

qgibbs: $(call objects,$(qgibbs))
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

dimer: $(call objects,$(dimer))
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

cluster: $(call objects,$(cluster))
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS) 

mccluster: $(call objects,$(mccluster))
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS) 

pimc: $(call objects,$(pimc))
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS) 

gibbs4.ps: $(gibbs4)
	a2ps -1 --borders=no -f 10 -o $@ $^

vgw.o: vgw0.f90 rhss0.f90 interaction_lists.f90 potential_energy.f90

vgwfm.o: vgw0fm.f90 rhssfm.f90
vgwspfm.o: vgw0spfm.f90 rhssspfm.f90 

clean:
	$(RM) *.o *.mod *.mod.F90 *.opari.inc opari.rc opari.tab.c *.mod.F

.PHONY: deps debug clean
