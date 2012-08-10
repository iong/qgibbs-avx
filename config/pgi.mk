CC=pgcc
FC=pgf95

COPT:=-fast
FOPT:=$(COPT)

CDBG:=-O0 -g -Ktrap=fp
FDBG:=-O0 -g -Ktrap=fp -Mbounds  

FFLAGS+=-mp=bind
LDFLAGS+=-mp=bind
