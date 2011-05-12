CC=pgcc
FC=pgf95
OPTFLAGS=-fast
DBGFLAGS=-O0 -g -Ktrap=fp
FDBG:=-Mbounds  
FFLAGS:=-mp=bind
LDFLAGS:=$(LDFLAGS) $(FFLAGS) -mp=bind
LAPACK=-lblas
