CC:=gcc
FC:=gfortran
OPTFLAGS=-O3 
DBGFLAGS=-O0 -ggdb -Wall

CFLAGS=
FFLAGS=-ffree-line-length-0 -ffixed-line-length-0 -fimplicit-none  -fopenmp
FDBG:=-fbounds-check -ffpe-trap=invalid,zero,overflow,denormal \
	-finit-real=SNAN -finit-integer=-1
LDFLAGS:=-fopenmp

LAPACK:=-framework vecLib
