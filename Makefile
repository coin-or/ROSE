# Rose main Makefile
# Author: Leo Liberti
# Source: GNU Make
# License: (C) Leo Liberti, all rights reserved. Published under the 
#          Common Public License
# History: 050215 work started

CXX		= c++
#CXXFLAGS      = -g -DDEBUG -pg -DPROFILER
#CXXFLAGS      = -g -DDEBUG -DMEMDEBUG -DVERBOSEDEBUG 
#CXXFLAGS	 = -g 
CXXFLAGS	 = -O2 
CC		= cc
CFLAGS		= $(CXXFLAGS)
#FC		= gfortran
#FFLAGS		= $(CXXFLAGS)
MAKE		= make
CP		= cp
RM		= rm
EV3LIB		= Ev3/libEv3.a
AMPLLIB         = external/libamplsolver.a
## only used for static linking if libg2c not available
#F2CLIB		= libf2c.o
OBJS		= problem.o utils.o solver.o newsolver.o \
                  solver_tabu.o solver_vns.o \
                  solver_gomory.o solver_limitedbranch.o solver_localbranch.o\
                  solver_analyser.o solver_Rprodbincont.o solver_Rprint.o \
                  solver_Rsmith.o solver_Rcopy.o solver_Rdisaggr.o \
	          solver_Roa.o solver_Rsymmgroup.o solver_Rprintmod.o \
                  solver_Rvinci.o solver_Rporta.o \
                  solver_Rconvexifier.o solver_RQuarticConvex.o \
                  solver_Rprintdat.o solver_Rcdd.o \
		  solver_Rrelaxation.o solver_subgradient.o \
                  solver_Rconvexifiermod.o solver_Rprintncvxdiscr.o
MAIN		= rose.o
EXE		= rose
LIBS		= $(OBJS) $(EV3LIB)
CXXFLAGS        += -DGLPK48
INCLUDES	= -Iinclude -Iamplsolver
#FORTRANLIB      = -lgfortran

all:	writeversion dir $(AMPLLIB) $(EV3LIB) \
        $(OBJS) $(MAIN) cleanexe\
        $(EXE) nl2ros 

writeversion:
	rm -f rose.o
	./mkversion
	./mksolverlist

dir:
	test -d external/ || mkdir external

$(EXE):  $(LIBS) $(MAIN)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(LIBS) $(MAIN) -ldl $(AMPLLIB)

#libf2cobj: $(F2CLIB)

$(EV3LIB): Ev3/
	$(MAKE) -C Ev3 CXXFLAGS="$(CXXFLAGS) -DHAS_STRINGSTREAM" 

$(AMPLLIB):
	./get.ASL
	cd amplsolver/ ; ./configurehere
	$(MAKE) -C amplsolver
	mv -fv amplsolver/amplsolver.a external/libamplsolver.a

nl2ros: $(AMPLLIB) nl2ros.cxx
	$(CXX) $(INCLUDES) $(CXXFLAGS) -o nl2ros nl2ros.cxx -Iamplsolver -ldl $(AMPLLIB)

cleanexe:
	$(RM) -f $(EXE)

clean:	cleanexe
	$(RM) -f $(OBJS) $(MAIN)
	$(RM) -f nl2ros

Ev3clean:
	$(MAKE) -C Ev3 distclean

libf2cclean:
	$(RM) -rf $(F2CLIB)

amplsolverclean:
	$(MAKE) -C amplsolver clean

distclean: clean Ev3clean libf2cclean 
	$(RM) -f *~ core* redcongraph.dot
	$(RM) -rf amplsolver/
	$(RM) -rf external/

%.o:    %.cxx %.h
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) -o $@ $<

%.o:    %.cc
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) -o $@ $<

%.o:    %.c
	$(CC) -c $(CFLAGS) -o $@ $<
