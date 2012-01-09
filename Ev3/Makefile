# Expression v3 main Makefile
# RCS2 "$Id: Makefile,v 1.11 2005/09/08 23:37:13 liberti Exp liberti $"

CXX	 = c++
CXXFLAGS = -O3 -DHAS_STRINGSTREAM
OBJS	= tree.o operand.o expression.o exprauxiliary.o exprexternal.o parser.o fastexpression.o
LIBS	= $(OBJS)
AR	= ar
RM	= rm

all:	libEv3.a

libEv3.a: test
	$(AR) r libEv3.a $(OBJS)

test:	smithmain.cxx $(OBJS)
	$(CXX) $(CXXFLAGS) -o test $(LIBS) smithmain.cxx
#	$(CXX) $(CXXFLAGS) -o test $(LIBS) convexifiermain.cxx

clean:	
	$(RM) -rf $(OBJS) test libEv3.a

distclean:	clean
	$(RM) -rf *~ core

mrproper:	distclean

%.o:	%.cxx %.h
	$(CXX) -c $(CXXFLAGS) -o $@ $<

%.o:	%.cxx expression.h
	$(CXX) -c $(CXXFLAGS) -o $@ $<
