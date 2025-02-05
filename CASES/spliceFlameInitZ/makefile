#include /usr/local/include/cantera/Cantera.mak
include $(LEMCANTERA)/Cantera.mak

EXE	= spliceFlameZ
SRC	= spliceFlame.cpp
LIB	= -lLEM -lyaml-cpp -lmpi -lmpi_cxx

CXX	 = g++
CXXFLAGS  = -Wall -pedantic -g -O2 -std=c++11
debugFlag = -Wall -pedantic -g -ggdb -O0 -std=c++11

INCL2	= -I/home/abhi/lesLem2_0/trunk/lean-LEM-for-combustion \
	  $(CANTERA_INCLUDES) \
	  -I/home/abhi/dlignell-odt2.0/source/yaml/include \
	  -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include

LDLIB2	= -L/home/abhi/lesLem2_0/trunk/lean-LEM-for-combustion \
		  -Wl,-rpath=$(LEMLIB)\
	  -L/usr/include/sundials \
	  -L/home/abhi/dlignell-odt2.0/source/yaml/lib -L/usr/lib/x86_64-linux-gnu/openmpi/lib  -L/usr//lib


OBJS	= spliceFlame.o

all: $(EXE)

$(EXE): $(OBJS) makefile
	@echo "\n...making executable \n"
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(EXE) $(INCL2) $(LDLIB2) -lm $(LIB) $(CANTERA_LIBS)


$(OBJS): $(SRC)
	@echo "\n...making object \n"
	$(CXX) $(CXXFLAGS)  -c $(SRC) $(INCL2)

debug: CXXFLAGS = $(debugFlag)
debug: all

clean:
	rm -f $(OBJS) $(EXE)
Clear:
	rm -f *.dat  *.png *.gif *.pdf
