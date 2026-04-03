CXX = mpic++


BLITZ_INCLUDE = -I$(HOME)/ESC412/tools/include
BOOST_INCLUDE = -I$(BOOST_ROOT)/include


FFTW_DIR = /user-environment/linux-zen2/fftw-3.3.10-zlmokkafjsvo5lcy5fazu2ssriey3adk
FFTW_INCLUDE  = -I$(FFTW_DIR)/include


LDFLAGS = -L$(FFTW_DIR)/lib
LDLIBS = -lfftw3f_threads -lfftw3f

CXXFLAGS = $(BLITZ_INCLUDE) $(BOOST_INCLUDE) $(FFTW_INCLUDE) -Wno-deprecated-declarations -O3 -std=c++17 -fopenmp

assign: assign.o test.o
	$(CXX) $(CXXFLAGS) -o assign assign.o test.o $(LDFLAGS) $(LDLIBS)

assign.o: assign.cpp tipsy.h assign.h
	$(CXX) $(CXXFLAGS) -c assign.cpp

test.o: test.cpp test.h assign.h tipsy.h
	$(CXX) $(CXXFLAGS) -c test.cpp

clean:
	rm -f assign assign.o test.o