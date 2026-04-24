CXX = mpic++


BLITZ_INCLUDE = -I$(HOME)/ESC412/tools/include
BOOST_INCLUDE = -I$(BOOST_ROOT)/include

FFTW_DIR = /user-environment/linux-zen2/fftw-3.3.10-zlmokkafjsvo5lcy5fazu2ssriey3adk
FFTW_INCLUDE  = -I$(FFTW_DIR)/include


OMP = -fopenmp
CXXFLAGS = $(BLITZ_INCLUDE) $(BOOST_INCLUDE) $(FFTW_INCLUDE) $(OMP) -Wno-deprecated-declarations -O3 -std=c++17


LDFLAGS = -L$(FFTW_DIR)/lib $(OMP)
LDLIBS = -lfftw3f_mpi -lfftw3f_threads -lfftw3f


assign : assign.o tipsy.h Makefile
	$(CXX) $(CXXFLAGS) -o assign assign.o $(LDFLAGS) $(LDLIBS)

assign.o : assign.cpp tipsy.h aweights.h Makefile
	$(CXX) $(CXXFLAGS) -c assign.cpp

clean:
	rm -f assign assign.o