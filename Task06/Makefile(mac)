CXX = g++

BLITZ_INCLUDE=-I/opt/homebrew/include
BOOST_INCLUDE=-I/opt/homebrew/include
OMP_INCLUDE=-I/opt/homebrew/opt/libomp/include
FFTW_INCLUDE=-I/opt/homebrew/include

CXXFLAGS = $(BLITZ_INCLUDE) $(BOOST_INCLUDE) $(OMP_INCLUDE) $(FFTW_INCLUDE) -Wno-deprecated-declarations -O3 -std=c++17 -Xpreprocessor -fopenmp

LDFLAGS = -L/opt/homebrew/opt/libomp/lib -L/opt/homebrew/lib
LDLIBS = -lomp -lfftw3f_omp -lfftw3f

assign: assign.o test.o
	$(CXX) -o assign assign.o test.o $(LDFLAGS) $(LDLIBS)

assign.o: assign.cpp tipsy.h assign.h
	$(CXX) $(CXXFLAGS) -c assign.cpp

test.o: test.cpp test.h assign.h tipsy.h
	$(CXX) $(CXXFLAGS) -c test.cpp

clean:
	rm -f assign assign.o test.o