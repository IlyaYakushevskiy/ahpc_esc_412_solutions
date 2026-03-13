CXX = g++

BLITZ_INCLUDE=-I /opt/homebrew/include
BOOST_INCLUDE=-I /opt/homebrew/include
OMP_INCLUDE=-I /opt/homebrew/opt/libomp/include

CXXFLAGS = $(BLITZ_INCLUDE) $(BOOST_INCLUDE) $(OMP_INCLUDE) -Wno-deprecated-declarations -O3 -std=c++17 -Xpreprocessor -fopenmp

LDFLAGS = -L /opt/homebrew/opt/libomp/lib
LDLIBS = -lomp

assign: assign.o
	$(CXX) -o assign assign.o $(LDFLAGS) $(LDLIBS)

assign.o: assign.cpp tipsy.h

clean:
	rm -f assign assign.o