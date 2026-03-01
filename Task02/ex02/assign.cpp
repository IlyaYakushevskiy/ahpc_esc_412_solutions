// This uses features from C++17, so you may have to turn this on to compile
#include <fstream>
#include <cstdint>
#include <stdlib.h>
#include "blitz/array.h"
#include "tipsy.h"
using namespace blitz;

int main(int argc, char *argv[]) {
    if (argc<=1) {
        std::cerr << "Usage: " << argv[0] << " tipsyfile.std [grid-size]"
                  << std::endl;
        return 1;
    }

    int nGrid = 100;
    if (argc>2) nGrid = atoi(argv[2]);

    std::ifstream io(argv[1],std::ifstream::binary);
    if (!io) {
        std::cerr << "Unable to open tipsy file " << argv[1] << std::endl;
        return errno;
    }

    tipsy::header h;
    if (!io.read(reinterpret_cast<char*>(&h),sizeof(h))) {
        std::cerr << "error reading header" << std::endl;
        return errno;
    }

    // ilya: big endian - most sign bit first , small endian (e.g. ARM , modern)- LSB first

    tipsy::swap(h); // Don't forget to write this function in tipsy.h

    // Load particle positions and masses
    std::uint64_t N = h.nDark;
    std::cerr << "Loading " << N << " particles" << std::endl;
    Array<float,2> r(N,3);
    Array<float,1> m(N);

    // Load the particles
    tipsy::dark d;
    for(int i=0; i<N; ++i) {
        if (!io.read(reinterpret_cast<char*>(&d),sizeof(d))) {
            std::cerr << "error reading particle" << std::endl;
            return errno;
        }
        tipsy::swap(d); // Don't forget to write this function in tipsy.h
        // Save the position
        // r(?,?) = 
        // m(?) = 
        r(i, 0) = d.pos[0];
        r(i, 1) = d.pos[1];
        r(i, 2) = d.pos[2];

        m(i) = d.mass;
    }

    // Create Mass Assignment Grid
    Array<float,3> grid(nGrid,nGrid,nGrid);

    grid = 0;
    for(int pn=0; pn<N; ++pn) {
        float x = r(pn,0);
        float y = r(pn,1);
        float z = r(pn,2);

        // Convert x, y and z into a grid position i,j,k such that
        // 0 <= i < nGrid
        // 0 <= j < nGrid
        // 0 <= k < nGrid

        int i = static_cast<int>((x + 0.5) * nGrid);
        int j = static_cast<int>((y + 0.5) * nGrid);
        int k = static_cast<int>((z + 0.5) * nGrid);

        if (i >= nGrid) i = nGrid - 1; if (i < 0) i = 0;
        if (j >= nGrid) j = nGrid - 1; if (j < 0) j = 0;
        if (k >= nGrid) k = nGrid - 1; if (k < 0) k = 0;

        // Deposit the mass onto grid(i,j,k)
        // grid(i,j,k) += ??
        grid(i,j,k) += m(pn);
    }


    // Calculate projected density
    // - create a 2D map and initialize to zero
    // - loop over the 3D grid
    // - if the density of the grid projected onto the map is greater than the current value
    //   - update the current value


    // Write out the 2D map
    // Read into Python and plot
    Array<float, 2> projected(nGrid, nGrid);
    projected = 0;

    for (int i = 0; i < nGrid; ++i) {
        for (int j = 0; j < nGrid; ++j) {
            float max_val = 0;
            for (int k = 0; k < nGrid; ++k) {
                if (grid(i, j, k) > max_val) {
                    max_val = grid(i, j, k);
                }
            }
            projected(i, j) = max_val;
        }
    }

    //save in binary 
    std::ofstream out("projection.dat", std::ios::binary);
    out.write(reinterpret_cast<const char*>(projected.data()), 
              projected.size() * sizeof(float));
    out.close();

    std::cout << "saved to projection.dat" << std::endl;
}
