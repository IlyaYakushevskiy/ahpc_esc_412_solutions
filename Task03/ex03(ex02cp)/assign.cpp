// This uses features from C++17, so you may have to turn this on to compile
#include <fstream>
#include <cstdint>
#include <stdlib.h>
#include "blitz/array.h"
#include "tipsy.h"
#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>


using namespace blitz;

int main(int argc, char *argv[]) {
    if (argc<=1) {
        std::cerr << "Usage: " << argv[0] << " tipsyfile.std [grid-size]"
                  << std::endl;
        return 1;
    }

    int nGrid = 100; //default 
    if (argc>2) nGrid = atoi(argv[2]);

    std::string scheme = "PCS";
    if (argc>3) scheme = argv[3];

    auto start_read = std::chrono::high_resolution_clock::now();
    
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

    auto end_read = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_read = end_read - start_read;
    auto start_mass = std::chrono::high_resolution_clock::now();
    
    Array<float,3> grid(nGrid,nGrid,nGrid);
    grid = 0;


    auto calc_NGP = [](float s) { return (std::abs(s) < 0.5f) ? 1.0f : 0.0f; };
    auto calc_CIC = [](float s) { return (std::abs(s) < 1.0f) ? (1.0f - std::abs(s)) : 0.0f; };
    auto calc_TSC = [](float s) {
        float abs_s = std::abs(s);
        if (abs_s < 0.5f) return 0.75f - abs_s * abs_s;
        if (abs_s < 1.5f) return 0.5f * std::pow(1.5f - abs_s, 2.0f);
        return 0.0f;
    };

    auto calc_PCS = [](float s) {
        float abs_s = std::abs(s);
        if (abs_s < 1.0f) return (1.0f / 6.0f) * (4.0f - 6.0f * abs_s * abs_s + 3.0f * std::pow(abs_s, 3.0f));
        if (abs_s < 2.0f) return (1.0f / 6.0f) * std::pow(2.0f - abs_s, 3.0f);
        return 0.0f;
    };

    
    if (scheme == "NGP") {
        #pragma omp parallel for
        for(int pn=0; pn<N; ++pn) {
            float rx = (r(pn,0) + 0.5f) * nGrid;
            float ry = (r(pn,1) + 0.5f) * nGrid;
            float rz = (r(pn,2) + 0.5f) * nGrid;
            int istart = std::floor(rx - 0.0f);
            int jstart = std::floor(ry - 0.0f);
            int kstart = std::floor(rz - 0.0f);

            for (int i = istart; i < istart + 1; ++i) {
                for (int j = jstart; j < jstart + 1; ++j) {
                    for (int k = kstart; k < kstart + 1; ++k) {
                        float W = calc_NGP((i + 0.5f) - rx) * calc_NGP((j + 0.5f) - ry) * calc_NGP((k + 0.5f) - rz);
                        int iw = (i % nGrid + nGrid) % nGrid, jw = (j % nGrid + nGrid) % nGrid, kw = (k % nGrid + nGrid) % nGrid;
                        #pragma omp atomic
                        grid(iw, jw, kw) += m(pn) * W;
                    }
                }
            }
        }
    } 
    else if (scheme == "CIC") {
        #pragma omp parallel for
        for(int pn=0; pn<N; ++pn) {
            float rx = (r(pn,0) + 0.5f) * nGrid;
            float ry = (r(pn,1) + 0.5f) * nGrid;
            float rz = (r(pn,2) + 0.5f) * nGrid;
            int istart = std::floor(rx - 0.5f);
            int jstart = std::floor(ry - 0.5f);
            int kstart = std::floor(rz - 0.5f);

            for (int i = istart; i < istart + 2; ++i) {
                for (int j = jstart; j < jstart + 2; ++j) {
                    for (int k = kstart; k < kstart + 2; ++k) {
                        float W = calc_CIC((i + 0.5f) - rx) * calc_CIC((j + 0.5f) - ry) * calc_CIC((k + 0.5f) - rz);
                        int iw = (i % nGrid + nGrid) % nGrid, jw = (j % nGrid + nGrid) % nGrid, kw = (k % nGrid + nGrid) % nGrid;
                        #pragma omp atomic
                        grid(iw, jw, kw) += m(pn) * W;
                    }
                }
            }
        }
    }
    else if (scheme == "TSC") {
        #pragma omp parallel for
        for(int pn=0; pn<N; ++pn) {
            float rx = (r(pn,0) + 0.5f) * nGrid;
            float ry = (r(pn,1) + 0.5f) * nGrid;
            float rz = (r(pn,2) + 0.5f) * nGrid;
            int istart = std::floor(rx - 1.0f);
            int jstart = std::floor(ry - 1.0f);
            int kstart = std::floor(rz - 1.0f);

            for (int i = istart; i < istart + 3; ++i) {
                for (int j = jstart; j < jstart + 3; ++j) {
                    for (int k = kstart; k < kstart + 3; ++k) {
                        float W = calc_TSC((i + 0.5f) - rx) * calc_TSC((j + 0.5f) - ry) * calc_TSC((k + 0.5f) - rz);
                        int iw = (i % nGrid + nGrid) % nGrid, jw = (j % nGrid + nGrid) % nGrid, kw = (k % nGrid + nGrid) % nGrid;
                        #pragma omp atomic
                        grid(iw, jw, kw) += m(pn) * W;
                    }
                }
            }
        }
    }
    else if (scheme == "PCS") {
        #pragma omp parallel for
        for(int pn=0; pn<N; ++pn) {
            float rx = (r(pn,0) + 0.5f) * nGrid;
            float ry = (r(pn,1) + 0.5f) * nGrid;
            float rz = (r(pn,2) + 0.5f) * nGrid;
            int istart = std::floor(rx - 1.5f);
            int jstart = std::floor(ry - 1.5f);
            int kstart = std::floor(rz - 1.5f);

            for (int i = istart; i < istart + 4; ++i) {
                for (int j = jstart; j < jstart + 4; ++j) {
                    for (int k = kstart; k < kstart + 4; ++k) {
                        float W = calc_PCS((i + 0.5f) - rx) * calc_PCS((j + 0.5f) - ry) * calc_PCS((k + 0.5f) - rz);
                        int iw = (i % nGrid + nGrid) % nGrid, jw = (j % nGrid + nGrid) % nGrid, kw = (k % nGrid + nGrid) % nGrid;
                        #pragma omp atomic
                        grid(iw, jw, kw) += m(pn) * W;
                    }
                }
            }
        }
    }
    else {
        std::cerr << "wrong scheme" << std::endl;
        return 1;
    }


    auto end_mass = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_mass = end_mass - start_mass;


    auto start_proj = std::chrono::high_resolution_clock::now();
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

    auto end_proj = std::chrono::high_resolution_clock::now();

    //task 3.1
    std::chrono::duration<double> time_proj = end_proj - start_proj;

    std::cout << "Reading file took " << time_read.count() << " s\n";
    std::cout << "Mass assignment took " << time_mass.count() << " s\n";
    std::cout << "Projection took " << time_proj.count() << " s\n";

    //save in binary 
    std::ofstream out("projection.dat", std::ios::binary);
    out.write(reinterpret_cast<const char*>(projected.data()), 
              projected.size() * sizeof(float));
    out.close();

    std::cout << "saved to projection.dat" << std::endl;
}
