// This uses features from C++17, so you may have to turn this on to compile
#include <fstream>
#include <cstdint>
#include <stdlib.h>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <exception>
#include <string>
#include "tipsy.h"
#include "assign.h"
#include <new>

std::chrono::high_resolution_clock::time_point getTime() {
	return std::chrono::high_resolution_clock::now();
}

float getDeltaTime(std::chrono::high_resolution_clock::time_point startTime,
                   std::chrono::high_resolution_clock::time_point endTime) {
	return std::chrono::duration<float>(endTime - startTime).count();
}

void assignMassToGrid(Array<float,3>& grid, Array<float,2>& r, Array<float,1>& m, int N, int nGrid, float x_min, float x_max, float y_min, float y_max, float z_min, float z_max, MAS scheme) {
	switch (scheme) {
		case MAS::NGP:
			assignMassToGridImpl<MAS::NGP>(grid, r, m, N, nGrid, x_min, x_max, y_min, y_max, z_min, z_max);
			break;

		case MAS::CIC:
			assignMassToGridImpl<MAS::CIC>(grid, r, m, N, nGrid, x_min, x_max, y_min, y_max, z_min, z_max);
			break;

		case MAS::TSC:
			assignMassToGridImpl<MAS::TSC>(grid, r, m, N, nGrid, x_min, x_max, y_min, y_max, z_min, z_max);
			break;

		case MAS::PCS:
			assignMassToGridImpl<MAS::PCS>(grid, r, m, N, nGrid, x_min, x_max, y_min, y_max, z_min, z_max);
			break;
	}
}

void runAllTests() {
    std::cout << "Tests not implemented yet." << std::endl;
}

int main(int argc, char *argv[]) {
	std::chrono::high_resolution_clock::time_point startTime;
	std::chrono::high_resolution_clock::time_point endTime;
	float deltaTime = 0.0f;

	if (argc == 2 && std::string(argv[1]) == "--test") {
		try {
			runAllTests();
			return 0;
		} catch (const std::exception& e) {
			std::cerr << "Test run failed: " << e.what() << std::endl;
			return 2;
		}
	}

	if (argc <= 3) {
		std::cerr << "Usage: " << argv[0]
		          << " tipsyfile.std [grid-size] [mass-assignment-scheme: 0|NGP, 1|CIC, 2|TSC, 3|PCS]\n"
		          << "   or: " << argv[0] << " --test"
		          << std::endl;
		return 1;
	}

	int nGrid = atoi(argv[2]);
	MAS scheme = (MAS)atoi(argv[3]); // asci to integer

	startTime = getTime();
	std::ifstream io(argv[1], std::ifstream::binary);
	if (!io) {
		std::cerr << "Unable to open tipsy file " << argv[1] << std::endl;
		return errno;
	}

	tipsy::header h;
	if (!io.read(reinterpret_cast<char*>(&h), sizeof(h))) {
		std::cerr << "error reading header" << std::endl;
		return errno;
	}
	tipsy::swap(h);

	std::uint64_t N = h.nDark;
	std::cerr << "Loading " << N << " particles" << std::endl;
	Array<float,2> r(N,3);
	Array<float,1> m(N);

	tipsy::dark d;
	float x_min = __FLT_MAX__;
	float y_min = __FLT_MAX__;
	float z_min = __FLT_MAX__;
	float x_max = -__FLT_MAX__;
	float y_max = -__FLT_MAX__;
	float z_max = -__FLT_MAX__;

	for (int i = 0; i < N; ++i) {
		if (!io.read(reinterpret_cast<char*>(&d), sizeof(d))) {
			std::cerr << "error reading particle" << std::endl;
			return errno;
		}
		tipsy::swap(d);

		r(i,0) = d.pos[0];
		r(i,1) = d.pos[1];
		r(i,2) = d.pos[2];

		if (d.pos[0] < x_min) x_min = d.pos[0];
		if (d.pos[1] < y_min) y_min = d.pos[1];
		if (d.pos[2] < z_min) z_min = d.pos[2];
		if (d.pos[0] > x_max) x_max = d.pos[0];
		if (d.pos[1] > y_max) y_max = d.pos[1];
		if (d.pos[2] > z_max) z_max = d.pos[2];

		m(i) = d.mass;
	}

	endTime = getTime();
	deltaTime = getDeltaTime(startTime, endTime);
	std::cout << "Reading file took " << std::fixed << std::setprecision(7) << deltaTime << " s" << std::endl;

    startTime = getTime();
    
    size_t total_cells = static_cast<size_t>(nGrid) * nGrid * (nGrid + 2);
    float *grid_data = new (std::align_val_t(64)) float[total_cells]; //512 bits externally 
    Array<float,3> grid(grid_data, blitz::shape(nGrid, nGrid, nGrid + 2), blitz::deleteDataWhenDone); // defining blitz array with policy 
    grid = 0.0f; 
    assignMassToGrid(grid, r, m, N, nGrid, x_min, x_max, y_min, y_max, z_min, z_max, scheme);

	endTime = getTime();
	deltaTime = getDeltaTime(startTime, endTime);
	std::cout << "Mass assignment took " << std::fixed << std::setprecision(7) << deltaTime << " s" << std::endl;

	startTime = getTime();
	Array<float,2> projected(nGrid,nGrid);
	for (int i = 0; i < nGrid; i++) {
		for (int j = 0; j < nGrid; j++) {
			float max_val = __FLT_MIN__;
			for (int k = 0; k < nGrid; k++) {
				if (grid(i,j,k) > max_val) {
					max_val = grid(i,j,k);
				}
			}
			projected(i,j) = max_val;
		}
	}
	endTime = getTime();
	deltaTime = getDeltaTime(startTime, endTime);
	std::cout << "Projection took " << std::fixed << std::setprecision(7) << deltaTime << " s" << std::endl;

	std::string suffix;
	switch (scheme) {
		case MAS::NGP: suffix = "_NGP"; break;
		case MAS::CIC: suffix = "_CIC"; break;
		case MAS::TSC: suffix = "_TSC"; break;
		case MAS::PCS: suffix = "_PCS"; break;
	}

	const std::string outName =
	    "density" + suffix + "_" + std::to_string(nGrid) + ".bin";
	std::ofstream outFile(outName, std::ios::binary);
	if (!outFile) {
		std::cerr << "Error: Could not open file for writing." << std::endl;
		return 1;
	}

	outFile.write(reinterpret_cast<const char*>(&nGrid), sizeof(nGrid));
	outFile.write(reinterpret_cast<const char*>(projected.data()), nGrid * nGrid * sizeof(float));

	return 0;
}