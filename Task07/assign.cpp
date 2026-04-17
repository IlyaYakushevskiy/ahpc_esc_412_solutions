// This uses features from C++17, so you may have to turn this on to compile
#include "tipsy.h"
#include "assign.h"
#include "test.h"
#include <fstream>
#include <cstdint>
#include <stdlib.h>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <exception>
#include <string>
#include <new>
#include <complex>
#include <omp.h> 
#include <vector>
#include <mpi.h>
#include <fftw3-mpi.h>


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



int main(int argc, char *argv[]) {
    fftwf_mpi_init();
    // task 7.1
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    

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

	//task 6.4
	startTime = getTime();
	std::ifstream header_io(argv[1], std::ifstream::binary);
    if (!header_io) {
        std::cerr << "Unable to open tipsy file " << argv[1] << std::endl;
        return errno;
    }

    tipsy::header h;
    if (!header_io.read(reinterpret_cast<char*>(&h), sizeof(h))) { 
        std::cerr << "error reading header" << std::endl;
        return errno;
    }
    tipsy::swap(h);
    header_io.close();

    std::uint64_t N = h.nDark;
    std::cerr << "Loading " << N << " particles" << std::endl;


    //task 7.3

    std::uint64_t chunk_size = N / size;
    std::uint64_t remainder = N % size;

    std::uint64_t local_count = chunk_size + (rank < remainder ? 1 : 0);
    std::uint64_t start_idx = rank * chunk_size + std::min(static_cast<std::uint64_t>(rank), remainder);
    std::uint64_t end_idx = start_idx + local_count - 1;

    blitz::Range particle_range(start_idx, end_idx);
    Array<float,2> r(particle_range, blitz::Range(0, 2));
    Array<float,1> m(particle_range);

    float local_xmin = __FLT_MAX__, local_ymin = __FLT_MAX__, local_zmin = __FLT_MAX__;
    float local_xmax = -__FLT_MAX__, local_ymax = -__FLT_MAX__, local_zmax = -__FLT_MAX__;

    if (local_count > 0) {
        std::ifstream local_io(argv[1], std::ifstream::binary);
        std::streamoff offset = sizeof(tipsy::header) + static_cast<std::streamoff>(start_idx) * sizeof(tipsy::dark);
        local_io.seekg(offset, std::ios::beg);

        tipsy::dark d;

        for (int i = r.lbound(0); i <= r.ubound(0); ++i) {
            local_io.read(reinterpret_cast<char*>(&d), sizeof(d));
            tipsy::swap(d);

            r(i,0) = d.pos[0];
            r(i,1) = d.pos[1];
            r(i,2) = d.pos[2];
            m(i) = d.mass;

            if (d.pos[0] < local_xmin) local_xmin = d.pos[0];
            if (d.pos[1] < local_ymin) local_ymin = d.pos[1];
            if (d.pos[2] < local_zmin) local_zmin = d.pos[2];
            if (d.pos[0] > local_xmax) local_xmax = d.pos[0];
            if (d.pos[1] > local_ymax) local_ymax = d.pos[1];
            if (d.pos[2] > local_zmax) local_zmax = d.pos[2];
        }
        local_io.close();
    }

    //Map Reduce to find global mins and maxs
    float local_mins[3] = {local_xmin, local_ymin, local_zmin};
    float local_maxs[3] = {local_xmax, local_ymax, local_zmax};
    
    float global_mins[3];
    float global_maxs[3];

    MPI_Allreduce(local_mins, global_mins, 3, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(local_maxs, global_maxs, 3, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

    float x_min = global_mins[0];
    float y_min = global_mins[1];
    float z_min = global_mins[2];
    float x_max = global_maxs[0];
    float y_max = global_maxs[1];
    float z_max = global_maxs[2];

    endTime = getTime();
    deltaTime = getDeltaTime(startTime, endTime);
    if (rank == 0) {
        std::cout << "Reading file took " << std::fixed << std::setprecision(7) << deltaTime << " s" << std::endl;
    }

    endTime = getTime();
	


	deltaTime = getDeltaTime(startTime, endTime);
	if (rank == 0) {
        std::cout << "Reading file took " << std::fixed << std::setprecision(7) << deltaTime << " s" << std::endl;
    }

    startTime = getTime();
    
    //task 7

    size_t full_cells = static_cast<size_t>(nGrid) * nGrid * (nGrid + 2);
    float *full_grid_data = new (std::align_val_t(64)) float[full_cells]; 
    Array<float,3> full_grid_padded(full_grid_data, blitz::shape(nGrid, nGrid, nGrid + 2), blitz::deleteDataWhenDone);
    full_grid_padded = 0.0f; 
    Array<float,3> full_grid_view = full_grid_padded(blitz::Range::all(), blitz::Range::all(), blitz::Range(0, nGrid - 1));

    ptrdiff_t local_n0, local_0_start;
    ptrdiff_t alloc_local = fftwf_mpi_local_size_3d(nGrid, nGrid, nGrid / 2 + 1, MPI_COMM_WORLD, &local_n0, &local_0_start);

    float *slab_data = new (std::align_val_t(64)) float[alloc_local * 2];
    int padded_z = 2 * (nGrid / 2 + 1);

    Array<float,3> slab_padded(slab_data, blitz::shape(local_n0, nGrid, padded_z), blitz::deleteDataWhenDone);
    slab_padded = 0.0f; 
    slab_padded.reindexSelf(blitz::TinyVector<int, 3>(local_0_start, 0, 0));

    Array<float,3> slab_view = slab_padded(
        blitz::Range(local_0_start, local_0_start + local_n0 - 1), 
        blitz::Range::all(), 
        blitz::Range(0, nGrid - 1)
    );

    std::complex<float>* complex_data = reinterpret_cast<std::complex<float>*>(slab_data);
    Array<std::complex<float>, 3> kdata(complex_data, 
                                        blitz::shape(local_n0, nGrid, nGrid / 2 + 1), 
                                        blitz::neverDeleteData);
    kdata.reindexSelf(blitz::TinyVector<int, 3>(local_0_start, 0, 0));

    assignMassToGrid(full_grid_view, r, m, N, nGrid, x_min, x_max, y_min, y_max, z_min, z_max, scheme);

    Array<int, 1> rank_start(size);
    Array<int, 1> rank_count(size);
    int start_int = static_cast<int>(local_0_start);
    int count_int = static_cast<int>(local_n0);

    MPI_Allgather(&start_int, 1, MPI_INT, rank_start.data(), 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&count_int, 1, MPI_INT, rank_count.data(), 1, MPI_INT, MPI_COMM_WORLD);
    // task 7.7
    
    int slice_size = nGrid * (nGrid + 2); // Size of one 2D slice (including padding)
    
    for (int r_target = 0; r_target < size; ++r_target) {
        int reduce_count = rank_count(r_target) * slice_size;
        
        // sendbuf points to the start of the target rank's slab inside our full grid
        float* sendbuf = &full_grid_padded(rank_start(r_target), 0, 0); 
        
        // Reduce this specific slab onto the target rank's slab_data
        MPI_Reduce(sendbuf, slab_data, reduce_count, MPI_FLOAT, MPI_SUM, r_target, MPI_COMM_WORLD);
    }

    delete[] full_grid_data;

	

    float local_mass = blitz::sum(slab_view);
    float total_mass = 0.0f;
    MPI_Allreduce(&local_mass, &total_mass, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        
    float n_grid_cube = static_cast<float>(nGrid * nGrid * nGrid);

    
    for (int i = local_0_start; i < local_0_start + local_n0; i++) {
        for (int j = 0; j < nGrid; j++) {
            for (int k = 0; k < nGrid; k++) {
                slab_view(i,j,k) = (slab_view(i,j,k) * n_grid_cube / total_mass) - 1.0f;
            }
        }
    }    

    MPI_Barrier(MPI_COMM_WORLD);
    endTime = getTime();
    deltaTime = getDeltaTime(startTime, endTime);
    if (rank == 0) {
        std::cout << "Mass assignment and distribute took " << std::fixed << std::setprecision(7) << deltaTime << " s" << std::endl;
    }

    startTime = getTime();

    // task 7.9
    if (fftwf_init_threads() == 0) {
        std::cerr << "error: fftwf_init_threads failed." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int nthreads = omp_get_max_threads();
    fftwf_plan_with_nthreads(nthreads);
    
    Array<std::complex<float>, 3> kslab(complex_data, 
                                        blitz::shape(local_n0, nGrid, nGrid / 2 + 1), 
                                        blitz::neverDeleteData);
    kslab.reindexSelf(blitz::TinyVector<int, 3>(local_0_start, 0, 0));

    // use mpi plan 
    fftwf_plan plan = fftwf_mpi_plan_dft_r2c_3d(
        nGrid, nGrid, nGrid,
        slab_data,                                           
        reinterpret_cast<fftwf_complex*>(complex_data),      
        MPI_COMM_WORLD,
        FFTW_ESTIMATE                                        
    );

    fftwf_execute(plan);

    fftwf_destroy_plan(plan);
    fftwf_cleanup_threads();

    // task 7.10
    float k_max = static_cast<float>(nGrid);
    float log_k_max = std::log(k_max);
    int n_bins = 80; 
    

    std::vector<float> local_fPower(n_bins, 0.0f);
    std::vector<int> local_nPower(n_bins, 0);
    std::vector<float> local_k_sum(n_bins, 0.0f); 
    
    int iNyquist = nGrid / 2;

    for(auto ii = kslab.begin(); ii != kslab.end(); ++ii) {
        auto pos = ii.position();
        auto bin = [iNyquist, nGrid](int k) { return k <= iNyquist ? k : k - nGrid; };
        
        auto kx = bin(pos[0]);
        auto ky = bin(pos[1]);
        auto kz = pos[2];

        float mag_k = std::sqrt(kx * kx + ky * ky + kz * kz); 

        if (mag_k > 0.0f) {
            std::complex<float> delta_k = *ii; 
            float power = std::norm(delta_k);

            int bin_idx = static_cast<int>((std::log(mag_k) / log_k_max) * n_bins);

            if (bin_idx >= 0 && bin_idx < n_bins) {
                local_fPower[bin_idx] += power;
                local_k_sum[bin_idx] += mag_k; 
                local_nPower[bin_idx] += 1; 
            }
        }
    }
            
   // reduce binc to rank 0 
    std::vector<float> fPower(n_bins, 0.0f);
    std::vector<int> nPower(n_bins, 0);
    std::vector<float> k_sum(n_bins, 0.0f);

    MPI_Reduce(local_fPower.data(), fPower.data(), n_bins, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_k_sum.data(), k_sum.data(), n_bins, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_nPower.data(), nPower.data(), n_bins, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    endTime = getTime();
    deltaTime = getDeltaTime(startTime, endTime);
    if (rank == 0) {
        std::cout << "FFT and Binning took " << std::fixed << std::setprecision(7) << deltaTime << " s" << std::endl;
    }

   //writing to bin files 
    std::string suffix;
    switch (scheme) {
        case MAS::NGP: suffix = "_NGP"; break;
        case MAS::CIC: suffix = "_CIC"; break;
        case MAS::TSC: suffix = "_TSC"; break;
        case MAS::PCS: suffix = "_PCS"; break;
    }

    if (rank == 0) {
        const std::string pkName = "power_log" + suffix + "_" + std::to_string(nGrid) + ".txt";
        std::ofstream pkFile(pkName);
        if (pkFile) {
            for (int b = 0; b < n_bins; ++b) {
                if (nPower[b] > 0) {
                    float average_power = fPower[b] / static_cast<float>(nPower[b]);
                    float average_k = k_sum[b] / static_cast<float>(nPower[b]);
                    pkFile << average_k << " " << std::scientific << average_power << "\n";
                }
            }
            pkFile.close();
            std::cout << "power spectrum saved to " << pkName << std::endl;
        }
    }

    // distributed projection 
    startTime = getTime();
    
    Array<float,2> local_projected(nGrid, nGrid);
    local_projected = __FLT_MIN__;

    for (int i = local_0_start; i < local_0_start + local_n0; i++) {
        for (int j = 0; j < nGrid; j++) {
            float max_val = __FLT_MIN__;
            for (int k = 0; k < nGrid; k++) {
                if (slab_view(i,j,k) > max_val) {
                    max_val = slab_view(i,j,k);
                }
            }
            local_projected(i,j) = max_val;
        }
    }
    // reducing to get full piccture 
    Array<float,2> global_projected(nGrid, nGrid);
    MPI_Reduce(local_projected.data(), global_projected.data(), nGrid * nGrid, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    endTime = getTime();
    deltaTime = getDeltaTime(startTime, endTime);
    
    if (rank == 0) {
        std::cout << "Projection took " << std::fixed << std::setprecision(7) << deltaTime << " s" << std::endl;

        const std::string outName = "density" + suffix + "_" + std::to_string(nGrid) + ".bin";
        std::ofstream outFile(outName, std::ios::binary);
        if (outFile) {
            outFile.write(reinterpret_cast<const char*>(&nGrid), sizeof(nGrid));
            outFile.write(reinterpret_cast<const char*>(global_projected.data()), nGrid * nGrid * sizeof(float));
        }
    }

    MPI_Finalize();
    return 0;
}