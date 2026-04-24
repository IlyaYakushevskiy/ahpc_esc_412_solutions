// This uses features from C++17, so you may have to turn this on to compile
#include <iostream>
#include <fstream>
#include <cstdint>
#include <chrono>
#include <locale>
#include <new>
#include <complex>
#include <stdlib.h>
#include "mpi.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "blitz/array.h"
#include "fftw3-mpi.h"
#include "tipsy.h"
#include "aweights.h"
#include <algorithm>
#include <cmath>
#include <numeric>

using namespace blitz;
using hrc = std::chrono::high_resolution_clock;
using duration = std::chrono::duration<double>;

// A separate version is created for each different "Order".
// This allows the compiler to optimize the process for each of the four orders
template<int Order=1>
void assign_mass(Array<float,3> &grid, Array<float,2> &R,Array<float,1> &M) {
    auto nGrid = grid.rows(); // total number of particles
    // C++ Lambda to apply the periodic wrap of the grid index
    auto wrap = [nGrid](int i) {
        if (i<0) i+=nGrid;
        else if (i>=nGrid) i-=nGrid;
        return i;
    };
    #pragma omp parallel for
    for(int pn=R.lbound(0); pn<=R.ubound(0); ++pn) {
        float x = R(pn,0);
        float y = R(pn,1);
        float z = R(pn,2);
        float m = M(pn);
        AssignmentWeights<Order,float> Hx((x+0.5f)*nGrid),Hy((y+0.5f)*nGrid),Hz((z+0.5f)*nGrid);
        for(auto i=0; i<Order; ++i) {
            for(auto j=0; j<Order; ++j) {
                for(auto k=0; k<Order; ++k) {
                    #pragma omp atomic
                    grid(wrap(Hx.i+i),wrap(Hy.i+j),wrap(Hz.i+k)) += m * Hx.H[i]*Hy.H[j]*Hz.H[k];
                }
            }
        }
    }
}

int main(int argc, char *argv[]) {
    int thread_support;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&thread_support);
    if (thread_support < MPI_THREAD_FUNNELED) {
        cerr << "Insufficient MPI thread support -- Funneled required" << std::endl;
        return 1;
    }
    int irank, nrank; 
    MPI_Comm_rank(MPI_COMM_WORLD,&irank);
    MPI_Comm_size(MPI_COMM_WORLD,&nrank);

    std::locale::global(std::locale("")); 
    std::cerr.imbue(std::locale()); 
    if (argc<=2) {
        std::cerr << "Usage: " << argv[0] << " tipsyfile.std grid-size[/[L]bin-count] [order]" << std::endl;
        return 1;
    }

#ifdef _OPENMP
    fftwf_init_threads();
#endif

    int nGrid = atoi(argv[2]);
    const auto iNyquist = nGrid/2;
    int nBins = iNyquist;
    bool bLog = false;
    auto p = strchr(argv[2],'/');
    if (p) {
        *p++ = '\0';
        if (*p == 'L') {
            bLog = true;
            ++p;
        }
        nBins = atoi(p);
    }

    int iOrder = 1;
    if (argc>3) iOrder = atoi(argv[3]);

    
    auto k_nz = nGrid/2 + 1;
    ptrdiff_t local_n0, local_0_start, complex_count;
    complex_count = fftwf_mpi_local_size_3d(nGrid,nGrid,k_nz,MPI_COMM_WORLD,&local_n0,&local_0_start);

    int start = local_0_start;
    Array<int,1> rank_start(nrank);
    MPI_Allgather(&start,1,MPI_INT,rank_start.data(),1,MPI_INT,MPI_COMM_WORLD);

    int count = local_n0;
    Array<int,1> rank_count(nrank);
    MPI_Allgather(&count,1,MPI_INT,rank_count.data(),1,MPI_INT,MPI_COMM_WORLD);

    // Task 8.1: Create map from slab to rank
    Array<int, 1> slab_to_rank(nGrid);
    for (int r = 0; r < nrank; ++r) { 
        for (int i = rank_start(r); i < rank_start(r) + rank_count(r); ++i) {
            slab_to_rank(i) = r; 
        }
    }

    // reading file 
    auto t0 = hrc::now();
    std::ifstream io(argv[1],std::ifstream::binary);
    if (!io) {
        std::cerr << "Unable to open tipsy file " << argv[1] << std::endl;
        return errno;
    }
    constexpr std::streamsize buffer_size = 8 * 1024 * 1024;
    char* buffer = new char[buffer_size];
    io.rdbuf()->pubsetbuf(buffer, buffer_size);

    tipsy::header h;
    if (!io.read(reinterpret_cast<char*>(&h),sizeof(h))) {
        std::cerr << "error reading header" << std::endl;
        return errno;
    }

    std::uint64_t N = h.nDark;
    if (irank==0) std::cerr << "Loading " << N << " particles" << std::endl;

    auto nper = (N + nrank - 1) / nrank;
    auto beg = nper * irank;
    auto end = nper * (irank + 1);
    if (beg>N) beg = N;
    if (end>N) end = N;
    auto local_count = end - beg; 

    // Task 8.1Combine r and m arrays
    Array<float, 2> r_m(local_count, 4);
    Array<float, 2> r = r_m(Range::all(), Range(0, 2));
    Array<float, 1> m = r_m(Range::all(), 3);

    io.seekg( sizeof(tipsy::header) + beg*sizeof(tipsy::dark));
    tipsy::dark d;
    for(int i=0; i<local_count; ++i) {
        if (!io.read(reinterpret_cast<char*>(&d),sizeof(d))) {
            perror(argv[1]); abort();
        }
        r(i,0) = d.pos[0];
        r(i,1) = d.pos[1];
        r(i,2) = d.pos[2];
        m(i) = d.mass;
    }
    delete [] buffer;

    duration dt = hrc::now() - t0;
    auto rate = (sizeof(tipsy::header) + N*sizeof(tipsy::dark)) / dt.count() / 1024 / 1024;
    if (irank==0) std::cerr << "File reading took " << std::setw(9) << dt.count() << " seconds (" << rate <<" MB/s)." << std::endl;

    //Task 8.2: Sort particles 
    struct ParticleData {
        float x, y, z, m;
    };

    ParticleData* p_data = reinterpret_cast<ParticleData*>(r_m.data());

    auto get_target_rank = [nGrid, &slab_to_rank](float x) {
        int start_slab = static_cast<int>(std::floor(x * nGrid)); 
        if (start_slab < 0) start_slab += nGrid;
        else if (start_slab >= nGrid) start_slab -= nGrid;
        return slab_to_rank(start_slab);
    };

    auto compare_particles = [&get_target_rank](const ParticleData& p1, const ParticleData& p2) {
        return get_target_rank(p1.x) < get_target_rank(p2.x); 
    };

    std::sort(p_data, p_data + local_count, compare_particles);

    //task 8.4 

    Array<int, 1> scounts(nrank);
    scounts = 0; 

    for (int i = 0; i < local_count; ++i) {
        int target = get_target_rank(p_data[i].x);
        scounts(target)++;
    }
    if (irank == 0) {
        std::cerr << "scounts: ";
        for (int r = 0; r < nrank; ++r) {
            std::cerr << scounts(r) << " ";
        }
        std::cerr << std::endl;
    }

    //task 8.4

    Array<int, 1> soffset(nrank);
    std::exclusive_scan(scounts.data(), scounts.data() + nrank, soffset.data(), 0);
    // if (irank == 0) {
    //     std::cerr << "rank 0 scounts: ";
    //     for (int r = 0; r < nrank; ++r) std::cerr << scounts(r) << " ";
    //     std::cerr << "rank 0 soffset: ";
    //     for (int r = 0; r < nrank; ++r) std::cerr << soffset(r) << " ";
    //     std::cerr << std::endl;
    // }

    //task 8.5
    for (int r = 0; r < nrank; ++r) {
        for (int i = 0; i < scounts(r); ++i) {
            int p_idx = soffset(r) + i;
            int actual_target = get_target_rank(p_data[p_idx].x);
            
            if (actual_target != r) {
                std::cerr << " error on Rank " << irank 
                          << ": Particle at memory index " << p_idx 
                          << " belongs to rank " << actual_target 
                          << ", but is sorted into the offset block for rank " << r << "!" 
                          << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
    }

    if (irank == 0) {
        std::cerr << "particle sorting and memory offsets verified successfully." << std::endl;
    }

    //task 8.6

    Array<int, 1> rcounts(nrank);

    
    MPI_Alltoall(scounts.data(), 1, MPI_INT, 
                 rcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    Array<int, 1> roffset(nrank);

    std::exclusive_scan(rcounts.data(), rcounts.data() + nrank, roffset.data(), 0);

    //task 8.7 
    int total_recv_count = blitz::sum(rcounts);

    Array<float, 2> r_m2(total_recv_count, 4);
    Array<float, 2> r2 = r_m2(Range::all(), Range(0, 2));
    Array<float, 1> m2 = r_m2(Range::all(), 3);
    Array<int, 1> scounts_float(nrank), soffset_float(nrank);
    Array<int, 1> rcounts_float(nrank), roffset_float(nrank);
    
    for(int r = 0; r < nrank; ++r) {
        scounts_float(r) = scounts(r) * 4;
        soffset_float(r) = soffset(r) * 4;
        rcounts_float(r) = rcounts(r) * 4;
        roffset_float(r) = roffset(r) * 4;
    }
    MPI_Alltoallv(r_m.data(), scounts_float.data(), soffset_float.data(), MPI_FLOAT,
                  r_m2.data(), rcounts_float.data(), roffset_float.data(), MPI_FLOAT,
                  MPI_COMM_WORLD);
                  
    if (irank == 0) {
        std::cerr << "Particle exchange complete." << std::endl;
    }

    
   // Mass assignment
    t0 = hrc::now();

    auto n_floats = size_t(1) * nGrid * nGrid * 2*k_nz; 
    float *data = new (std:: align_val_t (64)) float [n_floats]; 
    Array<float,3> raw_grid(data,shape(nGrid,nGrid,2*k_nz),deleteDataWhenDone);
    Array<float,3> grid(raw_grid(Range(0,nGrid-1),Range(0,nGrid-1),Range(0,nGrid-1)));
    Array <std::complex<float>,3> kgrid(reinterpret_cast<std::complex<float>*>(data),shape(nGrid,nGrid,k_nz),neverDeleteData);

    raw_grid = 0;

    if (irank==0) std::cerr << "Assigning mass to the grid using order " << iOrder <<std::endl;
    switch(iOrder) {
        case 1: assign_mass<1>(grid,r,m); break;
        case 2: assign_mass<2>(grid,r,m); break;
        case 3: assign_mass<3>(grid,r,m); break;
        case 4: assign_mass<4>(grid,r,m); break;
        default: std::cerr << "Invalid order " << iOrder << std::endl;
    }
    
    dt = hrc::now() - t0;
    if (irank==0) std::cerr << "Mass assignment took " << std::setw(9) << dt.count() << " seconds." << std::endl;
    if (irank==0) std::cerr << "Total mass assigned is " << std::setw(9) << blitz::sum(grid) << std::endl;

    // FFT and binning 
    float *slab_data = new (std:: align_val_t (64)) float [2*complex_count]; 
    Array<float,3> raw_slab(slab_data,shape(local_n0,nGrid,2*k_nz),deleteDataWhenDone);
    raw_slab.reindexSelf(TinyVector<int,3>(local_0_start,0,0)); 
    Array<float,3> slab(raw_slab(Range(local_0_start,local_0_start+local_n0-1),Range(0,nGrid-1),Range(0,nGrid-1)));

    Array <std::complex<float>,3> kslab(reinterpret_cast<std::complex<float>*>(slab_data),shape(local_n0,nGrid,k_nz),neverDeleteData);
    kslab.reindexSelf(TinyVector<int,3>(local_0_start,0,0)); 

#ifdef _OPENMP
    fftwf_plan_with_nthreads(omp_get_max_threads());
#endif
    auto plan = fftwf_mpi_plan_dft_r2c_3d(nGrid,nGrid,nGrid,
            slab.dataFirst(),
            reinterpret_cast<fftwf_complex*>(kslab.dataFirst()),
            MPI_COMM_WORLD,FFTW_ESTIMATE);

    int slab_size = nGrid * 2*k_nz;
    for(auto root=0; root<nrank; ++root) {
        MPI_Reduce(&grid(rank_start(root),0,0),slab.data(),slab_size * rank_count(root),
            MPI_FLOAT, MPI_SUM, root, MPI_COMM_WORLD);
    }

    float local_mass = blitz::sum(slab);
    float total_mass;
    MPI_Allreduce(&local_mass, &total_mass, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    float diRhoBar = ((1.0f*nGrid*nGrid*nGrid)/total_mass);
    slab = slab * diRhoBar - 1.0f;
    slab /= (nGrid*nGrid*nGrid);

    t0 = hrc::now();
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    dt = hrc::now() - t0;
    if (irank==0) std::cerr << "FFT took " << std::setw(9) << dt.count() << " seconds." << std::endl;

    Array<double,1> ak(nBins);  
    Array<double,1> pk(nBins);  
    Array<long,1> nk(nBins);  
    ak = 0; pk = 0; nk = 0;

    if (irank==0) std::cerr << "Using " << nBins << (bLog?" logarithmic":" linear") << " bins." << std::endl;

    AssignmentWindow W(nGrid,iOrder);
    for(auto ii=kslab.begin(); ii!=kslab.end(); ++ii) {
        auto pos = ii.position();
        auto bin = [iNyquist,nGrid](int k) {return k<=iNyquist ? k : k-nGrid;};
        auto kx = bin(pos[0]);
        auto ky = bin(pos[1]);
        auto kz = pos[2];
        float k = sqrt(kx*kx + ky*ky + kz*kz);
        *ii *= W[std::abs(kx)] * W[std::abs(ky)] * W[kz]; 
        int i;
        if (bLog) i = int(log(k) / log(iNyquist) * nBins);
        else  i = int(k / iNyquist * nBins);
        if (i>0 && i < nBins) {
            ak(i) += k;
            pk(i) += std::norm(*ii);
            nk(i) += 1;
        }
    }

    Array<double,1> sum_ak(nBins);
    Array<double,1> sum_pk(nBins);
    Array<long,1> sum_nk(nBins);

    MPI_Reduce(ak.data(), sum_ak.data(), nBins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(pk.data(), sum_pk.data(), nBins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(nk.data(), sum_nk.data(), nBins, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (irank==0)  {
        sum_ak = sum_ak / sum_nk;
        sum_pk = sum_pk / sum_nk;
        for(auto i=1; i<nBins; ++i) {
            if (sum_nk(i)) {
                printf("%.10g %.10g %ld\n", sum_ak(i), sum_pk(i),sum_nk(i));
            }
        }
    }

    MPI_Finalize();
}