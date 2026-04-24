## Task 8

```
uenv start --view=modules,spack prgenv-gnu/25.6:v2
module load cray-mpich boost
module load fftw/3.3.10
./assign $SCRATCH/B100.00100 100 1

sbatch run.job # put the run command here 
```

From last task feedback : I will try to include more testing of my code this time, including plots etc. e.g. plot power spectrum from .txt files

For this Task I adapted the code from official solutions of Task 7

These are results of code implemented in task 7, I'll keep it as a baseline and reference: 

```
iyakushe@eiger-ln004:~/ahpc_esc_412_solutions/Task07> ./assign $SCRATCH/B100.00100 100 PCS
Loading 1000000 particles
Reading file took 6.2895842 s
Reading file took 6.2896409 s
Mass assignment and distribute took 0.0874053 s
FFT and Binning took 2.0238338 s
power spectrum saved to power_log_NGP_100.txt
Projection took 0.0004315 s
```

# Excercise 1 

Here we create kinda look up table, I'm marking a places in my code with comments so you can see it under task 8.1 in assign.cpp
# Excercise 2

when we sort r_m, the data in r and m is sorted too , we can keep assign_mass<1>(grid,r,m). 
; changed arrays to 

 Array<float, 1> m = r_m(Range::all(), 3);

 now sorting in mass assignment grid 

 std::sort(p_data, p_data + local_count, compare_particles);


```
Loading 1000000 particles
File reading took   6.11412 seconds (5.61525 MB/s).
Assigning mass to the grid using order 1
Mass assignment took 0.0619709 seconds.
Total mass assigned is      0.32
FFT took   2.02606 seconds.
Using 50 linear bins.
1.367171259 2.531723988e-05 17
2.385326409 2.398413419e-05 41
3.401583682 2.519325305e-05 89
4.412521673 2.202985551e-05 129
5.454827237 1.708678763e-05 225
6.434946136 1.477104853e-05 253
7.436757396 1.346303408e-05 393
.
.
.
48.48399907 2.491962632e-06 14869
49.48505401 2.510178456e-06 15641
iyakushe@eiger-ln004:~/ahpc_esc_412_solutions/Task08> 
```

# Excercise 3
added simple loop reusing lambda function;  made rank 0 to print scounts, got scounts: 1000000 , works for now with 1 mpi rank 

# Excercise 4

so this scan is effectively a sum. 
added std::exclusive_scan(scounts.data(), scounts.data() + nrank, soffset.data(), 0);

printing gives us for one rank (will test more rank later): 
rank 0 scounts: 1000000 rank 0 soffset: 0 

# Excercise 5 
successsfully calculated index of particles



# Excercise 6 

we prepped rcount and roffset for future tasks 

# Excercise 7 

completed particle exchange, will test in next exercise 
# Excercise 8 

I must admit, last commit was a bit after deadline. I'm sorry, I got a bit preassured with my bachelor thesis and needed more time for those task. 