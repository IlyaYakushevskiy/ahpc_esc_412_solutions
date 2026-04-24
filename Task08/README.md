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

# Excercise 1 

Here we create kinda look up table, I'm marking a places in my code with comments so you can see it under task 8.1 in assign.cpp
# Excercise 2

when we sort r_m, the data in r and m is sorted too , we can keep assign_mass<1>(grid,r,m). 
; changed arrays to 

 Array<float, 1> m = r_m(Range::all(), 3);

 now sorting in mass assignment grid 

 std::sort(p_data, p_data + local_count, compare_particles);

# Excercise 1 
# Excercise 1 
# Excercise 1 
# Excercise 1 
# Excercise 1 
# Excercise 1 