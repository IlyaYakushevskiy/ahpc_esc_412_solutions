### Task 7 

## Notes 

I'm gonna run this task on eiger, in pdf of task6 there are instructions how to run on eiger. tests are done with : (perhaps later I'll need slurm jobfile )


```
uenv start --view=modules,spack prgenv-gnu/25.6:v2
module load cray-mpich boost
module load fftw/3.3.10
./assign $SCRATCH/B100.00100 100 PCS

UPD: now I simply run : sbatch run.job 
```

OpenMP test run with B1000: 

iyakushe@eiger-ln001:~/ahpc_esc_412_solutions/Task07> ./assign $SCRATCH/B1000.00100 100 PCS
Loading 1000000000 particles
Reading file took 10.8694878 s
Mass assignment took 1.4989183 s
FFT and Binning took 2.1743073 s
power spectrum saved to power_log_NGP_100.txt
Projection took 0.0002507 s

Also this morning it takes very long to access $SCRATCH so probably it's eiger issue. Reading should be a bit faster and I tested with last assignments code.  

## Exercise 1 

MPI_Init_thread added 
## Exercise 2

done

created run.job

now one can see outputs in Task07/job_outputs/assign_mpi.log, which I make deliberlately visible 

## Exercise 3

done
also changed assignMassToGridImpl()

## Exercise 4 

done 
using MPI_IN_PLACE

## Exercise 5
done 
now using fftw3-mpi.h

## Exercise 6
done 
## Exercise 7
done 
## Exercise 8
done 
## Exercise 9
done 
## Exercise 10
done 

however I couldn't debug segmentation fault 