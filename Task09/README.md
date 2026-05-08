## Task 09 

```
uenv start --view=modules,spack prgenv-gnu/25.6:v2
module load cray-mpich boost
module load fftw/3.3.10
./assign $SCRATCH/B100.00100 100 1

sbatch run.job # put the run command here 
```

Last time feedback: 

"You mention that you tested your code and included plots, but unfortunately, I cannot find them. You also do not describe your approach in sufficient detail, which makes it difficult to understand what you have done.

If you are struggling with Git, I encourage you to attend the exercise sessions so we can help you with it.

Furthermore, please make sure to explain your approaches thoroughly. In a formal setting, such as a thesis or a paper, your work could be rejected due to a lack of clarity and documentation."

- Ok.

For clarity I'll be doing 1 commit per Exercise, I write in comments what I do but I'll also comment on it redundatly in this doc 

# Exercise 1 

So we define those ghost cells as artificial extension of local cell to solve the problem of density assignements for  particles on the boundary. 

- How many slabs? 
    - g=iOrder−1

- Should you do this for kslab (the complex view after the FFT)? HINT: Do you need/want a ghost region after the FFT?
    - Ghost cells virtually dont really change mathematical system like padding so we dont need them after FFt 

In this commit I: 
1) Instantiated raw_slab as base physical array. Reindexed its origin left by the ghost radius to explicitly align the local memory indices with the global MPI spatial coordinates.
2) Extracted two distinct spatial views from the base array: ghost which spans the full computational footprint, including boundary overlaps and slab which is strictly bounded to the unpadded interior domain.
3) Gave the same raw memory pointer to construct kslab. Restricted its shape and origin to the strict unpadded interior, as frequency space does not utilize spatial domain overlaps.

I've then implemented a testing code in my main and ran it on euler. 
And tested with srun ./assign $SCRATCH/B100.00100 100 3 
# Exercise 2 

# Exercise 3

# Exercise 4 

# Exercise 5 