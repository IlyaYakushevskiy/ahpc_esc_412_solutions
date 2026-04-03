## Task 6 

### Excercise 1

following the guidelines, --enable-openmp, 

    I will run with omp and without: 

- 100 particles without OMP
    Loading 1000000 particles
    Reading file took 0.0425707 s
    Mass assignment took 0.0104646 s
    FFT and Binning took 0.0094540 s
    power spectrum saved to power_log_NGP_100.txt
    Projection took 0.0010447 s

    With OMP: 

    Loading 1000000 particles
    Reading file took 0.0352381 s
    Mass assignment took 0.0115774 s
    FFT and Binning took 0.0079458 s
    power spectrum saved to power_log_NGP_100.txt
    Projection took 0.0009310 s

- On 500 particles (over several runs)

without 
FFT and Binning took 0.0092483 s
FFT and Binning took 0.0122685 s
FFT and Binning took 0.0094190 s

FFT and Binning took 0.0067426 s
FFT and Binning took 0.0071324 s
FFT and Binning took 0.0064683 s

so multithreading is not really big of a change so far 
Perhaps would have cleaner numbers if I separated serial binning part.
### Exercise 2 