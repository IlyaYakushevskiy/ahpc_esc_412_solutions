## Task 4

so first of all, I switched to the approach from the solutions, including assign.h file for better modularity and assign.cpp with much better weight handling outside of the looop 

### Exercise 1 
    ```
    size_t total_cells = static_cast<size_t>(nGrid) * nGrid * nGrid;
    float *grid_data = new (std::align_val_t(64)) float[total_cells]; //512 bits externally 
    Array<float,3> grid(grid_data, blitz::shape(nGrid, nGrid, nGrid), blitz::deleteDataWhenDone); // defining blitz array with policy 
    grid = 0.0f; 
    assignMassToGrid(grid, r, m, N, nGrid, x_min, x_max, y_min, y_max, z_min, z_max, scheme);
    ```

### Excercise 2 

modified assignMassToGrid() and grid_data with +2 padding

### Excercise 3 

ranges are inclusive of the index, so we slice from 0 to nGrid - 1 , add 
```
Array<float,3> grid_view = grid_padded(blitz::Range::all(), blitz::Range::all(), blitz::Range(0, nGrid - 1));
```
now project with grid_view 

### Excecise 4 

so we gotta initialize complex grid now for  FFT,  so we cast (float) grid to complex float.

Following the hint we use new policy for new array kdata neverDeleteData to avoid double detete[] call which would lead to SegFault


### Excercise 5 

since fftw_complex is defined as float[2] , we can safely cast it to float 

Loading 1000000 particles
Reading file took 0.0278137 s
Mass assignment took 0.0091642 s
FFT took 0.0064633 s
Projection took 0.0008679 s