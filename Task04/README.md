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