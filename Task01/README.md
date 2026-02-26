# HOMEWORK 1 

## Ex1 

Hi! I'll be writing my answers in this file. 

1) - `blitz::` is the namespace
- <float, 3>` specifies the template parameters
-  `(10,8,6)` passes arguments to the constructor 

2) zhis line creates a 3-dimensional array named `data` that holds `float` data types

int main() {
    std::cout << "data" << std::endl;
    return 0;
}

## Ex2 

2^36 particles * 56 bytes  = 3584GB
3584 GB / 64 GB = 56 nodes

## Ex3 

1) god --endian big -j 0 -N 8 -t f8 B100.00100
0000000       1.0000000000000018 // time of the file 
0000010
2) god --endian big -j 8 -N 4 -t d4 B100.00100 //od doesnt work for mac 
0000010     1000000 // amount of particles 
0000014


3) are above

## Ex4

1) 2) god --endian big -j 3596 -N 16 -t f4 B100.00100
0007014         3.2e-07 //mass      0.49921852     -0.49564773      -0.4898032
//the rest are positions coordinates 


## Ex5

