# Exercise 3 

## Task 1 

Reading file took 0.0426287 s
Mass assignment took 0.0149476 s
Projection took 0.00899075 s

## Task 2 

here we expect algebraic scaling ofc, incidentally the times I get are not deterministic, for some reason my system gives different amount of resourses then, if I were to compare times properly I should compare avg. 

algebraic scaling always implies log-log plot is usefull for linear relationships, the lin-log would be used for exponential scaling 

100^3 : 

Loading 1000000 particles
Reading file took 0.0285946 s
Mass assignment took 0.00291283 s
Projection took 0.000845333 s

200^3 : 

Loading 8000000 particles
Reading file took 0.186947 s
Mass assignment took 0.0428932 s
Projection took 0.00829096 s

500^3 : 

Loading 125000000 particles
Reading file took 3.21436 s
Mass assignment took 1.07015 s
Projection took 0.145411 s

![alt text](ex03(ex02cp)/time-plots-task3-2.png)

## Task 3 

I defined the interpolation schemes as cpp lambda functions ins assign.cpp . Now we can specify scheme by passing 3rd argument as a string e.g. "./assign B100.00100 CIC"


## Task 4 

Great success! Although reading file is a slower, mass assignment loop is faster. Unfortunately, I my initial universal loop couln't be paralelised so I had to rewrite the loop for each scheme. 
Loading 125000000 particles
Reading file took 3.52853 s
Mass assignment took 0.861399 s
Projection took 0.15141 s
saved to projection.dat


--------------------
my notes: 
Cubic splines come from piecewise interpolation by polynimials of 4th order (deg = 3 ), incidentally looks like RK4 but ain't the one. Each particle is shaped as B spline, for density estimation we sum up contribution of each partcile in the cell. 