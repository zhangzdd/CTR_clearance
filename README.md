# CTR_clearance
Please run zeroClearance.m to generate a initial guess, then run main2017.m to solve for model with clearance (with Dr. Ha model), or run ellip_sol_new (with Zhouyu model). The zeroClearance_pierre.m version used ode45 for solving the IVP, thus making the customization of Ni (discretization coefficient) very inconvenient.

## findShape2: 
Takes a curvature array to compute tube centerline shape.

## main2017: 
Implemented based on https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7989797

## zeroClearance:  
Implemented based on https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=5509311
