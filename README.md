# CTR_clearance
Run zerosClearance.m to generate a initial guess, then run main2017.m to solve for model with clearance.
findShape2: takes a curvature array to compute tube centerline shape.
main2017: Implemented based on https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7989797
zeroClearance:  Implemented based on https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=5509311
Current problem: The loop for optimizing J_lambda was not working properly. J_lambda is supposed to increase for each step, but it starts to decrease in the very first step.
