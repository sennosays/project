dt = 0.1; 
dx = 0.1; 
include("implicit_solver.jl")

fp = open("results_01_01.dat");
write(fp,length(n)); 
write(fp,n); 
close(fp); 
