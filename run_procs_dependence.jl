println("The number of workers to be used is ",nworkers()); 

D = 1.0; 
xi = 0.0; 
xf = 1.0; 

tf = 0.1; 

npoints = 100; 
walkers = 7500; 
this_time = @elapsed (p_mc_sol = pmap_MC_crank_solver(u0_test,xi,xf,tf,D,npoints,walkers)[1]);

println("Time required was ",this_time);  
