function npoints_dependence(nsamples = 5,walkers = 1000) 
    println("The number of workers to be used is ",nworkers()); 
    D = 1.0; 
    xi = 0.0; 
    xf = 1.0; 
    
    tf = 0.1; 
    
    npoints = 100; 
    times = Array(Float64,4,nsamples); 
    errors = Array(Float64,4,nsamples); 
    
    for n in 1:nsamples 
        println("Starting with size of ",npoints); 
        times[1,n] = @elapsed ((e_sol, e_x) = explicit_solver(u0_test,xi,xf,tf,D,npoints));
        times[2,n] = @elapsed (c_sol = crank_solver(u0_test,xi,xf,tf,D,npoints,200)[1]);
        times[3,n] = @elapsed (m_mc_sol = map_MC_crank_solver(u0_test,xi,xf,tf,D,npoints,walkers)[1]);
        times[4,n] = @elapsed (p_mc_sol = pmap_MC_crank_solver(u0_test,xi,xf,tf,D,npoints,walkers)[1]);
        
        exact_sol = map(a->exact_solution_test(a,tf),e_x); 
        
        errors[1,n] = maximum(abs(exact_sol .- e_sol)); 
        errors[2,n] = maximum(abs(exact_sol .- c_sol)); 
        errors[3,n] = maximum(abs(exact_sol .- m_mc_sol)); 
        errors[4,n] = maximum(abs(exact_sol .- p_mc_sol)); 
        
        println("times: \n",times[1:4,n]'); 
        println("errors: \n ",errors[1:4,n]'); 
        
        npoints += 100; 
    end
end

npoints_dependence(5,5000);
