D = 1.0; 
xi = 0.0; 
xf = 1.0; 
npoints = 100; 

tf = 0.05; 

walkers = 60000; 

@time mc_sol, mc_x = MC_crank_solver(u0_test,xi,xf,tf,D,npoints,walkers);
@time m_mc_sol, m_mc_x = map_MC_crank_solver(u0_test,xi,xf,tf,D,npoints,walkers)
@time p_mc_sol, p_mc_x = pmap_MC_crank_solver(u0_test,xi,xf,tf,D,npoints,walkers); 
