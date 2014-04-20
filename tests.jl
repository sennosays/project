#Function that tests for NaN's, Inf's, and negative numbers in a solution
#Inputs: 
#	this_sol -- array of Floats which describes the solution to a PDE
#	positive -- bool which determes if we should check if the solution only contians positive elements 
function test_real(this_sol::Array,positive = false)
    len = length(this_sol)
    index = 1; 
    sol_nan = false; 
    sol_inf = false; 
    sol_neg = !positive;
    is_number = 3; 
    
    #determine if the solution is real 
    while is_number == 3 && index <= len
        sol_nan = !isnan(this_sol[index]); 
        sol_inf = !isinf(this_sol[index]); 
        
        if positive 
            sol_neg = (this_sol[index] >= 0.0);
            println(this_sol[index],' ',sol_neg);
        end
        
        index += 1; 
        
        is_number = sol_nan + sol_inf + sol_neg; 
    end
    
    #if the loop ends early determine and print what the problem is 
    #and where it is located in the array 
    if index <= len 
        index -= 1; 
        
        if !sol_nan
            println("Solution at index ",index, " is NaN"); 
        elseif !sol_inf
            println("Solution at index ",index, " is inf"); 
        elseif !sol_neg
            println("Solution at index ",index, " is negative"); 
        else
            println("Unknown error at index ", index); 
        end
    elseif !positive 
        println("Solution is comprised of real numbers"); 
    else
        println("Solution is comprised of real positive numbers"); 
    end
end


#compare two solutions to a specified absolute tolerance 
function compare_solutions(sol_1::Array,sol_2::Array,ab_tol = 1e-5); 
    if size(sol_1) != size(sol_2)
        println("Solutions are not of the same size"); 
    end
    
    ab_errors = abs(sol_1 .- sol_2); 

    max_error, max_error_loc = findmax(ab_errors); 
    
    if(max_error > ab_tol)
        println("Solutions do not match to an absolute tolerance of ",ab_tol);
        println("at the solution location ",max_error_loc); 
    else
        println("Solutions agree to an absolute tolerance of ",ab_tol); 
    end
end

#compare the solutions generated using different methods with the exact solution
#Inputs:
#	tf 	-- final time of the solution 
#	npoints -- number of spatial points in the solution grid 
#	walkers -- number of walkers used in MC solvers
#	tol 	-- absolute tolerance between the approximate and exact solution 
function compare_solvers(tf::Float64,npoints = 20,walkers = 1000,tol::Float64 = 1e-5)
    @assert(tf > 0.0); 
    
    D = 1.0; 
    xi = 0.0; 
    xf = 1.0;
    
    e_sol, e_x = explicit_solver(u0_test,xi,xf,tf,D,npoints);
    i_sol, i_x = implicit_solver(u0_test,xi,xf,tf,D,npoints,200); 
    c_sol, c_x = crank_solver(u0_test,xi,xf,tf,D,npoints,200);
    m_mc_sol, m_mc_x = map_MC_crank_solver(u0_test,xi,xf,tf,D,npoints,walkers);
    p_mc_sol, p_mc_x = pmap_MC_crank_solver(u0_test,xi,xf,tf,D,npoints,walkers);
    
    exact_sol = map(a->exact_solution_test(a,tf),e_x); 
    
    println("Comparing exact solution with explicit solver"); 
    compare_solutions(exact_sol,e_sol,tol);
    
    println("Comparing exact solution with implicit solver"); 
    compare_solutions(exact_sol,i_sol,tol); 
    
    println("Comparing exact solution with Crank-Nicolson solver"); 
    compare_solutions(exact_sol,c_sol,tol);
    
    println("Comparing exact solution with serial Monte Carlo solver"); 
    compare_solutions(exact_sol,m_mc_sol,tol); 
    
    println("Comparing exact solution with parallel Monte Carlo solver"); 
    compare_solutions(exact_sol,p_mc_sol,tol);
    
end
       
