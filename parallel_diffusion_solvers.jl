@everywhere function parallel_explicit_stepper!(num_steps::Int,d_u::DArray,this_lambda::Float64)
    local_index = findfirst(d_u.pmap,myid()); 
    local_range = (d_u.indexes[local_index])[1]; 
    local_array = localpart(d_u); 
    len = length(local_array);
    @assert(len == length(local_range)); 
    new_u = Array(Float64,len);
    
    temp = 1.0-2.0*this_lambda; 
    
    left_boundary = false; 
    right_boundary = false; 
    middle = false; 
    
    #figure out where this piece of the array is in the distributed array
    #and how to get neighboring points 
    if local_index == 1
        left_boundary = true; 
    elseif myid() == d_u.pmap[end]
        right_boundary = true; 
    elseif !left_boundary & !right_boundary
        middle = true; 
    else
        println("Error Determining Array Location"); 
    end

    for j in 1:num_steps
    	#@sync begin 
              if middle 
	      	 for i in 1:len 
		     index = local_range[i];
		     new_u[i] = this_lambda*(d_u[index+1] + d_u[index-1]) + temp*d_u[index] 
		 end
              elseif left_boundary 
              	 new_u[1] = this_lambda*d_u[2] + temp*d_u[1];  
       	       	 for i in 2:len
		     index = local_range[i]; 
		     new_u[i] = this_lambda*(d_u[index+1] +d_u[index-1]) + temp*d_u[index];   
		 end
              elseif right_boundary 
       	      	 new_u[len] = this_lambda*d_u[local_range[len-1]] + temp*d_u[local_range[len]];
       	         for i in 1:len-1
		     index = local_range[i];
		     new_u[i] = this_lambda*(d_u[index+1] + d_u[index-1]) + temp*d_u[index];
		 end
              else
		 println("Error with Bounds"); 
              end 
	#end 
        @sync begin  
            localpart(d_u)[:] = new_u;
        end 
    end 
end

function parallel_explicit_solver(u0,this_xi::Float64,this_xf::Float64,this_tf::Float64,this_D::Float64,nx::Int = 20)
    @assert(this_xf > this_xi); 
    @assert(this_tf > 0.0); 
    @assert(this_D > 0.0); 

    x = linspace(this_xi,this_xf,nx+2); 
    dx = x[2] - x[1]; 
    @assert(dx == (xf-xi)/(nx+1));
    M = iround(tf*this_D/(0.49*dx^2)); 
    dt = tf/M; 
    lambda = dt*this_D/dx^2; 
    @assert(lambda < 0.5); 
    u = map(u0,x[2:end-1]); 
    distributed_u = distribute(u);
     r = RemoteRef(); 
    for p in procs(distributed_u) 	
    	r = @spawnat p parallel_explicit_stepper!(M,distributed_u,lambda);
    end
   
    wait(r);  

    u = convert(Array{Float64,1},distributed_u); 
 
    prepend!(u,[0.0]); 
    append!(u,[0.0]); 

    return u, x, map(u0,x);
end

function pmap_MC_crank_solver(u0,xi::Float64,xf::Float64,tf::Float64,D::Float64,nx::Int = 50,num_walkers::Int = 1000)
    @assert(xf>xi); 
    @assert(tf>0.0); 
    @assert(D > 0.0); 

    x = linspace(xi,xf,nx+2); 
    dx = x[2] - x[1]; 
    @assert(dx == (xf-xi)/(nx+1));
    nt = iceil(tf*D/dx^2); 
    dt = tf/nt; 
    lambda = dt*D/dx^2;
    @assert(lambda < 1.0); 
    
    initial_condition = map(u0,x);
    solution = Array(Float64,nx+2);  
    probs = zeros(Float64,nx+2); 
    probability = (1.0-lambda)/(1.0+lambda); 
    solution[1] = 0.0; 
    solution[end] = 0.0; 

    solution[2:end-1] = pmap(a -> sum(initial_condition.*random_walker(a,1,nx+2,nt,num_walkers,probability)),2:nx+1);
    return solution, x; 
end	 
