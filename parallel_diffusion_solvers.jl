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
       # left_point_index = length((d_u.indexes[local_index-1])[1]); 
    elseif !left_boundary & !right_boundary
        middle = true; 
        #left_point_index = length((d_u.indexes[local_index-1])[1]); 
    else
        println("Error Determining Array Location"); 
    end
   
    for j in 1:num_steps
    	@sync begin 
        if middle 
           new_u = [this_lambda*(d_u[i+1] + d_u[i-1]) + temp*d_u[i] for i in local_range]; 
        elseif left_boundary 
               new_u[1] = this_lambda*d_u[2] + temp*d_u[1]; 
	       #new_u[len] = this_lambda*(d_u[len+1] + d_u[len-1]) + temp*d_u[len]; 
       	       new_u[2:len] = [this_lambda*(d_u[i+1] +d_u[i-1]) + temp*d_u[i] for i in local_range[2:len]];   
        elseif right_boundary 
       	       new_u[len] = this_lambda*d_u[local_range[len-1]] + temp*d_u[local_range[len]];
       	        new_u[1:len-1] = [this_lambda*(d_u[i+1] + d_u[i-1]) + temp*d_u[i] for i in local_range[1:len-1]];
        else
		println("Error with Bounds"); 
        	#return;
        end 
	end 
        @sync begin  
            localpart(d_u)[:] = copy(new_u);
        end 
    end
  
    return localpart(d_u); 
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
    
    refs = Array(RemoteRef,nworkers());
    counter = 1;  
    for p in procs(distributed_u) 	
    	refs[counter] = @spawnat p parallel_explicit_stepper!(M,distributed_u,lambda);
	counter += 1; 
    end
    
    for k in 1:nworkers();
    	this_range = (distributed_u.indexes[k])[1]; 
    	u[this_range] = fetch(refs[k])[1];     
    end  

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
