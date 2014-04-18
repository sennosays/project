@everywhere function parallel_explicit_stepper!(num_steps::Int,d_u,local_u::Array{Float64,1},this_lambda::Float64)
    println("this is the newest explicit stepper"); 
    this_range = (d_u.indexes[myid()-1])[1];
    len = length(this_range);
    @assert(len == length(local_u)); 
    new_u = Array(Float64,len); 
    
    left_boundary = false; 
    right_boundary = false; 
    middle = false; 
    if this_range[1] == 1;
        left_boundary = true; 
    elseif this_range[len] == length(d_u); 
        right_boundary = true; 
	left_point_index = length((d_u.indexes[myid()-2])[1]); 
    elseif !left_boudary & !right_boundary 
        midde = true;
        left_point_index = length((d_u.indexes[myid()-2])[1]); 
    else
        println("Error determining place in array"); 
        return 
    end
    
    temp = 1.0-2.0*this_lambda; 

    for j in 1:num_steps
        if middle 
            left_point = remotecall_fetch(myid()-1,getindex,localpart(d_u),left_point_index); 
            right_point = remotecall_fetch(myid()+1,getindex,localpart(d_u),1); 
            new_u[1] = this_lambda*(local_u[2]+left_point) + temp*local_u[1]; 
            new_u[len] = this_lambda*(local_u[len-1]+right_point) + temp*local_[len]; 
        elseif left_boundary 
            right_point = remotecall_fetch(myid()+1,getindex,localpart(d_u),1); 
            new_u[1] = this_lambda*local_u[2] + temp*local_u[1]; 
            new_u[len] = this_lambda*(local_u[len-1] + right_point) + temp*local_u[len]; 
        elseif right_boundary
	    r = @spawnat myid()-1 getindex(localpart(d_u),left_point_index); 
	    left_point = fetch(r); 
            new_u[1] = this_lambda*(local_u[2] + left_point) + temp*local_u[1]; 
            new_u[len] = this_lambda*local_u[len-1] + temp*local_u[len];
        else
            println("Error with Bounds"); 
            return
        end
        new_u[2:len-1] = [this_lambda*(local_u[i-1] + local_u[i+1]) + temp*local_u[i] for i in 2:len-1];
        
        for k in 1:len 
            local_u[k] = new_u[k]; 
        end
    end
    local_u;   
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
    #for j in 0:M-1
    #    new_u[1] = lambda*old_u[2] + temp*old_u[1]; 
    #    for i in 2:nx-1
    #        new_u[i] = lambda*(old_u[i+1]+old_u[i-1]) + temp*old_u[i]; 
    #    end
    #    new_u[nx] = lambda*old_u[nx-1] + temp*old_u[nx];
    #    old_u = copy(new_u); 
    #end 
    refs = Array(RemoteRef,2); 
    counter = 1; 
    for p in procs(distributed_u) 
        refs[counter] = @spawnat p parallel_explicit_stepper!(1000,distributed_u,localpart(distributed_u),lambda);
	counter += 1; 
    end

    println(fetch(refs[2])); 


    #plot(x[2:end-1],distributed_u); 

    #u = convert(Array{Float64,1},distributed_u); 

    prepend!(u,[0.0]); 
    append!(u,[0.0]); 

    return u, x;
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
