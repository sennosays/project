function u0_test(x::Float64)
    return sin(2pi*x)
end


function exact_solution_test(x::Float64,t::Float64)
    return exp(-4*pi^2*t)*sin(2*pi*x);
end

function explicit_solver(u0,this_xi::Float64,this_xf::Float64,this_tf::Float64,this_D::Float64,nx::Int = 20)
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
    
    temp = 1.0-2.0*lambda;
    diag = temp*ones(nx); 
    off_diag = lambda*ones(nx-1); 
    
    A = spdiagm((off_diag,diag,off_diag),[-1,0,1]); 
    new_u = Array(Float64,nx); 
    old_u = map(u0,x[2:end-1]);
    
    tic()
    for j in 0:M-1
        new_u = A*old_u; 
        old_u = new_u; 
    end
    time = toq(); 
    
    prepend!(old_u,[0.0]); 
    append!(old_u,[0.0]);
    
    println("average time for serial solver is ",time/M); 
    return old_u, x;
end

function implicit_solver(u0,xi::Float64,xf::Float64,tf::Float64,D::Float64,nx::Int = 20,nt::Int=20)
    @assert(xf>xi); 
    @assert(tf>0.0); 
    @assert(D > 0.0); 

    x = linspace(xi,xf,nx+2); 
    dx = x[2] - x[1]; 
    @assert(dx == (xf-xi)/(nx+1)); 
    dt = tf/nt; 
    lambda = dt*D/dx^2; 
    
    temp = 1.0+2.0*lambda;
    diag = temp*ones(nx); 
    off_diag = -lambda*ones(nx-1); 
    
    B = spdiagm((off_diag,diag,off_diag),[-1,0,1]); 
    B = full(B); 
    inverse_B = inv(B); 
    new_u = Array(Float64,nx); 
    old_u = map(u0,x[2:end-1]); 
    
    for j in 1:nt
        new_u = inverse_B*old_u; 
        old_u = new_u; 
    end
    prepend!(old_u,[0.0]); 
    append!(old_u,[0.0]); 
    return old_u, x; 
    
end

function crank_solver(u0,xi::Float64,xf::Float64,tf::Float64,D::Float64,nx::Int = 20,nt::Int=20)
    @assert(xf>xi); 
    @assert(tf>0.0); 
    @assert(D > 0.0); 

    x = linspace(xi,xf,nx+2); 
    dx = x[2] - x[1]; 
    @assert(dx == (xf-xi)/(nx+1)); 
    dt = tf/nt; 
    lambda = dt*D/dx^2;
    
    I = speye(nx); 
    
    temp = 1.0-2.0*lambda;
    diag = temp*ones(nx); 
    off_diag = lambda*ones(nx-1); 
    
    A = spdiagm((off_diag,diag,off_diag),[-1,0,1]); 
    temp = 1.0+2.0*lambda; 
    diag = temp*ones(nx); 
    B = spdiagm((-off_diag,diag,-off_diag),[-1,0,1]); 
    B = full(B); 
    inverse_opB = inv(I+B); 
    
    new_u = Array(Float64,nx); 
    old_u = map(u0,x[2:end-1]); 
    
    for j in 1:nt
        new_u = inverse_opB*((I+A)*old_u); 
        old_u = new_u; 
    end
    prepend!(old_u,[0.0]); 
    append!(old_u,[0.0]); 
    return old_u, x; 
end

function MC_crank_solver(u0,xi::Float64,xf::Float64,tf::Float64,D::Float64,nx::Int = 50,num_walkers::Int = 1000)
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
    counts = zeros(Int,nx+2); 
    probs = zeros(Float64,nx+2); 
    probability = (1.0-lambda)/(1.0+lambda); 
    solution[1] = 0.0;
    solution[end] = 0.0;

    for i = 2:nx+1
        boundary_hits = 0; 
        fill!(probs,0.0); 
        fill!(counts,0.0); 
        for n = 1:num_walkers
            position_index = i; 
            time_index = nt; 
            
            while((time_index > 0) && (position_index > 1) && (position_index < nx+2))
                if(rand() < probability)
                    time_index -= 1;
                else 
                    r = rand(); 
                    if r < 0.25 
                        position_index += 1; 
                    elseif r < 0.5 
                        position_index -= 1;
                    elseif r < 0.75
                        position_index += 1; 
                        time_index -= 1; 
                    else
                        position_index -= 1;
                        time_index -= 1; 
                    end
                end
            end
            if time_index > 0
                boundary_hits += 1; 
            else
                counts[position_index] += 1; 
            end
        end
        @assert(sum(counts) + boundary_hits == num_walkers);
        probs = counts./num_walkers; 
        solution[i] = sum(probs.*initial_condition); 
    end

    return solution, x; 
end

function map_MC_crank_solver(u0,xi::Float64,xf::Float64,tf::Float64,D::Float64,nx::Int = 50,num_walkers::Int = 1000)
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

    solution[2:end-1] = map(a -> sum(initial_condition.*random_walker(a,1,nx+2,nt,num_walkers,probability)),2:nx+1);

  

    return solution, x; 
end

function random_walker(initial_position::Int, left_bound::Int, right_bound::Int, initial_time::Int,n_walkers::Int,this_probability::Float64)
	 @assert(left_bound < initial_position < right_bound);
	 @assert(initial_time > 1); 
	 @assert(0.0 < this_probability < 1.0); 
	 counts = zeros(Int, right_bound-left_bound + 1); 
	 boundary_hits = 0; 

	 for n in 1:n_walkers
	    position_index = initial_position; 
            time_index = initial_time; 
            
            while((time_index > 0) && (position_index > left_bound) && (position_index < right_bound))
                if(rand() < this_probability)
                    time_index -= 1;
                else 
                    r = rand(); 
                    if r < 0.25 
                        position_index += 1; 
                    elseif r < 0.5 
                        position_index -= 1;
                    elseif r < 0.75
                        position_index += 1; 
                        time_index -= 1; 
                    else
                        position_index -= 1;
                        time_index -= 1; 
                    end
                end
            end
            if time_index > 0
                boundary_hits += 1; 
            else
                counts[position_index] += 1; 
            end
        end
	@assert(sum(counts) + boundary_hits == n_walkers);
        return counts./n_walkers; 
end
