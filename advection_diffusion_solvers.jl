u0(this_x::Float64,this_xg::Float64 = 50.0,this_sigma::Float64=5.0) = exp(-0.5*(this_x-this_xg)^2/this_sigma^2);

function energy_spectrum(this_E::Float64,min_E::Float64,max_E::Float64)
    @assert(min_E >= 0.0); 
    @assert(max_E > min_E); 
    @assert(this_E >= 0.0); 
    
    if this_E <= 0.5*min_E
        return 0.0;
    elseif this_E <= min_E 
        return exp(-(this_E-min_E)^2)/this_E/this_E; 
    elseif this_E < max_E 
        return 1/this_E/this_E;
    else
        return exp(-(this_E-max_E)^2)/this_E/this_E;
    end
end

function calc_diffusion_constant(this_E::Float64, s::Float64 = 0.4)
    @assert(this_E >= 0.0); 
    return 5.0*minimum((1.0,this_E^(-s))); 
end
    
calc_energy_loss_rate(this_E::Float64, this_r::Float64 = 0.5) = 0.1*this_E^this_r; 

function smooth_grid!(this_grid::Array{Float64,2},tol::Float64)
    @assert(tol>0.0); 
    this_max, this_max_loc = findmax(abs(this_grid));
    counter = 0; 
    while(this_max > tol)
        this_grid[this_max_loc] = 0.0; 
        this_max, this_max_loc = findmax(abs(this_grid)); 
        counter += 1; 
    end
    println(counter," points were needed to be smoothed")
end

function calc_MC_probabilities!(this_probabilities::Array,this_alpha::Array{Float64,1}, this_beta::Array{Float64,1})
    num = length(this_alpha); 
    @assert(num == length(this_beta));
    @assert(num == length(this_probabilities[:,1]));
    @assert(length(this_probabilities[1,:]) == 3);

    normalization = 1.0+this_alpha[1] + this_beta[2]
    this_probabilities[1,:] = [this_alpha[1]/normalization, (0.5 - 0.5* this_alpha[1] + this_beta[2])/normalization, (0.5 - 0.5*this_alpha[1])/normalization];

    for k in 2:num - 1
        normalization = 1 + this_alpha[k] + this_beta[k+1] - this_beta[k-1];
        this_probabilities[k,:] = [this_alpha[k]/normalization, (0.5 - 0.5* this_alpha[k] + this_beta[k+1])/normalization, (0.5 - 0.5*this_alpha[k] - this_beta[k-1])/normalization];
    end

    normalization = 1.0 + this_alpha[num] - this_beta[num-1];
    this_probabilities[num,:] = [this_alpha[num]/normalization, (0.5 - 0.5* this_alpha[num])/normalization, (0.5 - 0.5*this_alpha[num] - this_beta[num-1])/normalization];
    this_probabilities *= 0.5;

    for k in 2:num-1
        if abs((4.0*this_probabilities[k,1] + 2.0*this_probabilities[k,2] + 2.0*this_probabilities[k,3]) - 1.0) > 0.0001
            println("the offending probability is ",4.0*this_probabilities[k,1] + 2.0*this_probabilities[k,2] + 2.0*this_probabilities[k,3]); 
            println("at index ",k);
        end
        @assert(abs((4.0*this_probabilities[k,1] + 2.0*this_probabilities[k,2] + 2.0*this_probabilities[k,3]) - 1.0) < 0.0001)
    end
end

function MC_walker(this_i::Int, this_k:: Int,this_probabilities::Array,this_initial_condition::Array,this_nt::Int, this_num_walkers::Int = 1000)
    this_nx = length(this_initial_condition[:,1]);
    this_nE = length(this_initial_condition[1,:]);
    @assert(this_i < this_nx);
    @assert(this_k < this_nE);
    counts = zeros(Int,this_nx,this_nE);
    probs = zeros(Float64,this_nx,this_nE);
    boundary_hits = 0;


    for n = 1:this_num_walkers
        position_index = this_i;
        energy_index = this_k;
        time_index = this_nt;
            
        while((time_index > 0) && (this_nx-1 > position_index > 2) && (this_nE-1 > energy_index > 2))
            r1 = rand()
            if(r1 < 4.0*this_probabilities[this_k,1])
                r2 = rand();
                if r2 < 0.25
                    position_index += 1;
                elseif r2 < 0.5
                    position_index += 1;
                    time_index -= 1;
                elseif r2 < 0.75
                    position_index -= 1;
                else
                    position_index -= 1;
                    time_index -= 1;
                end
            elseif(r1 < (4.0*this_probabilities[this_k,1] + 2.0*this_probabilities[this_k,2]))
                r2 = rand();
                energy_index += 1;
                if(r2 < 0.5)
                    time_index -= 1;
                end
            else
                r2 = rand();
                energy_index -= 1;
                if(r2 < 0.5);
                    time_index -= 1;
                end
            end
        end
        if time_index > 0
            boundary_hits += 1;
        else
            counts[position_index,energy_index] += 1;
        end
    end
    @assert(sum(counts) + boundary_hits == this_num_walkers);
    probs = counts./this_num_walkers;
            
            
    return sum(probs.*this_initial_condition);
end

function MC_advection_diffusion_solver(xi::Float64,xf::Float64,nx::Int,Ei::Float64,Ef::Float64,nE::Int,spec_E_min::Float64,spec_E_max::Float64,tf::Float64,num_walkers::Int = 1000)
    srand(1620);
    @assert(xf > xi); 
    @assert(nx > 3); 
    @assert(Ef > Ei); 
    @assert(nE > 3); 
    @assert(spec_E_max > spec_E_min)
    @assert(Ef > spec_E_max); 
    @assert(spec_E_min > Ei); 
    @assert(tf > 0.0); 
    
    dx = (xf-xi)/(nx-1);
    x = [xi+dx*i for i in 0:nx-1];

    dE = (Ef - Ei)/(nE-1); 
    E = [Ei + dE*k for k in 0:nE-1];

    alpha = map(calc_diffusion_constant,E);
    beta = -map(calc_energy_loss_rate,E);

    nt = iceil(tf*maximum(alpha./dx^2 + beta/dE));
    dt = tf/nt;

    alpha *= (dt/dx^2);
    beta *= (dt*0.5/dE);

    @assert(minimum(0.5.-0.5*alpha - beta) > 0)

    initial_condition = [energy_spectrum(a,spec_E_min,spec_E_max)*u0(b,0.5*(xf+xi),1.0) for b=x, a=E];
    solution = zeros(Float64,nx,nE);
    
    probabilities = Array(Float64,nE,3); 
    calc_MC_probabilities!(probabilities,alpha, beta);

    for i = 2:nx-1
        if( i%(0.1*nx) == 0)
            println(100*i/nx,"% completed");
        end
        for k = 2:nE-1
            solution[i,k] = MC_walker(i,k,probabilities,initial_condition, nt, num_walkers)
        end
    end
    
    #solution[1,:] = 0.0;
    #solution[end,:] = 0.0;
    #solution[:,1] = 0.0;
    #solution[:,end] = 0.0;
    println(100.0,"% completed");
    
    return solution, x, E, initial_condition;
    

end

function exp_advection_diffusion_solver(xi::Float64,xf::Float64,nx::Int,Ei::Float64,Ef::Float64,nE::Int,spec_E_min::Float64,spec_E_max::Float64,tf::Float64,nt::Int)
    @assert(xf > xi); 
    @assert(nx > 3); 
    @assert(Ef > Ei); 
    @assert(nE > 3); 
    @assert(spec_E_max > spec_E_min)
    @assert(Ef > spec_E_max); 
    @assert(spec_E_min > Ei); 
    @assert(tf > 0.0); 
    
    dx = (xf-xi)/(nx-1);
    x = [xi+dx*i for i in 0:nx-1];

    dE = (Ef - Ei)/(nE-1); 
    E = [Ei + dE*k for k in 0:nE-1];

    alpha = map(calc_diffusion_constant,E);
    beta = map(calc_energy_loss_rate,E);

    dt = tf/nt;

    alpha *= (dt/dx^2);
    beta *= (dt*0.5/dE);

    old_u = [energy_spectrum(a,spec_E_min,spec_E_max)*u0(b,0.5*(xf+xi),1.0) for b=x, a=E];
    original_u = copy(old_u); 
    new_u = Array(Float64,nx,nE);
    
    alpha = map(calc_diffusion_constant,E); 
    alpha *= (dt/dx^2); 
    beta = map(calc_energy_loss_rate,E); 
    beta *= (dt*0.5/dE); 
    
    w = (1.0 .- 2.0*alpha).^2 + 4*beta.^2; 
    println("the maximum of w^2 is ",maximum(w)); 
    @assert(maximum(w) < 1.0); 
    
    for j = 0:nt
        [new_u[i,k] = alpha[k]*(old_u[i+1,k]+old_u[i-1,k]) + (1.0-2.0*alpha[k])*old_u[i,k] + beta[k+1]*old_u[i,k+1] - 
            beta[k-1]*old_u[i,k-1] for k = 2:nE-1, i = 2:nx-1];
        old_u = new_u;
        old_u[1,:] = 0.0; 
        old_u[end,:] = 0.0;
    end 
    
    return old_u, x, E, original_u;
end

function get_x_index(min_x::Float64, max_x::Float64, x::Float64, num_x::Int)
    @assert(max_x > x > min_x); 
    dx = (max_x - min_x)/(num_x-1); 

    return iround((x-min_x)/dx);
end

function get_E_index(min_E::Float64, max_E::Float64, E::Float64, num_E::Int)
    @assert(max_E > E > min_E);
    dE = (max_E - min_E)/(num_E-1); 
    
    return iround((E - min_E)/dE); 
end