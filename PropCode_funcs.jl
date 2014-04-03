const age = 13.82     #current age of the universe [Gyr] 
const H0 = 67.11      #Hubble constant at the present time [km/s/Mpc]
const H0_gyr = 0.001022689*H0; #Hubble constant at the present time [Gyr^-1]
const inverse_H0 = 1.0/H0_gyr #1/H0 [Gyr] 
const Omegam = 0.3175 #fraction of matter in the unverise taken from Planck 
const Omegavac = 1.0 - Omegam; 
const c = 299800;  #speed of light [km/s]  
const coverH0 = c/67.11;  #[Mpc]
const alpha = 9.7222e-4;  #converstion from Gyr/Mpc to s/km 
const inverse_alpha = 1.0/alpha; 

#calculates and returns the derivative of time wrt redshift 
#takes a redshift as input 
function dtdz (this_z::Float64)
    @assert(this_z >= 0.0)
    this_opz = 1+this_z; 
    temp_result = this_opz*sqrt(Omegam*this_opz*this_opz*this_opz + Omegavac); 
    return -inverse_H0/temp_result; 
end

function dzdt(this_z::Float64)
    @assert(this_z >= 0.0);
    this_opz = 1.0+this_z; 
    return -H0_gyr*this_opz*sqrt(Omegam*this_opz*this_opz*this_opz + Omegavac); 
end

function drdz (this_z::Float64)
    @assert(this_z >= 0.0)
    this_opz = 1.0+this_z; 
    return coverH0/sqrt(Omegam*this_opz*this_opz*this_opz + Omegavac);
end

function calc_dt(max_H::Float64,max_r::Float64,this_dr::Float64)
    denominator = (alpha*(c+H0*max_H*max_r)/this_dr)^2 + 9.0*(H0_gyr*max_H)^2; 
    return 3.0*H0_gyr*max_H/denominator; 
end

function calc_w(this_H::Float64,this_r::Float64,this_dt::Float64,this_dr::Float64)
    w_sqrd = (1.0 - 3.0*this_dt*H0_gyr*this_H)^2; 
    #println("term 1 is ",w_sqrd); 
    w_sqrd = w_sqrd + (alpha*this_dt*(c+H0*this_H*this_r)/this_dr)^2;
    #println("term 2 is ",w_sqrd); 
    return sqrt(w_sqrd); 
end


#TODO may want to make A a sparce matrix 
function make_A!(A_mat::Array{Float64,2},this_H::Float64,this_r::Array{Float64,1},this_dt::Float64,this_dr::Float64,source_position::Float64)
    this_size = size(A_mat); 
    num = length(this_r); 
    @assert(this_size == (num,num));
    @assert(num == this_size[1]); 
    @assert(this_H >= 1.0);
    @assert(this_dt > 0.0); 
    @assert(this_dr > 0.0); 
    
    A_mat[1,1] = (1.0-3.0*this_dt*H0_gyr*this_H + alpha*c*this_dt/this_dr); 
    A_mat[1,2] = -alpha*c*this_dt/this_dr; 
    
    for i in 2:num-1
        a =  (this_r[i] < source_position) ? -0.5*alpha*this_dt*(c+H0*this_H*this_r[i])/this_dr :  -0.5*alpha*this_dt*(H0*this_H*this_r[i]-c)/this_dr;  
        A_mat[i,i-1] = a; 
        A_mat[i,i] = 1.0-3.0*this_dt*H0_gyr*this_H;
        A_mat[i,i+1] = -a; 
    end
    
    A_mat[num,num-1] = alpha*this_dt*(c+H0*this_H*this_r[num])/this_dr; 
    A_mat[num,num] = 1.0-3.0*this_dt*H0_gyr*this_H - alpha*this_dt*(c+H0*this_H*this_r[num])/this_dr; 

end

function prop_step(old::Array{Float64,1},new::Array{Float64,1},mat::Array{Float64,2})
    @assert(length(old) == length(new)); 
    @assert(size(mat)[1] == length(new)); 
    
    return mat*old;
end

function prop_step(old::Array{Float64,1},new::Array{Float64,1},mat::Array{Float64,2},source::Array{Float64,2})
    @assert(length(old) == length(new)); 
    @assert(size(mat)[1] == length(new)); 
    @assert(length(old) == length(source)); 
    
    return mat*old + source;
end
