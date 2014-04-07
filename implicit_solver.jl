include("PropCode_funcs_test2.jl")

function calc_dt(this_dx::Float64,this_H::Float64,this_a::Float64)
    @assert(this_dx > 0.0); 
    @assert(this_H >= 1.0); 
    @assert(this_a <= 1.0 && this_a > 0.0); 
    
    denominator = 18*(H0_gyr*this_H)^2 + 2*(alpha*c/(this_a*this_dx))^2; 
    return 6.0*H0_gyr*this_H/denominator; 
end

function calc_w(this_H::Float64,this_lambda::Float64,this_dt::Float64)
    @assert(this_H >= 1.0); 
    @assert(this_dt > 0.0); 
    
    return (1.0-3.0*this_dt*H0_gyr*this_H)^2 + 4*this_lambda^2; 
end

function make_B!(this_B::Array{Float64},this_H::Float64,this_dt::Float64,this_lambda::Float64,this_si::Int)
    this_size = size(this_B); 
    @assert(this_size[1] == this_size[2]); 
    num = this_size[1]; 
    @assert(this_dt > 0.0); 
    @assert(this_H >= 1.0); 
    
    
    this_B[1,1] = 1.0 - 3.0*this_lambda+3.0*this_dt*H0_gyr*this_H; 
    this_B[1,2] = 2*this_lambda; 
    this_B[1,3] = -this_lambda; 
    for i in 2:num-1
        this_B[i,i-1] = (i < this_si) ? -this_lambda : this_lambda; 
        this_B[i,i] = (1.0 + 3.0*this_dt*H0_gyr*this_H); 
        this_B[i,i+1] = (i < this_si) ? this_lambda : -this_lambda; 
    end
    this_B[num,num-2] = -this_lambda; 
    this_B[num,num-1] = 2*this_lambda; 
    this_B[num,num] = 1.0 - 3.0*this_lambda+3.0*this_dt*H0_gyr*this_H;
end

z_g = 0.1;  #source redshift location 
Q0 = 50.0;   #source intensity [??] 

#comoving distance to the source and maximum distance considered[Mpc] 
x_g = quadgk(a-> (1.0+a)*dtdz(a),0.0,z_g)[1];
x_max = x_g+quadgk(a->(1.0+a)*dtdz(a),z_g,2.0*z_g)[1]; 
x_g = -1c*x_g/H0;
x_max = -c*x_max/H0;
x_max = 2.1*x_g; 

#cosmological time of the universe at begining of calculation [Gyr] 
t_min = age + quadgk(dtdz,0.0,z_g)[1]/H0_gyr;

dx = 0.1; 
H_max = sqrt(Omegam*(1.0+z_g)^3+Omegavac);
a_max = 1.0/(1.0+z_g); 

#IMPORTANT added this line for the implicit solver 
dt = 0.001; 

nt = ifloor((age-t_min)/dt); 
dt = (age-t_min)/nt;

nx = iceil(x_max/dx);
dx = x_max/nx; 
println("the number of x points is ",nx); 
println("this requires ",nx^2*sizeof(Float64)/1e9," Gbytes"); 
x = linspace(0.0,x_max,nx); 
t = linspace(t_min,age,nt);
H = Array(Float64,nt); 
a = Array(Float64,nt); 
n = Array(Float64,nx); 
Q = zeros(Float64,nx);

z_old = 0.0; 
H[end] = 1.0; 
a[end] = 1.0; 

for i in 1:nt-1
    z_new = z_old - dzdt(z_old)*dt; 
    H[end-i] = sqrt(Omegam*(1.0+z_new).^3+Omegavac); 
    a[end-i] = 1.0/(1.0+z_new); 
    z_old = z_new; 
end

@assert(abs(H[1] - H_max) < 0.01)
@assert(abs(a[1] - a_max) < 0.01)


#use implicit solver 
B = Array(Float64,nx,nx); 
new_n = Array(Float64,nx); 

#initialize n and Q 
n = zeros(Float64,nx); 

#get the index that locates the point closes to the source position 
ceil_err = abs(x[iceil(x_g/dx)] - x_g); 
floor_err = abs(x[ifloor(x_g/dx)] - x_g); 
si = (ceil_err < floor_err) ? iceil(x_g/dx) : ifloor(x_g/dx); 

Q[si] = Q0/dx^3;

for j in 1:nt-1
    make_B!(B,H[j+1],dt,-0.5*alpha*c*dt/(a[j+1]*dx),si); 
    new_n = inv(B)*(n+dt*Q./a[j+1]^3); 
    n = copy(new_n); 
end

