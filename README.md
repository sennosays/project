project
=======

Cosmic Ray Diffusion Project

This project examines the effect of cosmic ray diffusion and energy loss on the observed spectrum of cosmic rays (sort of).
For those uncomfortable with PDE solvers there are some simple examples in the Diffusion Examples notebook. The main body
of the code is on display in the notebook Advection Diffusion Test. All of the functions are in their own .jl files as explained in the notebooks. 

Note that to run mc_test.jl, run_npoints_dependence.jl, and tests.jl one should include the locations for parallel_diffusion_solvers.jl and diffusion_solvers.jl in the command line. For example: 

julia -p 3 -L diffusion_solvers.jl -L parallel_diffusion_solvers.jl run_npoints_dependence.jl 

Please send all comments or concerns to

nbs5044@psu.edu 
