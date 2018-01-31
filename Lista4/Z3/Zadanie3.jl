workspace();
include("../interp_polyn.jl")
using InterPolynomialModule;

#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#

#=
Two very simple test cases of naturalna method.
=#
#println(naturalna([-1.0,0.0,1.0,2.0],[3.0,-7.0,8.0,-6.0]))
fx=ilorazyRoznicowe([-1.0,0.0,1.0,2.0],[-1.0,0.0,-1.0,2.0])
println(fx);
println(naturalna([-1.0,0.0,1.0,2.0],fx))
