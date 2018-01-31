workspace();
include("../interp_polyn.jl")
using InterPolynomialModule;

#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#

#=
Simple test case of rysujNnfx
=#
f(x)=12+log(x)*3-x^2
rysujNnfx(f,0.5,15.0,2);
