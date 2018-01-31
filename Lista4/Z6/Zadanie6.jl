workspace();
include("../interp_polyn.jl")
using InterPolynomialModule;

#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#

f1(x)=abs(x)
f2(x)=1/(1+x^2)

n=15

#rysujNnfx(f1,-1.0,1.0,n);
rysujNnfx(f2,-5.0,5.0,n);
