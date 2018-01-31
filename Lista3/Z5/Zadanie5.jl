workspace();

include("../f_roots.jl")
using RootsModule;

#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#

#=
    Program created to find common point of two functions y=3x and y=e^x with biscetion method.
    δ,ϵ were given.
=#

f(x)=3x-e^x;
δ=1e-4;
ϵ=1e-4;
println("Bisekcja:\n",mbisekcji(f,-0.5,1.0,δ,ϵ),"\n",mbisekcji(f,1.0,3.0,δ,ϵ));
