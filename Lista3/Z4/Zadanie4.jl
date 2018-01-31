workspace();

include("../f_roots.jl")
using RootsModule;

#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#

#=
    Program created to find roots of function f(x)=sin(x)-(x/2.0)^2 using three diffrent methods : biscetion method, Newton method, Scatent method.
    f(x),δ,ϵ were given.
=#
f(x)= sin(x)-(x/2.0)^2;
df(x)= cos(x)-x/2.0;
δ=0.5e-5;
ϵ=0.5e-5;
println("Bisekcja:\n",mbisekcji(f,1.5,2.0,δ,ϵ));
println("Newton:\n",mstycznych(f,df,1.5,δ,ϵ,50))
println("Siecznych:\n",msiecznych(f,1.0,2.0,δ,ϵ,50))
