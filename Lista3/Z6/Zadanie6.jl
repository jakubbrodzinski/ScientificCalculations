workspace();

include("../f_roots.jl")
using RootsModule;

#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#

#=
    Program created to find roots of function of functions f1(x)=e^(1-x)-1 and f2(x)=x*e^(-x) using three diffrent methods : biscetion method, Newton method, Scatent method.
    f_1(x),f_2(x),δ,ϵ were given.
=#
f1(x)=e^(1-x)-1
f2(x)=x*e^(-x)
δ=0.5e-5;
ϵ=0.5e-5;

df1(x)=-e^(1-x);
df2(x)=e^(-x)-x*e^(-x);

println("Bisekcja:\n","f1:\t",mbisekcji(f1,0.0,5.0,δ,ϵ),"\nf2:\t",mbisekcji(f2,-1.0,2.0,δ,ϵ));
println("Newton:\n","f1:\t",mstycznych(f1,df1,-1.0,δ,ϵ,50),"\nf2:\t",mstycznych(f2,df2,-2.0,δ,ϵ,50))
println("Siecznych:\n","f1:\t",msiecznych(f1,0.0,-1.0,δ,ϵ,50),"\nf2:\t",msiecznych(f2,-2.0,-1.0,δ,ϵ,50))
