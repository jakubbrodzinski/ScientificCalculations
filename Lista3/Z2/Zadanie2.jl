workspace();

include("../f_roots.jl")
using RootsModule;

#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#

#=
    Program created to test Implementation of the Newton's method that is used to find function's roots.
    Before executing this program recommended is to look into 'm_bis.jl' file.
    In this case: f(x)=x^3-20 and g(x)=x^2-16
=#

f(x)=x^3-20;
df(x)=3x^2;
g(x)=x^2-16;
dg(x)=2x;
δ=-2e-5
ϵ=0.5e-5;
println(mstycznych(f,df,0.0,δ,ϵ,10));  # it should be (0,0,0,3)
δ=0.5e-5;
println(mstycznych(f,df,0.0,δ,ϵ,10)); # it should be (0,0,0,2)
println(mstycznych(f,df,2.0,δ,ϵ,2)); # it should be (x,x,x,1)
println(mstycznych(f,df,2.0,δ,ϵ,40)); # it should return actual result
println(mstycznych(g,dg,3.5,δ,ϵ,40)); # it should return actual result
