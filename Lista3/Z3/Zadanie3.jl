workspace();

include("../f_roots.jl")
using RootsModule;

#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#

#=
    Program created to test Implementation of the secant method that is used to find function's roots.
    Before executing this program recommended is to look into 'm_bis.jl' file.
    In this case: f(x)=x^3-20 and g(x)=x^2-16
=#

f(x)=x^3-20;
g(x)=x^2-16;
δ=-2e-5
ϵ=0.5e-5;
println(msiecznych(f,2.0,3.0,δ,ϵ,10));  # it should be (0,0,0,2)
δ=0.5e-5;
println(msiecznych(f,0.0,1.0,δ,ϵ,2)); # it should be (0,0,0,1)
println(msiecznych(f,2.0,3.0,δ,ϵ,10)); # it should return result that f(x) < ϵ
println(msiecznych(g,3.5,4.5,δ,ϵ,10)); # it should return actual result
