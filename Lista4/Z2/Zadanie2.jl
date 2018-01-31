workspace();
include("../interp_polyn.jl")
using InterPolynomialModule;

#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#

#=
Two very simple test cases of ilorazyRoznicowe method.
=#
x=[3.0,1.0,5.0,6.0];
f=[1.0,-3.0, 2.0, 4.0];
t=3.0
t2=5.0
fx=ilorazyRoznicowe(x,f);
res=warNewton(x,fx,t)
res2=warNewton(x,fx,t2)
println("wartosc wielomianu w punkcie ",t,"=",res);
println("wartosc wielomianu w punkcie ",t2,"=",res2);
