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
println(ilorazyRoznicowe([1.0,2.0],[1.0,2.0]));
println(ilorazyRoznicowe([3.0, 1.0, 5.0, 6.0], [1.0, -3.0, 2.0, 4.0]));
