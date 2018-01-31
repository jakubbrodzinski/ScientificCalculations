workspace();
include("../file_handling.jl")
include("./lu_module.jl")
using FileHandling
using LUGauss
#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#

function computeError(b_file::Vector{Float64},b_computed::Vector{Float64})
    return norm(b_computed-b_file)/norm(b_file);
end

(A,n,l)=loadA("./../A50k.txt");
(b,n2)=loadB("./../b50k.txt");

tic();
(L,U,b)=LUmatrixesDom(A,b,n,l);
(L,U)=LUmatrixes(A,n,l);
toc();
#saveSolution(b,"result")

b_computed=compute_b(A,n,l);
error=computeError(b,b_computed);
println();
#show(IOContext(STDOUT), "text/plain", full(L));
println();
#show(IOContext(STDOUT), "text/plain", full(U));
println();
