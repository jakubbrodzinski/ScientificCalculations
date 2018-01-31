workspace();
include("../file_handling.jl")
include("./gauss_module.jl")
using FileHandling
using GaussModule
#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#

function computeError(b_file::Vector{Float64},b_computed::Vector{Float64})
    return norm(b_computed-b_file)/norm(b_file);
end

(A,n,l)=loadA("./../A50k.txt");
(b,n2)=loadB("./../b50k.txt");
b_computed=compute_b(A,n,l);
error=computeError(b,b_computed);
tic();
comB=gaussDom(A,b,n,l);
toc();
saveSolution(comB,"result")
#show(IOContext(STDOUT), "text/plain", full(A));
println();
