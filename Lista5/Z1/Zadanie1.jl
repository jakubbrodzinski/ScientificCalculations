workspace();
include("../file_handling.jl")
include("../blocksys.jl")
include("../matrixgen.jl")
using FileHandling
using BlockSys
using matrixgen
using PyPlot
#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#

#(b,n2)=loadB("./../b10k.txt");
(A,n,l)=loadA("A.txt");
#(A,n,l)=loadA("./../Z3/Atest.txt")
b=compute_b(A,n);
tic();
#comX=gauss(A,b,n,l);
comX=gaussDom(A,b,n,l);
toc();
error=computeError(comX);
saveSolutionWithError(comX,error,"result")
#println(b);
#saveSolutionWithError(comX,error,"result1")
#show(IOContext(STDOUT), "text/plain", full(A));
#println(error);
