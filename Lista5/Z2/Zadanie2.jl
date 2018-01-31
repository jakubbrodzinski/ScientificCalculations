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
(A,n,l)=loadA("./../A50k.txt");
#(b,n2)=loadB("./../b16.txt");
#println(b);
b=compute_b(A,n)
#println(compute_b(A,n));
#show(IOContext(STDOUT), "text/plain", full(A));
#println();
tic();
#(L,U,b)=LUmatrixesDom(A,b,n,l);
(L,U)=LUmatrixes(A,n,l);
toc();
#saveSolution(b,"result")

#show(IOContext(STDOUT), "text/plain", full(L*U));
#println();
#show(IOContext(STDOUT), "text/plain", full(L));
#println();
#show(IOContext(STDOUT), "text/plain", full(U));
#println();
