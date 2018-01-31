workspace();
include("../file_handling.jl")
include("../blocksys.jl")
include("../matrixgen.jl")
using FileHandling
using BlockSys
using matrixgen
#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#

(A,n,l)=loadA("./../A50k.txt");
#(b,n)=loadB("./../b50k.txt");
#blockmat(1000, 4, 10.0, "Atest.txt")
#(A,n,l)=loadA("AA.txt");
#show(IOContext(STDOUT), "text/plain", full(A)); println();
b=compute_b(A,n);
#println();
tic();
(L,U,b2)=LUmatrixesDom(A,b,n,l);
x=findXWithLU(L,U,b2,n,l);
toc();
#show(IOContext(STDOUT), "text/plain", full(L*U)); println();
#println(b2);
#saveSolutionWithError(b2,0.0,"arek")
#(L,U)=LUmatrixes(A,n,l);
#show(IOContext(STDOUT), "text/plain", full(L));
#println();
#show(IOContext(STDOUT), "text/plain", full(U));
#println();
#x=findXWithLU(L,U,b2,n,l);
error=computeError(x);
saveSolutionWithError(x,error,"result")
println(error);
