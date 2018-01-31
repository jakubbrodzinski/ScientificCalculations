#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#

include("matcond.jl")
include("hilb.jl")

#=
    Function that calculates x from equasion Ax=b using first given algorithm (Gauss method, x=b/A=A\b). The goal is to check how much x will vary form its actual value [1,...,1].
    variables:
        A- matrix A from Ax=b
        n - matrix A is n x n
        b - variable from equasion
        x - variable from equasion, we want to calculate it
        result - calculated x with Gauss' method
        relErr - relative error that occured during calculating x from equasion (variable result)
=#
function findX_FirstAlgorithm(A::Array{Float64,2},n::Integer)
    x= fill(Float64(1.0),n);
    b=A*x;

    result=Array{Float64,1}(A\b);
    relErr=norm(result-x)/norm(x);

    println("cond= ",cond(A),"\t n=",n,"\trank=",rank(A));
    #println("x= ",result,"\t err= ",relErr)
    println("err= ",relErr)
end

#=
    Function that calculates x from equasion Ax=b using second given algorithm (inv(A)*b, inv(A) = A^(-1)). The goal is to check how much x will vary form its actual value [1,...,1].

    variables:
        A- matrix A from Ax=b
        n - matrix A is n x n
        b - variable from equasion
        x - variable from equasion, we want to calculate it
        result - calculated x with second algorithm
        relErr - relative error that occured during calculating x from equasion (variable result)
=#

function findX_SecondAlgorithm(A::Array{Float64,2}, n::Integer)
    x= fill(Float64(1.0),n);
    b=A*x

    result=Array{Float64,1}(inv(A)*b)
    relErr=norm(result-x)/norm(x);

    println("cond= ",cond(A),"\t n=",n,"\trank=",rank(A));
    #println("x= ",result,"\t err= ",relErr)
    println("err= ",relErr)
end

#=
    Firstly we display results for Hilbert's matrixes (that have diffrent sizes).
    Secondly we do the same with random n x n matrixes.
    We want to compare errors and results from both of our algorithms.
    Funtions "hilb()" and "matcond()" were linked to the list with excercises.
=#

println("Hilbert, firstAlgorithm:");
for n in [5,10,20]
    findX_FirstAlgorithm(hilb(n),n);
end
println("Hilbert, secondAlgorithm:");
for n in [5,10,20]
    findX_SecondAlgorithm(hilb(n),n);
end

println("Random matrix, firstAlgorithm:");
for n in [5,10,20]
    for c in [Float64(1.0),Float64(10.0),Float64(1.0e3),Float64(1.0e7),Float64(1.0e12),Float64(1.0e16)]
        findX_FirstAlgorithm(matcond(n,c),n);
        println("c= ",c);
    end
end
println("Random matrix, secondAlgorithm:");
for n in [5,10,20]
    for c in [Float64(1.0),Float64(10.0),Float64(1.0e3),Float64(1.0e7),Float64(1.0e12),Float64(1.0e16)]
        findX_SecondAlgorithm(matcond(n,c),n);
        println("c= ",c);
    end
end
