#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#

#=
    Function that calculate exact value of p_n in type of data given by user as first function's parameter.
    Variables:
        p_0 - first value of sequence p_n
        r - variable from p_n's equeaston
        n - number of iteration
        result - variable where we store p_(n-1) in every iteration
=#

function p(T,p_0::Float64,r::Float64,n::Integer)
    if(T != Float64)
        p_0=T(p_0);
        r=T(r);
    end
    result=T(p_0);
    for i in 1:n
        result=T(result+r*result*(T(1.0)-result));
    end
    return result;
end

#=
    Function that calculate slightly modified value of p_n.
    We were supposed to use floor function after 10th iteration.
    This function only works for Float32 since we don't needed to calculate modified p_n for diffrent type of data.
    Variables:
        p_0 - first value of sequence p_n
        r - variable from p_n's equeaston
        n - number of iteration
        result - variable where we store p_(n-1) in every iteration
=#

function p_modified(p_0::Float32,r::Float32,n::Integer)
    result=p_0;
    for i in 1:10
        result=Float32(result+r*result*Float32(1.0-result));
    end
    result=trunc(result,3);
    for i in 11:n
        result=Float32(result+r*result*Float32(1.0-result));
    end
    return result;
end

#=
    Firstly we display constants from p_n such as p_0, number of iteration or 'r',
    Then we show much huge diffrence we get even only slightly modifing p_n.
=#

println("1. p_0=0.01, r=3, n=40 (Float32)\np_n:");
println("result = ",p(Float32,0.01,3.0,40));
println("modified p_n:\nresult = ",p_modified(Float32(0.01),Float32(3.0),40));

#=
    We display diffrence between Flaot32 and Float64
=#

println("2. p_0=0.01, r=3, n=40\np_n:");
println("(Float32) result = ",p(Float32,0.01,3.0,40));
println("(Float64) result = ",p(Float64,0.01,3.0,40));
