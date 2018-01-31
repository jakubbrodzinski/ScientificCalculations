#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#
#= given function =#
function f(x::Float64)
    return Float64(sin(x)+ cos(Float64(3.0)x))
end

#= function calculates derivative of function f using approximatied equation =#
function df(x0::Float64,h::Float64)
    return Float64((f(x0+h)-f(x0))/h)
end

#= function calculates derivative of function f using exact equation =#
function exact_df(x0::Float64)
    return Float64(cos(x0)-Float64(3.0)sin(Float64(3.0)x0))
end

#= function calculates approxmiation error in our case of  calculating derivatives of function f. =#
function approximationError(df1::Float64,df2::Float64)
    return abs(df1-df2)
end

#= prints derivative of function f in x0 point =#
x0=Float64(1.0)
h=Float64(1.0)
for i = 0:54
    println("h=",h,"\ndf(",x0,")=",df(x0,h))
    h=h/Float64(2.0)
end

#= prints approximation error for diffrent h =#
h=Float64(1.0)
exactValue= exact_df(x0)
for i = 0:54
    println("h=",h,"\naproxError= ",approximationError(df(x0,h),exactValue))
    h=h/Float64(2.0)
end
