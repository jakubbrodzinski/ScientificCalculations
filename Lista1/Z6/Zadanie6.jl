#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#

#= function calculate f(x) and g(x) given in exercise =#

function f(x::Float64)
    return Float64(sqrt(x^2+1.0)-1.0)
end

function g(x::Float64)
    return Float64(x^2/Float64(x^2+1)+1)
end

#= presenting calculated values for diffrent arguemnts. it shows that results vary =#

x=8e-1
η=nextfloat(Float64(0.0))
while(x>η)
    println("f(",x,")=",f(x),"\tg(",x,")=",g(x))
    x=x/8.0
end
