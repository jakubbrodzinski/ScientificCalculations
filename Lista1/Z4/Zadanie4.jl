#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#

#=
    Function returns the smallest number bigger than 1 and smaller than 2 that fullfills the equation
    variables
        -x the value we are looking for
=#

function exercise4()
    x=nextfloat(1.0)
    while(x*(1.0/x) == 1.0 && x<2.0)
        x=nextfloat(x)
    end
    return x
end

#= prints value calculated by the exercise4 function in decimal and binary =#

println(exercise4())
println(bits(exercise4()))
