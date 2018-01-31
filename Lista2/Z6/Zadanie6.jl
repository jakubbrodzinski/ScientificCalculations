#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#
using PyPlot;
#=
    Equation : x_(n+1)= (x_n)^2 + c
    We experiment with this sequence using diffrent x_0 and constant c.
    We use PyPlot to present result on graphs.
=#

#=
    Function returns result of the sequence x_n.
    Variables:
        x_0 - value from that we start
        c - constant from equation
        n - number of iterations
        outputArray - outputArray[i] =x_{i-1}, array with values x_{i}, that is needed to cosntruct graph
=#
function x_n(x_0::Float64, c::Float64, n::Integer)
    outputArray=Array{Float64,1}(n+1);
    outputArray[1]=x_0;
    for i in 1:n
        outputArray[i+1]=outputArray[i]^2+c;
    end
    return outputArray;
end

#= Input data and number of iteration that we r looking for. =#
x_0_array = [ 1.0, 2.0, 1.99999999999999, 1.0, -1.0, 0.75, 0.25 ];
c_array = [ -2.0, -2.0, -2.0, -1.0, -1.0, -1.0, -1.0 ];
n=40;

for i in 1:7
    output=x_n(x_0_array[i],c_array[i],n)
    println("c= ",c_array[i],"\tx_0= ",x_0_array[i],"\tresult= ",output[n+1]);
    #plot(output)
end
