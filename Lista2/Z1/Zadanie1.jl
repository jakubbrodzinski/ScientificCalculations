#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#

x = [2.718281828, -3.141592654, 1.414213562, 0.577215664, 0.301029995]
y = [1486.2497, 878366.9879, -22.37492, 4773714.647, 0.000185049]

#= functions that calculate dot product of vectors that are stored in x,y arrays using 4 diffrent methods:
    -sumA sum dot product with normal order
    -sumB sum dot product from the last vector to the first one
    -sumC firstly sum dot products that are negative (from the smallest to the biggest) and then sum dot products that are positive (from biggest to smallest)
        then add those two values
    -sumD sum dot product after sorting them
=#

#= variables:
    S - acumulator
=#

function sumA(T)
    S=T(0.0)
    for i = 1:length(x)
        S=T(S+T(T(x[i])*T(y[i])))
    end
    return S
end

#=
    S-acumulator
=#

function sumB(T)
    S=T(0.0)
    for i = length(x):-1:1
        S=T(S+T(T(x[i])*T(y[i])))
    end
    return S
end

#=
    SubArray - array where dot products are stored (after its calculated)
    S_minus -acumulator for negative values
    S_plus - acumulator for positve values
    i - iterator
    S -acumulator
=#

function sumC(T)
    SubArray=T[]
    for i = 1:length(x)
        push!(SubArray,T(T(x[i])*T(y[i])))
    end
    sort!(SubArray)
    S_minus=T(0.0)
    i=1
    while i<=length(SubArray) && SubArray[i]<0.0
        S_minus=T(S_minus+SubArray[i])
        i+=1
    end
    S_plus=T(0.0)
    i=length(SubArray)
    while i>=1 && SubArray[i]>0.0
        S_plus=T(S_plus+SubArray[i])
        i-=1
    end
    return T(S_minus+S_plus)
end

#=
    SubArray - array where dot products are stored (after its calculated)
    i - iterator
    S -acumulator
=#

function sumD(T)
    SubArray=T[]
    for i = 1:length(x)
        push!(SubArray,T(T(x[i])*T(y[i])))
    end
    sort!(SubArray)
    S=T(0.0)
    for i=1:length(SubArray)
        S=T(S+SubArray[i])
    end
    return T(S)
end

#= Preseting results that we got =#
println("Proper result: ",Float64(-1.0065710700000010*(10.0^-11)) )

println("a) Float32 Result= ",sumA(Float32))
println("a) Float64 Result= ",sumA(Float64))

println("b) Float32 Result= ",sumB(Float32))
println("b) Float64 Result= ",sumB(Float64))

println("c) Float32 Result= ",sumC(Float32))
println("c) Float64 Result= ",sumC(Float64))

println("d) Float32 Result= ",sumD(Float32))
println("d) Float64 Result= ",sumD(Float64))
