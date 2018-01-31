#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#

#=
    Functions calculates in while loop mechanical ϵ for Flot16/Float32/Float64
    Mechanical ϵ - the smallest number that added to 1.0 gives more than 1.0.

    variables:
        myMachEps - calculated value
=#

function myMachEpsFloat16()
    myMachEps= Float16(1.0)
    while Float16(1.0)+myMachEps/Float16(2.0)>Float16(1.0)
        myMachEps=myMachEps/Float16(2.0)
    end
    return myMachEps
end

function myMachEpsFloat32()
    myMachEps= Float32(1.0)
    while Float32(1.0)+myMachEps/Float32(2.0)>Float32(1.0)
        myMachEps=myMachEps/Float32(2.0)
    end
    return myMachEps
end

function myMachEpsFloat64()
    myMachEps= Float64(1.0)
    while Float64(1.0)+myMachEps/Float64(2.0)>Float64(1.0)
        myMachEps=myMachEps/Float64(2.0)
    end
    return myMachEps
end

#=
    Functions calculates in while loop η for Flot16/Float32/Float64
    η - the smallest number that is bigger than 0.0

    variables:
        myEta - calculated value
=#

function myEtaFloat16()
    myEta=Float16(1.0)
    while myEta/Float16(2.0)>Float16(0.0)
            myEta=myEta/Float16(2.0)
    end
    return myEta
end

function myEtaFloat32()
    myEta=Float32(1.0)
    while myEta/Float32(2.0)>Float32(0.0)
            myEta=myEta/Float32(2.0)
    end
    return myEta
end

function myEtaFloat64()
    myEta=Float64(1.0)
    while myEta/Float64(2.0)>Float64(0.0)
            myEta=myEta/Float64(2.0)
    end
    return myEta
end

#=
    Functions calculates in while loop maximum value that is not ∞ for Flot16/Float32/Float64

    variables:
        myMax - calculated value
=#

function myMaxFloat16()
    myMax= prevfloat(Float16(1.0))
    while !isinf(myMax*Float16(2.0))
        myMax=myMax*Float16(2.0)
    end
    return myMax
end

function myMaxFloat32()
    myMax=prevfloat(Float32(1.0))
    while !isinf(myMax*Float32(2.0))
        myMax=myMax*Float32(2.0)
    end
    return myMax
end

function myMaxFloat64()
    myMax=prevfloat(Float64(1.0))
    while !isinf(myMax*Float64(2.0))
        myMax=myMax*Float64(2.0)
    end
    return myMax
end

#=
    Values calculated by earlier implemented functions are compared to the (proper) Julia's values
=#

println("type\tmyMachEps\t\t(Julia's machEps)")
println("Float16\t",myMachEpsFloat16(),"\t\t(",eps(Float16),")")
println("Float32\t",myMachEpsFloat32(),"\t\t(",eps(Float32),")")
println("Float64\t",myMachEpsFloat64(),"\t\t(",eps(Float64),")")
println()
println("type\tmyEta\t\t(Julia's eta)")
println("Float16\t",myEtaFloat16(),"\t\t(",nextfloat(Float16(0.0)),")")
println("Float32\t",myEtaFloat32(),"\t\t(",nextfloat(Float32(0.0)),")")
println("Float64\t",myEtaFloat64(),"\t(",nextfloat(Float64(0.0)),")")
println()
println("type\tmyMax\t\t(Julia's max)")
println("Float16\t",myMaxFloat16(),"\t\t(",realmax(Float16),")")
println("Float32\t",myMaxFloat32(),"\t(",realmax(Float32),")")
println("Float64\t",myMaxFloat64(),"\t(",realmax(Float64),")")
