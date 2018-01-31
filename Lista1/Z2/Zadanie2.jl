#=
    Jakub Brodziński
    229781@student.pwr.edu.pl
=#

#=
    Functions calculates in mechanical ϵ for Flot16/Float32/Float64 using Kahan's method
    Mechanical ϵ - the smallest number that added to 1.0 gives more than 1.0.
=#

function kahanEpsFloat16()
    return Float16(3.0)*(Float16(4.0)/Float16(3.0)-Float16(1.0))-Float16(1.0)
end

function kahanEpsFloat32()
    return Float32(3.0)*(Float32(4.0)/Float32(3.0)-Float32(1.0))-Float32(1.0)
end

function kahanEpsFloat64()
    return Float64(3.0)*(Float64(4.0)/Float64(3.0)-Float64(1.0))-Float64(1.0)
end

#=
    Values calculated by earlier implemented functions are compared to the (proper) Julia's values
=#

println("kahanEps (realEps)")
println(kahanEpsFloat16()," (",eps(Float16),")")
println(kahanEpsFloat32()," (",eps(Float32),")")
println(kahanEpsFloat64()," (",eps(Float64),")")
